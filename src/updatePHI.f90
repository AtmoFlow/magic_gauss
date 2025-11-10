module updatePhi_mod
   !
   ! This module handles the time advance of the phase field phi.
   ! It contains the computation of the implicit terms and the linear
   ! solves.
   !

   use omp_lib
   use precision_mod
   use truncation, only: n_r_max, lm_max, l_max
   use radial_data, only: n_r_icb, n_r_cmb, nRstart, nRstop
   use radial_functions, only: or2, rscheme_oc, r, or1
   use num_param, only: dct_counter, solve_counter
   use physical_parameters, only: pr, phaseDiffFac, stef, ktopphi, kbotphi
   use init_fields, only: phi_top, phi_bot
   use blocking, only: lo_map, lo_sub_map, llm, ulm, st_map
   use logic, only: l_full_sphere
   use parallel_mod, only: rank, chunksize, n_procs, get_openmp_blocks
   use radial_der, only: get_ddr, get_ddr_ghost, exch_ghosts, bulk_to_ghost
   use constants, only: zero, one, two
   use fields, only: work_LMloc
   use mem_alloc, only: bytes_allocated
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use dense_matrices
   use real_matrices
   use band_matrices
   use parallel_solvers, only: type_tri_par

   implicit none

   private

   !-- Local variables
   real(cp), allocatable :: rhs1(:,:,:)
   integer :: maxThreads
   class(type_realmat), pointer :: phiMat(:)
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: phiMat_fac(:,:)
#endif
   logical, public, allocatable :: lPhimat(:)
   type(type_tri_par), public :: phiMat_FD
   complex(cp), public, allocatable :: phi_ghost(:,:)

   public :: initialize_updatePhi, finalize_updatePhi, updatePhi, assemble_phase,  &
   &         get_phase_rhs_imp

contains

   subroutine initialize_updatePhi

      integer :: ll, n_bands
      integer, pointer :: nLMBs2(:)

         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

            allocate( type_densemat :: phiMat(nLMBs2(1+rank)) )

            do ll=1,nLMBs2(1+rank)
               call phiMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
            end do

#ifdef WITH_PRECOND_S
         allocate(phiMat_fac(n_r_max,nLMBs2(1+rank)))
         bytes_allocated = bytes_allocated+n_r_max*nLMBs2(1+rank)*SIZEOF_DEF_REAL
#endif

#ifdef WITHOMP
         maxThreads=omp_get_max_threads()
#else
         maxThreads=1
#endif
         allocate( rhs1(n_r_max,2*lo_sub_map%sizeLMB2max,0:maxThreads-1) )
         bytes_allocated = bytes_allocated + n_r_max*lo_sub_map%sizeLMB2max*&
         &                 maxThreads*SIZEOF_DEF_COMPLEX

      allocate( lPhimat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

   end subroutine initialize_updatePhi
!------------------------------------------------------------------------------
   subroutine finalize_updatePhi
      !
      ! This subroutine deallocates the matrices involved in the time-advance of
      ! phi.
      !

      integer, pointer :: nLMBs2(:)
      integer :: ll

      deallocate( lPhimat )
         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

         do ll=1,nLMBs2(1+rank)
            call phiMat(ll)%finalize()
         end do

#ifdef WITH_PRECOND_S
         deallocate(phiMat_fac)
#endif
         deallocate( rhs1 )
     
   end subroutine finalize_updatePhi
!------------------------------------------------------------------------------
   subroutine updatePhi(phi, dphidt, tscheme)
      !
      !  Updates the phase field
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      complex(cp),       intent(inout) :: phi(llm:ulm,n_r_max) ! Chemical composition
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables:
      integer :: l1,m1          ! degree and order
      integer :: lm1,lm         ! position of (l,m) in array
      integer :: nLMB2,nLMB
      integer :: nR             ! counts radial grid points
      integer :: n_r_out        ! counts cheb modes

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: threadid,iChunk,nChunks,size_of_last_chunk,lmB0

      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      nLMB=1+rank

      !-- Now assemble the right hand side and store it in work_LMloc
      call tscheme%set_imex_rhs(work_LMloc, dphidt)

      !$omp parallel default(shared)

      !$omp single
      call solve_counter%start_count()
      !$omp end single
      ! one subblock is linked to one l value and needs therefore once the matrix
      !$omp single
      do nLMB2=1,nLMBs2(nLMB)
         ! this inner loop is in principle over the m values which belong to the
         ! l value
         !$omp task default(shared) &
         !$omp firstprivate(nLMB2) &
         !$omp private(lm,lm1,l1,m1,threadid) &
         !$omp private(nChunks,size_of_last_chunk,iChunk)
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         ! This task treats one l given by l1
         l1=lm22l(1,nLMB2,nLMB)

         if ( .not. lPhimat(l1) ) then
#ifdef WITH_PRECOND_S
            call get_phiMat(tscheme,l1,phiMat(nLMB2),phiMat_fac(:,nLMB2))
#else
            call get_phiMat(tscheme,l1,phiMat(nLMB2))
#endif
             lPhimat(l1)=.true.
         end if

         do iChunk=1,nChunks
            !$omp task default(shared) &
            !$omp firstprivate(iChunk) &
            !$omp private(lmB0,lm,lm1,m1,nR,n_r_out) &
            !$omp private(threadid)
#ifdef WITHOMP
            threadid = omp_get_thread_num()
#else
            threadid = 0
#endif
            lmB0=(iChunk-1)*chunksize

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)

               if ( l1 == 0 ) then
                  rhs1(1,2*lm-1,threadid)      =phi_top
                  rhs1(1,2*lm,threadid)        =0.0_cp
                  rhs1(n_r_max,2*lm-1,threadid)=phi_bot
                  rhs1(n_r_max,2*lm,threadid)  =0.0_cp
               else
                  rhs1(1,2*lm-1,threadid)      =0.0_cp
                  rhs1(1,2*lm,threadid)        =0.0_cp
                  rhs1(n_r_max,2*lm-1,threadid)=0.0_cp
                  rhs1(n_r_max,2*lm,threadid)  =0.0_cp
               end if
               do nR=2,n_r_max-1
                  rhs1(nR,2*lm-1,threadid)= real(work_LMloc(lm1,nR))
                  rhs1(nR,2*lm,threadid)  =aimag(work_LMloc(lm1,nR))
               end do

#ifdef WITH_PRECOND_S
               rhs1(:,2*lm-1,threadid)=phiMat_fac(:,nLMB2)*rhs1(:,2*lm-1,threadid)
               rhs1(:,2*lm,threadid)  =phiMat_fac(:,nLMB2)*rhs1(:,2*lm,threadid)
#endif

            end do

            call phiMat(nLMB2)%solve(rhs1(:,2*(lmB0+1)-1:2*(lm-1),threadid), &
                 &                   2*(lm-1-lmB0))

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( m1 > 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     phi(lm1,n_r_out)=cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                     &                      rhs1(n_r_out,2*lm,threadid),kind=cp)
                  end do
               else
                  do n_r_out=1,rscheme_oc%n_max
                     phi(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                     &                       0.0_cp,kind=cp)
                  end do
               end if
            end do
            !$omp end task
         end do
         !$omp taskwait
         !$omp end task
      end do     ! loop over lm blocks
      !$omp end single
      !$omp taskwait
      !$omp single
      call solve_counter%stop_count(l_increment=.false.)
      !$omp end single

      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      !$omp do private(n_r_out,lm1) collapse(2)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=llm,ulm
            phi(lm1,n_r_out)=zero
         end do
      end do
      !$omp end do

      !$omp end parallel

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dphidt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_phase_rhs_imp(phi, dphidt, 1, tscheme%l_imp_calc_rhs(1),  &
              &                l_in_cheb_space=.true.)
      else
         call get_phase_rhs_imp(phi, dphidt, tscheme%istage+1,            &
              &                tscheme%l_imp_calc_rhs(tscheme%istage+1),  &
              &                l_in_cheb_space=.true.)
      end if

   end subroutine updatePhi
!------------------------------------------------------------------------------
   subroutine get_phase_rhs_imp(phi, dphidt, istage, l_calc_lin, l_in_cheb_space)
      !
      ! This subroutine computes the linear terms which enter the r.h.s. of the
      ! equation for phase field. This is the LM-distributed version.
      !

      !-- Input variables
      integer,             intent(in) :: istage
      logical,             intent(in) :: l_calc_lin
      logical, optional,   intent(in) :: l_in_cheb_space

      !-- Output variable
      complex(cp),       intent(inout) :: phi(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables
      complex(cp) :: dphi(llm:ulm,n_r_max)
      logical :: l_in_cheb
      integer :: n_r, lm, start_lm, stop_lm, l
      real(cp) :: dL
      integer, pointer :: lm2l(:)

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      lm2l(1:lm_max) => lo_map%lm2l

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr(phi, dphi, work_LMloc, ulm-llm+1,start_lm-llm+1,  &
           &       stop_lm-llm+1,n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(phi,ulm-llm+1,start_lm-llm+1, &
                            &                 stop_lm-llm+1)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      if ( istage == 1 ) then
         !$omp do
         do n_r=1,n_r_max
            dphidt%old(:,n_r,istage)=5.0_cp/6.0_cp*stef*pr*phi(:,n_r)
         end do
         !$omp end do
      end if

      if ( l_calc_lin ) then

         !$omp do private(n_r,lm,l,dL)
         do n_r=1,n_r_max
            do lm=llm,ulm
               l = lm2l(lm)
               dL = real(l*(l+1),cp)
               dphidt%impl(lm,n_r,istage)=phaseDiffFac*( work_LMloc(lm,n_r) + &
               &                                two*or1(n_r) * dphi(lm,n_r) - &
               &                                 dL*or2(n_r) *  phi(lm,n_r) )
            end do
         end do
         !$omp end do

      end if

      !$omp end parallel

   end subroutine get_phase_rhs_imp
!------------------------------------------------------------------------------
   subroutine assemble_phase(phi, dphidt, tscheme)
      !
      ! This subroutine is used to assemble the phase field when an
      ! IMEX-RK with an assembly stage is employed.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),       intent(inout) :: phi(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dphidt

      !-- Local variables
      integer :: lm, l, m, n_r
      integer, pointer :: lm2l(:), lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      call tscheme%assemble_imex(work_LMloc, dphidt)

      !$omp parallel default(shared)
      !$omp do private(n_r,lm,m)
      do n_r=2,n_r_max
         do lm=llm,ulm
            m = lm2m(lm)
            if ( m == 0 ) then
               phi(lm,n_r)=cmplx(real(work_LMloc(lm,n_r)),0.0_cp,cp) * &
               &           6.0_cp/5.0_cp/stef/pr
            else
               phi(lm,n_r)=work_LMloc(lm,n_r)*6.0_cp/5.0_cp/stef/pr
            end if
         end do
      end do
      !$omp end do

      !-- Boundary conditions
      if ( l_full_sphere) then
         if ( ktopphi == 1 ) then ! Dirichlet
            !$omp do private(lm,l)
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l == 0 ) then
                  call rscheme_oc%robin_bc(0.0_cp, one, cmplx(phi_top,0.0_cp,cp), &
                       &                   one, 0.0_cp, cmplx(phi_bot,0.0_cp,cp), &
                       &                   phi(lm,:))
               else
                  call rscheme_oc%robin_bc(0.0_cp, one, zero, 0.0_cp, one, &
                       &                   zero, phi(lm,:))
               end if
            end do
            !$omp end do
         else ! Neummann
            !$omp do private(lm,l)
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l == 0 ) then
                  call rscheme_oc%robin_bc(one, 0.0_cp, zero, one, 0.0_cp, &
                       &                   zero, phi(lm,:))
               else
                  call rscheme_oc%robin_bc(one, 0.0_cp, zero, 0.0_cp, one, &
                       &                   zero, phi(lm,:))
               end if
            end do
            !$omp end do
         end if

      else ! Spherical shell

         if ( ktopphi==1 .and. kbotphi==1 ) then
            !-- Boundary conditions: Dirichlet on both sides
            !$omp do private(lm,l)
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l == 0 ) then
                  call rscheme_oc%robin_bc(0.0_cp, one, cmplx(phi_top,0.0_cp,cp), &
                       &                   0.0_cp, one, cmplx(phi_bot,0.0_cp,cp), &
                       &                   phi(lm,:))
               else
                  call rscheme_oc%robin_bc(0.0_cp, one, zero, 0.0_cp, one, &
                       &                   zero, phi(lm,:))
               end if
            end do
            !$omp end do
         else if ( ktopphi==1 .and. kbotphi /= 1 ) then
            !$omp do private(lm,l)
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l == 0 ) then
                  call rscheme_oc%robin_bc(0.0_cp, one, cmplx(phi_top,0.0_cp,cp), &
                       &                   one, 0.0_cp, zero, phi(lm,:))
               else
                  call rscheme_oc%robin_bc(0.0_cp, one, zero, one, 0.0_cp, &
                       &                   zero, phi(lm,:))
               end if
            end do
            !$omp end do
         else if ( ktopphi/=1 .and. kbotphi == 1 ) then
            !$omp do private(lm,l)
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l == 0 ) then
                  call rscheme_oc%robin_bc(one, 0.0_cp, zero, 0.0_cp, one, &
                       &                   cmplx(phi_bot,0.0_cp,cp), phi(lm,:))
               else
                  call rscheme_oc%robin_bc(one, 0.0_cp, zero, 0.0_cp, one, &
                       &                   zero, phi(lm,:))
               end if
            end do
            !$omp end do
         else if ( ktopphi/=1 .and. kbotphi /= 1 ) then
            !-- Boundary conditions: Neuman on both sides
            !$omp do private(lm)
            do lm=llm,ulm
               call rscheme_oc%robin_bc(one, 0.0_cp, zero, one, 0.0_cp, zero, phi(lm,:))
            end do
            !$omp end do
         end if

      end if
      !$omp end parallel

      call get_phase_rhs_imp(phi, dphidt, 1, tscheme%l_imp_calc_rhs(1), .false.)

   end subroutine assemble_phase
!------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_phiMat(tscheme,l,phiMat,phiMat_fac)
#else
   subroutine get_phiMat(tscheme,l,phiMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrices
      !  phiMat(i,j) for the equation for phase field.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step
      integer,             intent(in) :: l

      !-- Output variables
      class(type_realmat), intent(inout) :: phiMat
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: phiMat_fac(n_r_max)
#endif

      !-- Local variables:
      integer :: info, nR_out, nR
      real(cp) :: dLh
      real(cp) :: dat(n_r_max,n_r_max)

      dLh=real(l*(l+1),kind=cp)

      !----- Boundary conditions:
      if ( ktopphi == 1 ) then ! Dirichlet
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      else ! Neumann
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         !dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         if ( l == 0 ) then
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         else
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         end if
      else
         if ( kbotphi == 1 ) then ! Dirichlet
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         else
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         end if
      end if

      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)      =0.0_cp
            dat(n_r_max,nR_out)=0.0_cp
         end do
      end if

      !----- Bulk points
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            dat(nR,nR_out)= rscheme_oc%rnorm * (                                &
            &               5.0_cp/6.0_cp*stef*pr* rscheme_oc%rMat(nR,nR_out) - &
            &  tscheme%wimp_lin(1)*phaseDiffFac*(rscheme_oc%d2rMat(nR,nR_out) + &
            &                     two*or1(nR)*    rscheme_oc%drMat(nR,nR_out) - &
            &                     dLh*or2(nR)*     rscheme_oc%rMat(nR,nR_out) ) )
         end do
      end do

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_S
      ! compute the linesum of each line
      do nR=1,n_r_max
         phiMat_fac(nR)=one/maxval(abs(dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*phiMat_fac(nR)
      end do
#endif

      !-- Array copy
      call phiMat%set_data(dat)

      !----- LU decomposition:
      call phiMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix phiMat!')

   end subroutine get_phiMat
 end module updatePhi_mod
