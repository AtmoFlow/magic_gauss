module updateS_mod
   !
   ! This module handles the time advance of the entropy s.
   ! It contains the computation of the implicit terms and the linear
   ! solves.
   !

   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, lm_max, l_max
   use radial_data, only: n_r_cmb, n_r_icb, nRstart, nRstop
   use radial_functions, only: orho1, or1, or2, beta, dentropy0, rscheme_oc,  &
       &                       kappa, dLkappa, dLtemp0, temp0, r, l_R
   use physical_parameters, only: opr, kbots, ktops, stef
   use num_param, only: dct_counter, solve_counter
   use init_fields, only: tops, bots
   use blocking, only: lo_map, lo_sub_map, llm, ulm, st_map
   use horizontal_data, only: hdif_S
   use logic, only: l_anelastic_liquid, l_phase_field, &
       &            l_full_sphere
   use parallel_mod
   use radial_der, only: get_ddr, get_dr, get_dr_Rloc, get_ddr_ghost, &
       &                 exch_ghosts, bulk_to_ghost
   use fields, only:  work_LMloc
   use constants, only: zero, one, two
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
   complex(cp), allocatable, public :: s_ghost(:,:)
   class(type_realmat), pointer :: sMat(:)
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: sMat_fac(:,:)
#endif
   logical, public, allocatable :: lSmat(:)

   type(type_tri_par), public :: sMat_FD
   real(cp), allocatable :: fd_fac_top(:), fd_fac_bot(:)

   integer :: maxThreads

   public :: initialize_updateS, updateS, finalize_updateS, assemble_entropy,  &
   &         finish_exp_entropy, get_entropy_rhs_imp

contains

   subroutine initialize_updateS
      !
      ! This subroutine allocates the arrays involved in the time-advance of the
      ! entropy/temperature equation.
      !

      integer, pointer :: nLMBs2(:)
      integer :: ll,n_bands


      nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

            allocate( type_densemat :: sMat(nLMBs2(1+rank)) )

            do ll=1,nLMBs2(1+rank)
               call sMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
            end do

#ifdef WITH_PRECOND_S
         allocate(sMat_fac(n_r_max,nLMBs2(1+rank)))
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
    
      allocate( lSmat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

   end subroutine initialize_updateS
!------------------------------------------------------------------------------
   subroutine finalize_updateS
      !
      ! Memory deallocation of updateS module
      !

      integer, pointer :: nLMBs2(:)
      integer :: ll

      deallocate( lSmat)
         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

         do ll=1,nLMBs2(1+rank)
            call sMat(ll)%finalize()
         end do
         deallocate(rhs1)

#ifdef WITH_PRECOND_S
         deallocate( sMat_fac )
#endif

   end subroutine finalize_updateS
!------------------------------------------------------------------------------
   subroutine updateS(s, ds, dsdt, phi, tscheme)
      !
      !  Updates the entropy field s and its radial derivative.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(in) :: phi(llm:ulm,n_r_max) ! Phase field

      !-- Input/output of scalar fields:
      complex(cp),       intent(inout) :: s(llm:ulm,n_r_max) ! Entropy
      type(type_tarray), intent(inout) :: dsdt
      !-- Output: ds
      complex(cp),       intent(out) :: ds(llm:ulm,n_r_max) ! Radial derivative of entropy

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
      call tscheme%set_imex_rhs(work_LMloc, dsdt)

      !$omp parallel default(shared)

      if ( l_phase_field ) then
         !-- Add the last remaining term to assemble St*\partial \phi/\partial t
         !$omp do private(nR,lm)
         do nR=1,n_r_max
            do lm=llm,ulm
               work_LMloc(lm,nR)=work_LMloc(lm,nR)+stef*phi(lm,nR)
            end do
         end do
         !$omp end do
      end if

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

         if ( .not. lSmat(l1) ) then
#ifdef WITH_PRECOND_S
            call get_sMat(tscheme,l1,hdif_S(l1),sMat(nLMB2),sMat_fac(:,nLMB2))
#else
            call get_sMat(tscheme,l1,hdif_S(l1),sMat(nLMB2))
#endif
            lSmat(l1)=.true.
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

               rhs1(1,2*lm-1,threadid)      = real(tops(l1,m1))
               rhs1(1,2*lm,threadid)        =aimag(tops(l1,m1))
               rhs1(n_r_max,2*lm-1,threadid)= real(bots(l1,m1))
               rhs1(n_r_max,2*lm,threadid)  =aimag(bots(l1,m1))

               do nR=2,n_r_max-1
                  rhs1(nR,2*lm-1,threadid)= real(work_LMloc(lm1,nR))
                  rhs1(nR,2*lm,threadid)  =aimag(work_LMloc(lm1,nR))
               end do

#ifdef WITH_PRECOND_S
               rhs1(:,2*lm-1,threadid)=sMat_fac(:,nLMB2)*rhs1(:,2*lm-1,threadid)
               rhs1(:,2*lm,threadid)  =sMat_fac(:,nLMB2)*rhs1(:,2*lm,threadid)
#endif
            end do

            call sMat(nLMB2)%solve(rhs1(:,2*(lmB0+1)-1:2*(lm-1),threadid), &
                 &                 2*(lm-1-lmB0))

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( m1 > 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     s(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                     &                     rhs1(n_r_out,2*lm,threadid),kind=cp)
                  end do
               else
                  do n_r_out=1,rscheme_oc%n_max
                     s(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                     &                     0.0_cp,kind=cp)
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
            s(lm1,n_r_out)=zero
         end do
      end do
      !$omp end do
      !$omp end parallel

      !-- Roll the arrays before filling again the first block
      call tscheme%rotate_imex(dsdt)

      !-- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_entropy_rhs_imp(s, ds, dsdt, phi, 1, tscheme%l_imp_calc_rhs(1), &
              &                   l_in_cheb_space=.true.)
      else
         call get_entropy_rhs_imp(s, ds, dsdt, phi, tscheme%istage+1,       &
              &                   tscheme%l_imp_calc_rhs(tscheme%istage+1), &
              &                   l_in_cheb_space=.true.)
      end if

   end subroutine updateS
!------------------------------------------------------------------------------
   subroutine finish_exp_entropy(w, dVSrLM, ds_exp_last)
      !
      ! This subroutine completes the computation of the advection term by
      ! computing the radial derivative (LM-distributed variant).
      !

      !-- Input variables
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVSrLM(llm:ulm,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: ds_exp_last(llm:ulm,n_r_max)

      !-- Local variables
      real(cp) :: dL
      integer :: n_r, lm, start_lm, stop_lm, l
      integer, pointer :: lm2l(:)

      lm2l(1:lm_max) => lo_map%lm2l

      !$omp parallel default(shared) private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)
      call get_dr( dVSrLM, work_LMloc, ulm-llm+1, start_lm-llm+1,  &
           &       stop_lm-llm+1, n_r_max, rscheme_oc, nocopy=.true. )
      !$omp barrier

      if ( l_anelastic_liquid ) then
         !$omp do private(n_r,l,lm,dL)
         do n_r=1,n_r_max
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l > l_R(n_r) ) cycle
               dL = real(l*(l+1),cp)
               ds_exp_last(lm,n_r)=orho1(n_r)*     ds_exp_last(lm,n_r) - &
               &        or2(n_r)*orho1(n_r)*        work_LMloc(lm,n_r) + &
               &       or2(n_r)*orho1(n_r)*dLtemp0(n_r)*dVSrLM(lm,n_r) - &
               &        dL*or2(n_r)*orho1(n_r)*temp0(n_r)*dentropy0(n_r)*&
               &                                             w(lm,n_r)
            end do
         end do
         !$omp end do
      else
         !$omp do private(n_r,l,dL,lm)
         do n_r=1,n_r_max
            do lm=llm,ulm
               l = lm2l(lm)
               if ( l > l_R(n_r) ) cycle
               dL = real(l*(l+1),cp)
               ds_exp_last(lm,n_r)=orho1(n_r)*(      ds_exp_last(lm,n_r)- &
               &                             or2(n_r)*work_LMloc(lm,n_r)- &
               &                    dL*or2(n_r)*dentropy0(n_r)*w(lm,n_r))
            end do
         end do
         !$omp end do
      end if
      !$omp end parallel

   end subroutine finish_exp_entropy
!-----------------------------------------------------------------------------
   subroutine get_entropy_rhs_imp(s, ds, dsdt, phi, istage, l_calc_lin, l_in_cheb_space)
      !
      ! This subroutine computes the linear terms that enters the r.h.s.. This is
      ! used with LM-distributed
      !

      !-- Input variables
      integer,           intent(in) :: istage
      logical,           intent(in) :: l_calc_lin
      logical, optional, intent(in) :: l_in_cheb_space
      complex(cp),       intent(in) :: phi(llm:ulm,n_r_max)

      !-- Output variable
      complex(cp),       intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: ds(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dsdt

      !-- Local variables
      integer :: n_r, lm, start_lm, stop_lm, l1
      logical :: l_in_cheb
      integer, pointer :: lm2l(:),lm2m(:)
      real(cp) :: dL

      if ( present(l_in_cheb_space) ) then
         l_in_cheb = l_in_cheb_space
      else
         l_in_cheb = .false.
      end if

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      !$omp parallel default(shared)  private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)

      !$omp single
      call dct_counter%start_count()
      !$omp end single
      call get_ddr(s, ds, work_LMloc, ulm-llm+1,start_lm-llm+1,  &
           &       stop_lm-llm+1,n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(s,ulm-llm+1,start_lm-llm+1, &
                            &                 stop_lm-llm+1)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single

      if ( istage == 1 ) then
         !$omp do private(n_r)
         do n_r=1,n_r_max
            dsdt%old(:,n_r,istage)=s(:,n_r)
         end do
         !$omp end do
         if ( l_phase_field ) then
            !$omp do private(n_r)
            do n_r=1,n_r_max
               dsdt%old(:,n_r,istage)=dsdt%old(:,n_r,istage)-stef*phi(:,n_r)
            end do
            !$omp end do
         end if
      end if

      if ( l_calc_lin ) then

         !-- Calculate explicit time step part:
         if ( l_anelastic_liquid ) then
            !$omp do private(n_r,lm,l1,dL)
            do n_r=1,n_r_max
               do lm=llm,ulm
                  l1 = lm2l(lm)
                  dL = real(l1*(l1+1),cp)
                  dsdt%impl(lm,n_r,istage)=  opr*hdif_S(l1)* kappa(n_r) *  (    &
                  &                                          work_LMloc(lm,n_r) &
                  &     + ( beta(n_r)+two*or1(n_r)+dLkappa(n_r) ) *  ds(lm,n_r) &
                  &                                     - dL*or2(n_r)*s(lm,n_r) )
               end do
            end do
            !$omp end do
         else
            !$omp do private(n_r,lm,l1,dL)
            do n_r=1,n_r_max
               do lm=llm,ulm
                  l1 = lm2l(lm)
                  dL = real(l1*(l1+1),cp)
                  dsdt%impl(lm,n_r,istage)=  opr*hdif_S(l1)*kappa(n_r) *   (       &
                  &                                        work_LMloc(lm,n_r)      &
                  &        + ( beta(n_r)+dLtemp0(n_r)+two*or1(n_r)+dLkappa(n_r) )  &
                  &                                              * ds(lm,n_r)      &
                  &        - dL*or2(n_r)                         *  s(lm,n_r) )
               end do
            end do
            !$omp end do
         end if

      end if
      !$omp end parallel

   end subroutine get_entropy_rhs_imp
!-----------------------------------------------------------------------------
   subroutine assemble_entropy(s, ds, dsdt, phi, tscheme)
      !
      ! This subroutine is used to assemble the entropy/temperature at assembly
      ! stages of IMEX-RK time schemes. This is used when LM is distributed.
      !

      !-- Input variable
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(in) :: phi(llm:ulm,n_r_max)

      !-- Output variables
      complex(cp),       intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: ds(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dsdt

      !-- Local variables
      integer :: lm, l1, m1, n_r
      integer, pointer :: lm2l(:), lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      call tscheme%assemble_imex(work_LMloc, dsdt)

      !$omp parallel default(shared)

      !-- In case phase field is used it needs to be substracted from work_LMloc
      !-- since time advance handles \partial/\partial t (T-St*Phi)
      if ( l_phase_field ) then
         !$omp do private(n_r,lm)
         do n_r=1,n_r_max
            do lm=llm,ulm
               work_LMloc(lm,n_r)=work_LMloc(lm,n_r)+stef*phi(lm,n_r)
            end do
         end do
         !$omp end do
      end if

      !$omp do private(n_r,lm,m1)
      do n_r=2,n_r_max
         do lm=llm,ulm
            m1 = lm2m(lm)
            if ( m1 == 0 ) then
               s(lm,n_r)=cmplx(real(work_LMloc(lm,n_r)),0.0_cp,cp)
            else
               s(lm,n_r)=work_LMloc(lm,n_r)
            end if
         end do
      end do
      !$omp end do

      !-- Get the boundary points using Canuto (1986) approach
      if ( l_full_sphere) then
         if ( ktops == 1 ) then ! Fixed entropy at the outer boundary
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               if ( l1 == 0 ) then
                  call rscheme_oc%robin_bc(0.0_cp, one, tops(l1,m1), one, 0.0_cp, &
                       &                   bots(l1,m1), s(lm,:))
               else
                  call rscheme_oc%robin_bc(0.0_cp, one, tops(l1,m1), 0.0_cp, one, &
                       &                   bots(l1,m1), s(lm,:))
               end if
            end do
            !$omp end do
         else ! Fixed flux at the outer boundary
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               if ( l1 == 0 ) then
                  call rscheme_oc%robin_bc(one, 0.0_cp, tops(l1,m1), one, 0.0_cp, &
                       &                   bots(l1,m1), s(lm,:))
               else
                  call rscheme_oc%robin_bc(one, 0.0_cp, tops(l1,m1), 0.0_cp, one, &
                       &                   bots(l1,m1), s(lm,:))
               end if
            end do
            !$omp end do
         end if
      else ! Spherical shell
         !-- Boundary conditions
         if ( ktops==1 .and. kbots==1 ) then ! Dirichlet on both sides
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               call rscheme_oc%robin_bc(0.0_cp, one, tops(l1,m1), 0.0_cp, one, &
                    &                   bots(l1,m1), s(lm,:))
            end do
            !$omp end do
         else if ( ktops==1 .and. kbots /= 1 ) then ! Dirichlet: top and Neumann: bot
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               call rscheme_oc%robin_bc(0.0_cp, one, tops(l1,m1), one, 0.0_cp, &
                    &                   bots(l1,m1), s(lm,:))
            end do
            !$omp end do
         else if ( kbots==1 .and. ktops /= 1 ) then ! Dirichlet: bot and Neumann: top
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               call rscheme_oc%robin_bc(one, 0.0_cp, tops(l1,m1), 0.0_cp, one, &
                    &                   bots(l1,m1), s(lm,:))
            end do
            !$omp end do
         else if ( kbots /=1 .and. kbots /= 1 ) then ! Neumann on both sides
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
               call rscheme_oc%robin_bc(one, 0.0_cp, tops(l1,m1), one, 0.0_cp, &
                    &                   bots(l1,m1), s(lm,:))
            end do
            !$omp end do
         end if
      end if
      !$omp end parallel

      !-- Finally call the construction of the implicit terms for the first stage
      !-- of next iteration
      call get_entropy_rhs_imp(s, ds, dsdt, phi, 1, tscheme%l_imp_calc_rhs(1), .false.)

   end subroutine assemble_entropy
!-------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_sMat(tscheme,l,hdif,sMat,sMat_fac)
#else
   subroutine get_sMat(tscheme,l,hdif,sMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrices
      !  sMat(i,j) and s0mat for the entropy equation.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step
      real(cp),            intent(in) :: hdif
      integer,             intent(in) :: l

      !-- Output variables
      class(type_realmat), intent(inout) :: sMat
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: sMat_fac(n_r_max)
#endif

      !-- Local variables:
      integer :: info,nR_out,nR
      real(cp) :: dLh
      real(cp) :: dat(n_r_max,n_r_max)

      dLh=real(l*(l+1),kind=cp)

      !----- Boundary conditions:
      if ( ktops == 1 ) then
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      else
         dat(1,:)=rscheme_oc%rnorm*rscheme_oc%drMat(1,:)
      end if

      if ( l_full_sphere ) then
         if ( l == 0 ) then
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%drMat(n_r_max,:)
         else
            dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)
         end if
      else
         if ( kbots == 1 ) then
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

      !----- Bulk points:
      if ( l_anelastic_liquid ) then
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               dat(nR,nR_out)= rscheme_oc%rnorm * (                         &
               &   rscheme_oc%rMat(nR,nR_out)-tscheme%wimp_lin(1)*opr*hdif* &
               &                kappa(nR)*(  rscheme_oc%d2rMat(nR,nR_out) + &
               &( beta(nR)+two*or1(nR)+dLkappa(nR) )*                       &
               &                              rscheme_oc%drMat(nR,nR_out) - &
               &      dLh*or2(nR)*             rscheme_oc%rMat(nR,nR_out) ) )
            end do
         end do
      else
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               dat(nR,nR_out)= rscheme_oc%rnorm * (                         &
               &   rscheme_oc%rMat(nR,nR_out)-tscheme%wimp_lin(1)*opr*hdif* &
               &                kappa(nR)*(  rscheme_oc%d2rMat(nR,nR_out) + &
               & ( beta(nR)+dLtemp0(nR)+                                    &
               &   two*or1(nR)+dLkappa(nR) )* rscheme_oc%drMat(nR,nR_out) - &
               &      dLh*or2(nR)*             rscheme_oc%rMat(nR,nR_out) ) )
            end do
         end do
      end if

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         dat(nR,1)      =rscheme_oc%boundary_fac*dat(nR,1)
         dat(nR,n_r_max)=rscheme_oc%boundary_fac*dat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_S
      ! compute the linesum of each line
      do nR=1,n_r_max
         sMat_fac(nR)=one/maxval(abs(dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*sMat_fac(nR)
      end do
#endif

#ifdef MATRIX_CHECK
      block

      integer :: i,j
      real(cp) :: rcond
      integer ::ipiv(n_r_max),iwork(n_r_max)
      real(cp) :: work(4*n_r_max),anorm,linesum
      real(cp) :: temp_Mat(n_r_max,n_r_max)
      integer,save :: counter=0
      integer :: filehandle
      character(len=100) :: filename

      ! copy the sMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "sMat_",l,"_",counter,".dat"
      open(newunit=filehandle,file=trim(filename))
      counter= counter+1

      do i=1,n_r_max
         do j=1,n_r_max
            write(filehandle,"(2ES20.12,1X)",advance="no") dat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_Mat=dat
      anorm = 0.0_cp
      do i=1,n_r_max
         linesum = 0.0_cp
         do j=1,n_r_max
            linesum = linesum + abs(temp_Mat(i,j))
         end do
         if (linesum  >  anorm) anorm=linesum
      end do
      !write(*,"(A,ES20.12)") "anorm = ",anorm
      ! LU factorization
      call dgetrf(n_r_max,n_r_max,temp_Mat,n_r_max,ipiv,info)
      ! estimate the condition number
      call dgecon('I',n_r_max,temp_Mat,n_r_max,anorm,rcond,work,iwork,info)
      write(*,"(A,I3,A,ES11.3)") "inverse condition number of sMat for l=",l," is ",rcond

      end block
#endif

      !-- Array copy
      call sMat%set_data(dat)

      !-- LU decomposition:
      call sMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix sMat!')

   end subroutine get_sMat
end module updateS_mod
