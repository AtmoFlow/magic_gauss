module updateV_mod
   !
   ! This module handles the time advance of the chemical composition v.
   ! It contains the computation of the implicit terms and the linear
   ! solves.
   !

   use omp_lib
   use precision_mod
   use truncation, only: n_r_max, lm_max, l_max
   use radial_data, only: n_r_icb, n_r_cmb, nRstart, nRstop
   use radial_functions, only: orho1, or1, or2, beta, rscheme_oc, r, l_R
   use num_param, only: dct_counter, solve_counter
   use init_fields, only: tope, bote
   use blocking, only: lo_map, lo_sub_map, llm, ulm, st_map
   use horizontal_data, only: hdif_v
   use parallel_mod, only: rank, chunksize, n_procs, get_openmp_blocks
   use radial_der, only: get_ddr, get_dr, get_dr_Rloc
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
   class(type_realmat), pointer :: vMat(:)
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: vMat_fac(:,:)
#endif
   logical, public, allocatable :: lvmat(:)

   public :: initialize_updatev, finalize_updatev, updatev, assemble_efield,  &
   &         finish_exp_efield, get_efield_rhs_imp

contains

   subroutine initialize_updatev

      integer :: ll, n_bands
      integer, pointer :: nLMBs2(:)

         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

            allocate( type_densemat :: vMat(nLMBs2(1+rank)) )

            do ll=1,nLMBs2(1+rank)
               call vMat(ll)%initialize(n_r_max,n_r_max,l_pivot=.true.)
            end do

#ifdef WITH_PRECOND_S
         allocate(vMat_fac(n_r_max,nLMBs2(1+rank)))
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
      

      allocate( lvmat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

   end subroutine initialize_updateV
!------------------------------------------------------------------------------
   subroutine finalize_updateV
      !
      ! This subroutine deallocates the matrices involved in the time-advance of
      ! v.
      !

      integer, pointer :: nLMBs2(:)
      integer :: ll

      deallocate( lvmat )
         nLMBs2(1:n_procs) => lo_sub_map%nLMBs2

         do ll=1,nLMBs2(1+rank)
            call vMat(ll)%finalize()
         end do

#ifdef WITH_PRECOND_S
         deallocate(vMat_fac)
#endif
         deallocate( rhs1 )

   end subroutine finalize_updateV
!------------------------------------------------------------------------------
   subroutine updateV(v, dv, Et_LMloc, tscheme)
      !
      !  Updates the chemical composition field s and its radial derivative.
      !

      !-- Input of variables:
      class(type_tscheme), intent(in) :: tscheme

      !-- Input/output of scalar fields:
      complex(cp),       intent(inout) :: v(llm:ulm,n_r_max) ! electric potential
!      type(type_tarray), intent(inout) :: dvdt
      !-- Output: dv
      complex(cp),       intent(out) :: dv(llm:ulm,n_r_max) ! Radial derivative of v
      complex(cp),       intent(in) :: Et_LMloc(llm:ulm,n_r_max) ! Radial derivative of v


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
!      call tscheme%set_imex_rhs(work_LMloc, dvdt)

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

         if ( .not. lvmat(l1) ) then
#ifdef WITH_PRECOND_S
            call get_vMat(tscheme,l1,hdif_v(l1),vMat(nLMB2),vMat_fac(:,nLMB2))
#else
            call get_vMat(tscheme,l1,hdif_v(l1),vMat(nLMB2))
#endif
             lvmat(l1)=.true.
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
               
               rhs1(1,2*lm-1,threadid)      = real(tope(l1,m1))
               rhs1(1,2*lm,threadid)        =aimag(tope(l1,m1))
               rhs1(n_r_max,2*lm-1,threadid)= real(bote(l1,m1))
               rhs1(n_r_max,2*lm,threadid)  =aimag(bote(l1,m1))
               do nR=2,n_r_max-1
                  rhs1(nR,2*lm-1,threadid)= real(Et_LMloc(lm1,nR))
                  rhs1(nR,2*lm,threadid)  =aimag(Et_LMloc(lm1,nR))
               end do

#ifdef WITH_PRECOND_S
               rhs1(:,2*lm-1,threadid)=vMat_fac(:,nLMB2)*rhs1(:,2*lm-1,threadid)
               rhs1(:,2*lm,threadid)  =vMat_fac(:,nLMB2)*rhs1(:,2*lm,threadid)
#endif

            end do

            call vMat(nLMB2)%solve(rhs1(:,2*(lmB0+1)-1:2*(lm-1),threadid), &
                 &                  2*(lm-1-lmB0))

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               lm1=lm22lm(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( m1 > 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     v(lm1,n_r_out)=cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                     &                     rhs1(n_r_out,2*lm,threadid),kind=cp)
                  end do
               else
                  do n_r_out=1,rscheme_oc%n_max
                     v(lm1,n_r_out)= cmplx(rhs1(n_r_out,2*lm-1,threadid), &
                     &                      0.0_cp,kind=cp)
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
            v(lm1,n_r_out)=zero
         end do
      end do
      !$omp end do

      !$omp end parallel

      !-- Roll the arrays before filling again the first block
!      call tscheme%rotate_imex(dvdt)

!      -- Calculation of the implicit part
      if ( tscheme%istage == tscheme%nstages ) then
         call get_efield_rhs_imp(v, dv, 1, tscheme%l_imp_calc_rhs(1),  &
              &                l_in_cheb_space=.true.)
      else
         call get_efield_rhs_imp(v, dv, tscheme%istage+1,          &
              &                tscheme%l_imp_calc_rhs(tscheme%istage+1),  &
              &                l_in_cheb_space=.true.)
      end if

   end subroutine updateV
!------------------------------------------------------------------------------
   subroutine finish_exp_efield(w, dVvrLM, dv_exp_last)
      !
      ! This subroutine completes the computation of the advection term which
      ! enters the composition equation by taking the radial derivative. This is
      ! the LM-distributed version.
      !

      !-- Input variables
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVvrLM(llm:ulm,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: dv_exp_last(llm:ulm,n_r_max)

      !-- Local variables
      real(cp) :: dLh
      integer :: n_r, start_lm, stop_lm, l, lm

      !$omp parallel default(shared) private(start_lm, stop_lm)
      start_lm=llm; stop_lm=ulm
      call get_openmp_blocks(start_lm,stop_lm)
      call get_dr( dVvrLM, work_LMloc, ulm-llm+1, start_lm-llm+1,  &
           &       stop_lm-llm+1, n_r_max, rscheme_oc, nocopy=.true. )
      !$omp barrier

      !$omp do private(lm,l)
      do n_r=1,n_r_max
         do lm=llm,ulm
            l = lo_map%lm2l(lm)
            if ( l > l_R(n_r) ) cycle
            dv_exp_last(lm,n_r)=orho1(n_r)*( dv_exp_last(lm,n_r)-   &
            &                         or2(n_r)*work_LMloc(lm,n_r) )
         end do
      end do
      !$omp end do

      !$omp end parallel

   end subroutine finish_exp_efield
!------------------------------------------------------------------------------
   subroutine get_efield_rhs_imp(v, dv, istage, l_calc_lin, l_in_cheb_space)
      !
      ! This subroutine computes the linear terms which enter the r.h.s. of the
      ! equation for composition. This is the LM-distributed version.
      !

      !-- Input variables
      integer,             intent(in) :: istage
      logical,             intent(in) :: l_calc_lin
      logical, optional,   intent(in) :: l_in_cheb_space

      !-- Output variable
      complex(cp),       intent(inout) :: v(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dv(llm:ulm,n_r_max)
!      type(type_tarray), intent(inout) :: dvdt

      !-- Local variables
      logical :: l_in_cheb
      integer :: n_r, lm, start_lm, stop_lm, l1
      real(cp) :: dL
      integer, pointer :: lm2l(:),lm2m(:)

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
      call get_ddr(v, dv, work_LMloc, ulm-llm+1,start_lm-llm+1,  &
           &       stop_lm-llm+1,n_r_max, rscheme_oc, l_dct_in=.not. l_in_cheb)
      if ( l_in_cheb ) call rscheme_oc%costf1(v,ulm-llm+1,start_lm-llm+1, &
                            &                 stop_lm-llm+1)
      !$omp barrier
      !$omp single
      call dct_counter%stop_count(l_increment=.false.)
      !$omp end single


      !$omp end parallel

   end subroutine get_efield_rhs_imp
!------------------------------------------------------------------------------
   subroutine assemble_efield(v, dv, dvdt, tscheme)
      !
      ! This subroutine is used to assemble the chemical composition when an
      ! IMEX-RK with an assembly stage is employed. Non-Dirichlet boundary
      ! conditions are handled using Canuto (1986) approach. This is the LM
      ! distributed version.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),       intent(inout) :: v(llm:ulm,n_r_max)
      complex(cp),       intent(out) :: dv(llm:ulm,n_r_max)
      type(type_tarray), intent(inout) :: dvdt

      !-- Local variables
      integer :: lm, l1, m1, n_r
      integer, pointer :: lm2l(:), lm2m(:)

      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      call tscheme%assemble_imex(work_LMloc, dvdt)

      !$omp parallel default(shared)
      !$omp do private(n_r,lm,m1)
      do n_r=2,n_r_max
         do lm=llm,ulm
            m1 = lm2m(lm)
            if ( m1 == 0 ) then
               v(lm,n_r)=cmplx(real(work_LMloc(lm,n_r)),0.0_cp,cp)
            else
               v(lm,n_r)=work_LMloc(lm,n_r)
            end if
         end do
      end do
      !$omp end do

         !-- Boundary conditions
            !$omp do private(lm,l1,m1)
            do lm=llm,ulm
               l1 = lm2l(lm)
               m1 = lm2m(lm)
!               call rscheme_oc%robin_bc(0.0_cp, one, topv(l1,m1), 0.0_cp, one, &
!                    &                   botv(l1,m1), v(lm,:))
            end do
            !$omp end do

      !$omp end parallel

      call get_efield_rhs_imp(v, dv, 1, tscheme%l_imp_calc_rhs(1), .false.)

   end subroutine assemble_efield
!------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_vMat(tscheme,l,hdif,vMat,vMat_fac)
#else
   subroutine get_vMat(tscheme,l,hdif,vMat)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrices
      !  vMat(i,j) for the equation for the chemical composition.
      !

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme        ! time step
      real(cp),            intent(in) :: hdif
      integer,             intent(in) :: l

      !-- Output variables
      class(type_realmat), intent(inout) :: vMat
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: vMat_fac(n_r_max)
#endif

      !-- Local variables:
      integer :: info, nR_out, nR
      real(cp) :: dLh
      real(cp) :: dat(n_r_max,n_r_max)

      dLh=real(l*(l+1),kind=cp)

      !----- Boundary conditions:
      dat(1,:)=rscheme_oc%rnorm*rscheme_oc%rMat(1,:)
      dat(n_r_max,:)=rscheme_oc%rnorm*rscheme_oc%rMat(n_r_max,:)


      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            dat(1,nR_out)      =0.0_cp
            dat(n_r_max,nR_out)=0.0_cp
         end do
      end if

      !----- Bulk points
      do nR_out=1,n_r_max
         do nR=2,n_r_max-1
            dat(nR,nR_out)= rscheme_oc%rnorm * (                           &
                    & rscheme_oc%d2rMat(nR,nR_out) + &
                            & two*or1(nR)*rscheme_oc%drMat(nR,nR_out) - &
                            & dLh*or2(nR)*rscheme_oc%rMat(nR,nR_out) )
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
         vMat_fac(nR)=one/maxval(abs(dat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         dat(nR,:) = dat(nR,:)*vMat_fac(nR)
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
      call vMat%set_data(dat)

      !----- LU decomposition:
      call vMat%prepare(info)
      if ( info /= 0 ) call abortRun('Singular matrix vMat!')

   end subroutine get_vMat
!-----------------------------------------------------------------------------
end module updateV_mod
