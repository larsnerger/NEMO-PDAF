!> Ensemble Initialisation
!!
!! This routine calls the routines for initialising the ensemble.
!! 
!! Separate calls are made for the 2D and 3D state variables to
!! allow for differences in how these variables are initialised.
!! 
!! The routine is called when the filter is initialized in
!! `PDAF_filter_init`.
!! 
!! The routine is called by all filter processes and
!! initializes the ensemble for the *PE-local domain*.
!! 
!!  **Calling Sequence**
!! 
!!  - Called from: `init_pdaf/PDAF_init` (PDAF module)
!!
subroutine init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

  use mod_kind_pdaf
  use mod_parallel_pdaf, &
       only: mype_filter
  use mod_assimilation_pdaf, &
       only: dim_state, type_ens_init, type_central_state, ensscale, &
       coupling_nemo, screen
  use mod_io_pdaf, &
       only: path_inistate, path_ens, file_ens, file_covar, &
             read_state_mv, &
             read_ens_mv_loop, read_ens, gen_ens_mv

  implicit none

! *** Arguments ***
  integer, intent(in) :: filtertype                     !< Type of filter to initialize
  integer, intent(in) :: dim_p                          !< PE-local state dimension
  integer, intent(in) :: dim_ens                        !< Size of ensemble
  real(pwp), intent(inout) :: state_p(dim_p)            !< PE-local model state
  !< It is not necessary to initialize the array 'state_p' for SEIK. 
  !< It is available here only for convenience and can be used freely.
  real(pwp), intent(inout) :: Uinv(dim_ens-1,dim_ens-1) !< Array not referenced for SEIK
  real(pwp), intent(out)   :: ens_p(dim_p, dim_ens)     !< PE-local state ensemble
  integer, intent(inout) :: flag                        !< PDAF status flag

! *** Local variables ***
  integer :: i, member              ! Counters
  real(pwp) :: ens_mean             ! Ensemble mean value
  real(pwp) :: inv_dim_ens          ! Inverse ensemble size


! ********************************
! *** Read ensemble from files ***
! ********************************

  if (type_ens_init == 0) then

     ! Read ensemble states as model snapshots from a single file

     if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Initialize ensemble from file holding model snapshots'

     call read_ens_mv_loop(path_ens, dim_p, dim_ens, coupling_nemo, ens_p)

  elseif (type_ens_init == 1) then
     if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Initialize ensemble from single ensemble file'

     call read_ens(trim(file_ens), dim_p, dim_ens, ens_p)

  elseif (type_ens_init == 2) then
     
     ! Real ensemble states as model snapshots from separate files

     if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Initialize ensemble from output files'

      call gen_ens_mv(1.0_pwp, .true., path_ens, dim_p, dim_ens, ens_p)

  elseif (type_ens_init == 3) then
     
     ! Real ensemble states as model snapshots from separate files

     if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Initialize ensemble by sampling from covariance matrix'

     CALL gen_ens_from_cov(trim(file_covar), dim_p, dim_ens, state_p, ens_p)

  end if

  ! Options to replace the central state of the ensemble (i.e. the ensemble mean)

  if (type_central_state == 1) then

     ! Read ensemble central state vector state_p
     if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Read central model state of ensemble from file'

     call read_state_mv(path_inistate, dim_p, 1, coupling_nemo, state_p)

  elseif (type_central_state == 2) then

     ! Obtain central state from model task 1 (set by model initialization)
     if (mype_filter==0) write (*,'(a,1x,a)') 'NEMO-PDAF', 'Collect central ensemble state from model'

     call collect_state_pdaf(dim_p, state_p)

  end if


  if (type_central_state==1 .or. type_central_state==2) then

     ! *** Replace ensemble mean and inflate ensemble perturbations ***

     inv_dim_ens = 1.0_pwp/real(dim_ens, kind=pwp)

     if (mype_filter==0) write (*,'(a, 1x,a)') 'NEMO-PDAF', 'Set central state of ensemble'

!$OMP PARALLEL DO private(k, ens_mean)
     do i = 1, dim_p
        ens_mean = 0.0_pwp
        do member = 1, dim_ens
           ens_mean = ens_mean + inv_dim_ens*ens_p(i, member)
        end do
        
        do member = 1, dim_ens
           ens_p(i, member ) = ensscale*(ens_p(i, member) - ens_mean) + state_p(i)
        end do
     end do
!$OMP END PARALLEL DO

  end if

contains
! ===================================================================================

!> Generate ensemble from reading covariance matrix from a file
!!
    subroutine gen_ens_from_cov(filename_cov, dim_p, dim_ens, state_p, ens_p)

      use pdaf_interfaces_module, only: PDAF_SampleEns
      use mod_io_pdaf, only: read_eof_cov
      use mod_parallel_pdaf, &
           only: mype=>mype_filter, abort_parallel


! *** Arguments ***
      character(*), intent(in) :: filename_cov          !< covariance filename
      integer, intent(in)      :: dim_p                 !< dimension of local state vector
      integer, intent(in)      :: dim_ens               !< ensemble size
      real(pwp), intent(inout) :: state_p(dim_p)        !< state vector
      real(pwp), intent(inout) :: ens_p(dim_p, dim_ens) !< ensemble array

! *** Local variables ***
      integer :: rank
      integer :: status_pdaf
      integer :: verbose_sampleens
      real(pwp), allocatable :: eofV(:, :)
      real(pwp), allocatable :: svals(:)
      logical :: readmean

      state_p = 0.0

      ! *****************************************
      ! *** Generate ensemble of model states ***
      ! *****************************************

      ! *** Rank of matrix is ensemble size minus one
      rank = dim_ens - 1

      ! allocate memory for temporary fields
      ALLOCATE(eofV(dim_p, rank))
      ALLOCATE(svals(rank))

      ! get eigenvalue and eigenvectors from file
      readmean = .false.
      call read_eof_cov(filename_cov, dim_state, dim_p, rank, state_p, eofV, svals, readmean)

      ! *** Generate full ensemble on filter-PE 0 ***
      verbose_sampleens = 0
      if (mype==0) then
         WRITE (*, '(a, 1x, a)') 'NEMO-PDAF', '--- generate ensemble using PDAF_SampleEns'
         WRITE (*, '(a, 3x, a)') &
              'NEMO-PDAF', '--- use 2nd order exact sampling'
         WRITE (*, '(a, 3x, a, i5)') 'NEMO-PDAF', '--- Ensemble size:  ', dim_ens
         WRITE (*, '(a, 3x, a, i5)') 'NEMO-PDAF', '--- number of EOFs: ', rank

         if (screen>0) verbose_sampleens = 1
      endif

      ! Use PDAF routine to generate ensemble from covariance matrix
      CALL PDAF_SampleEns(dim_p, dim_ens, eofV, svals, state_p, ens_p, verbose_sampleens, status_pdaf)

      if (status_pdaf /= 0) then
       write (*, '(/1x,a6,i3,a43,i4,a1/)') &
            'ERROR ', status_pdaf, &
            ' in sample ensemble of PDAF - stopping! (PE ', mype, ')'
       call abort_parallel()
      end if


      ! ****************
      ! *** clean up ***
      ! ****************

      DEALLOCATE(svals, eofV)

    end subroutine gen_ens_from_cov

end subroutine init_ens_pdaf
