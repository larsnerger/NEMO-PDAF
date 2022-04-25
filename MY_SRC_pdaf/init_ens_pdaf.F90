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
       only:  screen, type_ens_init, type_central_state, ensscale
  use mod_io_pdaf, &
       only: path_inistate, file_inistate_date1, file_inistate_date2, &
       path_ens, file_ens_date1, file_ens_date2, coupling_nemo, &
       read_state_mv, read_ens_mv_loop, read_ens, gen_ens_mv, &
       file_ens, ens_datelist

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
  integer :: i, k, member           ! Counters
  real(pwp) :: ens_mean             ! Ensemble mean value
  real(pwp) :: inv_dim_ens          ! Inverse ensemble size


! ********************************
! *** Read ensemble from files ***
! ********************************

  if (type_ens_init == 0) then

     ! Read ensemble states as model snapshots from a single file

     if (mype_filter==0) write (*,'(1x,a)') 'Initialize ensemble from file holding model snapshots'

     call read_ens_mv_loop(path_ens, file_ens_date1, file_ens_date2, dim_p, dim_ens, &
          coupling_nemo, ens_p)

  elseif (type_ens_init == 1) then
     if (mype_filter==0) write (*,'(1x,a)') 'Initialize ensemble from single ensemble file'

     call read_ens(trim(path_ens)//trim(file_ens), dim_p, dim_ens, ens_p)

  elseif (type_ens_init == 2) then
     
     ! Real ensemble states as model snapshots from separate files

     if (mype_filter==0) write (*,'(1x,a)') 'Initialize ensemble from output files'

     CALL gen_ens_mv(1.0, .true., ens_datelist, path_ens, dim_p, dim_ens, ens_p)

  end if

  ! Options to replace the central state of the ensemble (i.e. the ensemble mean)

  if (type_central_state == 1) then

     ! Read ensemble central state vector state_p
     if (mype_filter==0) write (*,'(1x,a)') 'Read central model state of ensemble from file'

     call read_state_mv(path_inistate, file_inistate_date1, file_inistate_date2, dim_p, 1, &
          coupling_nemo, state_p)

  elseif (type_central_state == 2) then

     ! Obtain central state from model task 1 (set by model initialization)
     if (mype_filter==0) write (*,'(1x,a)') 'Obtain central ensemble state from model'

     call collect_state_pdaf(dim_p, state_p)

  end if


  if (type_central_state==1 .or. type_central_state==2) then

     ! *** Replace ensemble mean and inflate ensemble perturbations ***

     inv_dim_ens = 1.0/real(dim_ens)

     if (mype_filter==0) write (*,'(1x,a)') 'Set central state of ensemble'

!$OMP PARALLEL DO private(k, ens_mean)
     do i = 1, dim_p
        ens_mean = 0.0
        do member = 1, dim_ens
           ens_mean = ens_mean + inv_dim_ens*ens_p(i, member)
        end do
        
        do member = 1, dim_ens
           ens_p(i, member ) = ensscale*(ens_p(i, member) - ens_mean) + state_p(i)
        end do
     end do
!$OMP END PARALLEL DO

  end if

end subroutine init_ens_pdaf
