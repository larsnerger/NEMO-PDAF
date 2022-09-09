!> Determining the Next Analysis Step
!!
!! The subroutine is called before each forecast phase
!! by `PDAF_get_state`. It has to initialize the number
!! of time steps until the next available observation
!! (`nsteps`). It indicates if the data assimilation process
!! is completed such that the ensemble loop in the model
!! routine can be exited.
!! 
!! The routine is called by all processes.
!! 
!!  **Calling Sequence**
!! 
!!  - Called from: `init_pdaf/PDAF_get_state` (as U_next_obs)
!! 
subroutine next_observation_pdaf(stepnow, nsteps, doexit, time)

  use mod_kind_pdaf
  use mod_assimilation_pdaf, &
       only: delt_obs
  use mod_parallel_pdaf, &
       only: mype_ens
  use in_out_manager, &
       only: nitend

  implicit none

! *** Arguments ***   
  integer, intent(in)  :: stepnow !< Number of the current time step
  integer, intent(out) :: nsteps  !< Number of time steps until next obs
  integer, intent(out) :: doexit  !< Whether to exit forecasting (1 for exit)
  real(pwp), intent(out) :: time  !< Current model (physical) time


! *******************************************************
! *** Set number of time steps until next observation ***
! *******************************************************

  ! Not used in this implementation
  time = 0.0

  if (stepnow + delt_obs <= nitend) then
     ! *** During the assimilation process ***
     nsteps = delt_obs   ! This assumes a constant time step interval
     doexit = 0          ! Not used in this implementation

     if (mype_ens == 0) write (*, '(a, i7, 3x, a, i7)') &
          'NEMO-PDAF', stepnow, 'Next observation at time step', stepnow + nsteps
  else
     ! *** End of assimilation process ***
     nsteps = 0          ! No more steps
     doexit = 1          ! Not used in this implementation

     if (mype_ens == 0) write (*, '(a, i7, 3x, a)') &
          'NEMO-PDAF', stepnow, 'No more observations - end assimilation'
  end if

end subroutine next_observation_pdaf
