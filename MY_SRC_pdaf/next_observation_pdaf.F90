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
  use mod_nemo_pdaf, &
       only: nitend, nit000
  use mod_iau_pdaf, &
       only: store_asm_step_pdaf, update_asm_step_pdaf


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
  doexit = 0

  if (stepnow + delt_obs <= nitend) then
     ! *** During the assimilation process ***
     
     if (stepnow == nit000 - 1) then
        ! First analysis step 
        nsteps = delt_obs-1     ! Analysis step one step before end of day
     else
        nsteps = delt_obs       ! Follow-up analysis steps daily
     end if

     if (mype_ens == 0) write (*, '(a, i7, 3x, a, i7)') &
          'NEMO-PDAF', stepnow, 'Next observation at time step', stepnow + nsteps

     ! Update analysis step information for NEMO-ASM
     if (stepnow == nit000 - 1) then
        ! At initial time - apply increments after first analysis step
        call store_asm_step_pdaf(stepnow+nsteps)
     else
        ! analysis step - apply increments after current analysis step
        call store_asm_step_pdaf(stepnow)
     end if

  else
     ! *** End of assimilation process ***

     ! Set nsteps so stepnow+delt_bs>nitend to ensure that PDAF calls distribute_state
     nsteps = delt_obs        

     if (mype_ens == 0) write (*, '(a, i7, 3x, a)') &
          'NEMO-PDAF', stepnow, 'No more observations - end assimilation and complete forecast'

     ! Update analysis step information for NEMO-ASM
     call store_asm_step_pdaf(stepnow)
  end if

end subroutine next_observation_pdaf
