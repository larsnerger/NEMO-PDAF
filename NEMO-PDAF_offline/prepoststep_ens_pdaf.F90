!> Controlling Pre- and Post-Processing of the PDAF output
!!
!!  - For global filters (e.g. SEIK), the routine is called
!! before the analysis and after the ensemble transformation.
!!  - For local filters (e.g. LSEIK), the routine is called
!! before and after the loop over all local analysis
!! domains.
!! 
!! The routine provides full access to the state
!! estimate and the state ensemble to the user.
!! Thus, user-controlled pre- and poststep
!! operations can be performed here. For example
!! the forecast and the analysis states and ensemble
!! covariance matrix can be analyzed, e.g. by
!! computing the estimated variances.
!! For the offline mode, this routine is the place
!! in which the writing of the analysis ensemble
!! can be performed.
!! 
!! If a user considers to perform adjustments to the
!! estimates (e.g. for balances), this routine is
!! the right place for it.
!! 
!! **Calling Sequence**
!! 
!!  - Called by: `PDAF_get_state` (as U_prepoststep) `PDAF_X_update` (as U_prepoststep)
!! 

subroutine prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
   state_p, Uinv, ens_p, flag)
   use mod_diagnostics_pdaf
   use mod_parallel_pdaf, &
        only: mype=>mype_filter
   use mod_assimilation_pdaf, &
        only: step_null
   use mod_statevector_pdaf, &
        only: n_fields
   use mod_io_pdaf, &
        only: save_state, save_var_time
   use mod_nemo_pdaf, &
        only: ndastp

   implicit none

   ! *** Arguments ***
   integer,   intent(in)    :: step           !< Current time step (negative for call after forecast)
   integer,   intent(in)    :: dim_p          !< PE-local state dimension
   integer,   intent(in)    :: dim_ens        !< Size of state ensemble
   integer,   intent(in)    :: dim_ens_p      !< PE-local size of ensemble
   integer,   intent(in)    :: dim_obs_p      !< Dimension of observation vector
   real(pwp), intent(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
   !< (The array 'state_p' is not generally not initialized in the case of SEIK.
   !< It can be used freely here.)
   real(pwp), intent(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
   real(pwp), intent(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
   integer,   intent(in)    :: flag                       !< PDAF status flag

   ! *** local variables ***
   real(pwp), allocatable :: state_tmp(:)          ! temporary state vector; holds state variances or increment

   real(pwp), allocatable :: rmse_est_p(:)         ! PE-local estimated RMS errors (ensemble standard deviations)
   real(pwp), allocatable :: rmse_est(:)           ! Global estimated RMS errors (ensmeble standard deviations)
   integer,   allocatable :: dimfield_p(:)         ! Local field dimensions
   integer,   allocatable :: dimfield(:)           ! Global field dimensions

   character(len=3)       :: forana                ! String indicating forecast or analysis
   character(len=8)       :: ndastp_str            ! String for model date


   ! **********************
   ! *** INITIALIZATION ***
   ! **********************

   if (step-step_null==0) then
      if (mype==0) write (*,'(a, i7,3x,a)') 'NEMO-PDAF', step, 'Analyze initial state ensemble'
      forana = 'ini'
   else if (step>0) then
      if (mype==0) write (*,'(a, 5x,a)') 'NEMO-PDAF', 'Analyze assimilated state ensemble'
      forana = 'ana'
   else
      if (mype==0) write (*,'(a, 5x,a)') 'NEMO-PDAF', 'Analyze forecast state ensemble'
      forana = 'for'
   end if

   ! Allocate fields
   allocate(state_tmp(dim_p))
   allocate(rmse_est_p(n_fields))
   allocate(rmse_est(n_fields))
   allocate(dimfield_p(n_fields))
   allocate(dimfield(n_fields))
   ! if (firsttime) call memcount(3,'r',dim_p)

   ! **************************************************************
   ! *** Perform prepoststep for ensemble filter.               ***
   ! *** The state and error information is completely in the   ***
   ! *** ensemble.                                              ***
   ! **************************************************************
   call get_meanstate(dim_ens, dim_p, ens_p, state_p)

   call get_variance(dim_ens, dim_p, ens_p, state_p, state_tmp)

   call get_rmse(state_tmp, forana, rmse_est_p, dimfield_p, rmse_est, dimfield)

   ! *******************
   ! *** File output ***
   ! *******************
   ! Set time string
   WRITE(ndastp_str,'(I8.8)') ndastp
   if (trim(save_var_time)=='ana' .or. trim(save_var_time)=='fcst' &
       .or. trim(save_var_time)=='both') then
       call write_variance(ndastp_str, forana, dim_p, state_tmp)
   end if

   if (save_state) then
     call write_state(ndastp_str, forana, dim_p, state_p, state_tmp)
   end if

   ! ********************
   ! *** finishing up ***
   ! ********************

    deallocate(state_tmp) 
    deallocate(rmse_est_p, rmse_est, dimfield_p, dimfield)

    firsttime = .false.
end subroutine prepoststep_ens_pdaf

subroutine prepoststep_ens_pdaf_offline(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
   state_p, Uinv, ens_p, flag)
   use mod_diagnostics_pdaf
   use mod_parallel_pdaf, &
        only: mype=>mype_filter
   use mod_io_pdaf, &
         only: save_state, save_var_time
   use mod_nemo_pdaf, &
        only: ndastp
   use mod_statevector_pdaf, &
        only: n_fields
   implicit none

   ! *** Arguments ***
   integer,   intent(in)     :: step                       !< Current time step (negative for call after forecast)
   integer,   intent(in)     :: dim_p                      !< PE-local state dimension
   integer,   intent(in)     :: dim_ens                    !< Size of state ensemble
   integer,   intent(in)     :: dim_ens_p                  !< PE-local size of ensemble
   integer,   intent(in)     :: dim_obs_p                  !< Dimension of observation vector
   real(pwp), intent(inout)  :: state_p(dim_p)             !< PE-local forecast/analysis state
   !< (The array 'state_p' is not generally not initialized in the case of SEIK.
   !< It can be used freely here.)
   real(pwp), intent(inout)  :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
   real(pwp), intent(inout)  :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
   integer,   intent(in)     :: flag                       !< PDAF status flag

   ! *** local variables ***
   real(pwp), allocatable :: state_tmp(:)    ! temporary state vector to hold model state variances or increment

   real(pwp), allocatable :: rmse_est_p(:)   ! PE-local estimated RMS errors (ensemble standard deviations)
   real(pwp), allocatable :: rmse_est(:)     ! Global estimated RMS errors (ensmeble standard deviations)
   integer,   allocatable :: dimfield_p(:) ! Local field dimensions
   integer,   allocatable :: dimfield(:)  ! Global field dimensions

   character(len=3)     :: forana ! String indicating forecast or analysis
   character(len=8)     :: ndastp_str            ! String for model date

   if (firsttime) then
      if (mype == 0) write (*, '(a, 5x, a)') 'NEMO-PDAF', 'Analyze forecasted state ensemble'
      forana = 'for'
   else
      if (mype == 0)write (*, '(a, 5x, a)') 'NEMO-PDAF', 'Analyze and write assimilated state ensemble'
      forana = 'ana'
   end if

   ! Allocate fields
   allocate(state_tmp(dim_p))
   allocate(rmse_est_p(n_fields))
   allocate(rmse_est(n_fields))
   allocate(dimfield_p(n_fields))
   allocate(dimfield(n_fields))
   ! if (firsttime) call memcount(3,'r',dim_p)

   ! **************************************************************
   ! *** Perform prepoststep for ensemble filter.               ***
   ! *** The state and error information is completely in the   ***
   ! *** ensemble.                                              ***
   ! **************************************************************
   call get_meanstate(dim_ens, dim_p, ens_p, state_p)

   call get_variance(dim_ens, dim_p, ens_p, state_p, state_tmp)

   call get_rmse(state_tmp, forana, rmse_est_p, dimfield_p, rmse_est, dimfield)

   ! *******************
   ! *** File output ***
   ! *******************
   ! Set time string
   WRITE(ndastp_str,'(I8.8)') ndastp
   if (trim(save_var_time)=='ana' .or. trim(save_var_time)=='fcst' &
       .or. trim(save_var_time)=='both') then
       call write_variance(ndastp_str, forana, dim_p, state_tmp)
   end if

   if (save_state) then
      call write_state(ndastp_str, forana, dim_p, state_p, state_tmp)
   end if

   if (firsttime) then
     if (.not. allocated(ens_p_bkg)) allocate(ens_p_bkg(dim_p, dim_ens))
     ens_p_bkg = ens_p
   else
     call write_increments(ndastp_str, dim_p, dim_ens, ens_p)
   end if

   ! ********************
   ! *** finishing up ***
   ! ********************

   deallocate(state_tmp) 
   deallocate(rmse_est_p, rmse_est, dimfield_p, dimfield)

   if ((.not. firsttime) .and. (allocated(ens_p_bkg))) deallocate(ens_p_bkg)
   firsttime = .false.
end subroutine prepoststep_ens_pdaf_offline