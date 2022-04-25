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

  use mpi
  use mod_kind_pdaf
  use mod_assimilation_pdaf, &
        only: step_null
  use mod_parallel_pdaf, &
       only: mype=>mype_filter, comm_filter, MPIerr
  use mod_statevector_pdaf, &
       only: n_fields, id, sfields
  use mod_io_pdaf, &
        only: save_state, save_var_time, file_PDAF_state, file_PDAF_variance, &
        write_field_mv
  use mod_nemo_pdaf, &
       only: ndastp
  
  implicit none

! *** Arguments ***
  integer, intent(in) :: step           !< Current time step (negative for call after forecast)
  integer, intent(in) :: dim_p          !< PE-local state dimension
  integer, intent(in) :: dim_ens        !< Size of state ensemble
  integer, intent(in) :: dim_ens_p      !< PE-local size of ensemble
  integer, intent(in) :: dim_obs_p      !< Dimension of observation vector
  real, intent(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
  !< (The array 'state_p' is not generally not initialized in the case of SEIK.
  !< It can be used freely here.)
  real, intent(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  real, intent(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
  integer, intent(in) :: flag           !< PDAF status flag


! *** local variables ***
  integer :: i, j, member              ! counters
  logical, save :: firsttime = .true.  ! Routine is called for first time?
  real :: invdim_ens                   ! Inverse ensemble size
  real :: invdim_ensm1                 ! Inverse of ensemble size minus 1
  real, allocatable :: rmse_est_p(:)   ! PE-local estimated RMS errors (ensemble standard deviations)
  real, allocatable :: state_tmp(:)    ! temporary state vector; holds state variances or increment
  integer,save :: writestep_var=1      ! Time index for file output of variance
  integer,save :: writestep_state=1    ! Time index for file output of state
  integer :: nsteps                    ! Number of steps written into file
  character(len=3) :: forana           ! String indicating forecast or analysis
  character(len=8) :: ndastp_str       ! String for model date 
  character(len=200) :: titleState, titleVar   ! Strings for file titles
  integer, allocatable :: dimfield_p(:) ! Local field dimensions
  integer, allocatable :: dimfield(:)  ! Global field dimensions
  real, allocatable :: rmse_est(:)     ! Global estimated RMS errors (ensmeble standard deviations)
  


  
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
!  if (firsttime) call memcount(3,'r',dim_p)

  ! Initialize numbers
  invdim_ens    = 1.0_8 / real(dim_ens,8)  
  invdim_ensm1  = 1.0_8 / real(dim_ens - 1,8)



! **************************************************************
! *** Perform prepoststep for ensemble filter.               ***
! *** The state and error information is completely in the   ***
! *** ensemble.                                              ***
! **************************************************************

  ! *** Compute mean state

  state_p = 0.0
  do member = 1, dim_ens
!$OMP PARALLEL DO  
     do i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i, member)
     end do
  end do
!$OMP PARALLEL DO  
  do i = 1, dim_p
    state_p(i) = invdim_ens * state_p(i)
  enddo


  ! *** Compute sampled variances ***
  state_tmp(:) = 0.0
  
  do member = 1, dim_ens
!$OMP PARALLEL DO  
     do j = 1, dim_p
        state_tmp(j) = state_tmp(j) &
             + (ens_p(j, member) - state_p(j)) &
             * (ens_p(j, member) - state_p(j))
     end do
  end do
!$OMP PARALLEL DO  
  do j = 1, dim_p
     state_tmp(j) = invdim_ensm1 * state_tmp(j)
  enddo
  

! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  allocate(rmse_est_p(n_fields))
  allocate(rmse_est(n_fields))
  allocate(dimfield_p(n_fields))
  allocate(dimfield(n_fields))
  rmse_est_p  = 0.0_8

  dimfield_p(:) = sfields(:)%dim
  call MPI_Allreduce(dimfield_p, dimfield, n_fields, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)

  ! total estimated mean-square error per field per process
  do j = 1, n_fields

     do i = 1+sfields(j)%off, sfields(j)%dim+sfields(j)%off
        rmse_est_p(j) = rmse_est_p(j) + state_tmp(i)
     enddo
     rmse_est_p(j) = rmse_est_p(j) / real(dimfield(j), 8)

  enddo

  ! Global sum of mean squared errors
  call MPI_Allreduce (rmse_est_p, rmse_est, n_fields, MPI_DOUBLE_PRECISION, MPI_SUM, &
       COMM_filter, MPIerr)

  ! Get global RMSE
  rmse_est = sqrt(rmse_est)


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
   if (mype == 0) then
      write (*, '(a,6x,a)') 'NEMO-PDAF', 'RMS errors according to sample variance'
      do i = 1, n_fields
         write (*,'(a,4x,a8,4x,a10,2x,es12.4)') &
              'NEMO-PDAF', 'RMSE-'//forana, trim(sfields(i)%variable), rmse_est(i)
      end do
   end if


! *******************
! *** File output ***
! *******************

   ! Set time string
   WRITE(ndastp_str,'(I8.8)') ndastp

   ! *** Write variance into nc file ***

   writevar: if (trim(save_var_time)=='ana' .or. trim(save_var_time)=='fcst' &
        .or. trim(save_var_time)=='both') then

      if (writestep_var==1) then

         if (mype == 0) then
            if (forana=='for') then
               write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write variance before analysis step'
            elseif (forana=='ana') then
               if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write variance after analysis'
            else
               if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write initial variance file'
            end if
         end if

         if (save_var_time=='both' .and. forana/='ini') then
            nsteps = 2
         else
            nsteps = 1
         end if

         titleVar='Ensemble variance'

         if (forana/='ini') then
            call write_field_mv(state_tmp, dim_p, trim(file_PDAF_variance)//'_'//trim(ndastp_str)//'.nc', &
                 titleVar, 1.0, nsteps, writestep_var)
         else
            call write_field_mv(state_tmp, dim_p, trim(file_PDAF_variance)//'_'//trim(ndastp_str)//'_ini.nc', &
                 titleVar, 1.0, nsteps, writestep_var)
         end if
         if (forana/='ini') writestep_var = writestep_var + 1

      elseif (writestep_var>1 .and. (trim(save_var_time)=='ana' .or. trim(save_var_time)=='both')) then

         if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write variance after analysis step'

         call write_field_mv(state_tmp, dim_p, trim(file_PDAF_variance)//'_'//trim(ndastp_str)//'.nc', &
              titleVar, 1.0, 2, writestep_var)

         writestep_var = 1
      end if
   end if writevar


   ! *** Write state into nc file ***
   if (save_State) then

      ! Store state in state_tmp to avoid changing state_p
      state_tmp = state_p

      titleState = 'Ensemble mean state'

      ! Write state file for viewing
      if (forana=='for') then
         if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write ensemble mean before analysis step'
         writestep_state = 1
      elseif (forana=='ana') then
         if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write ensemble mean after analysis step'
         writestep_state = 2
      else
         if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write ensemble mean at initial time'
         writestep_state = 1
      end if

      ! Write separate files for forecast and analysis
!      call write_field_mv(state_tmp, dim_p, trim(file_PDAF_state)//'_'//forana//'.nc', titleState, 1.0, 1, 1)
      ! Write forecast and analysis into the same file
      if (forana/='ini') then
         call write_field_mv(state_tmp, dim_p, trim(file_PDAF_state)//'_'//trim(ndastp_str)//'.nc', titleState, 1.0, 2, writestep_state)
      else
         call write_field_mv(state_tmp, dim_p, trim(file_PDAF_state)//'_'//trim(ndastp_str)//'_ini.nc', titleState, 1.0, 2, writestep_state)
      end if
   endif


! ********************
! *** finishing up ***
! ********************

   deallocate(state_tmp) 
   deallocate(rmse_est_p, rmse_est, dimfield_p, dimfield)

   firsttime = .false.


end subroutine prepoststep_ens_pdaf
