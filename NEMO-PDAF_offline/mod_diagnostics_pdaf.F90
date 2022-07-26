!> subroutines called by pre/post preocessing subroutine in PDAF
!! Autho
module mod_diagnostics_pdaf
use mod_kind_pdaf
implicit none

logical :: firsttime = .true.

integer :: writestep_var=1       ! Time index for file output of variance
integer :: writestep_state=1     ! Time index for file output of state

real(pwp), allocatable :: ens_p_bkg(:, :)     ! save background ensemble
                                              ! can be used to calculate increments
integer                :: allocflag = 0       ! Flag for memory counting

contains
   ! *** Compute mean state
   subroutine get_meanstate(dim_ens, dim_p, ens_p, state_p)
      integer,   intent(in)     :: dim_ens               !< Size of state ensemble
      integer,   intent(in)     :: dim_p                 !< PE-local state dimension
      real(pwp), intent(inout)  :: state_p(dim_p)        !< PE-local forecast/analysis state
      real(pwp), intent(inout)  :: ens_p(dim_p, dim_ens) !< PE-local state ensemble
      ! *** local variables ***
      integer                   :: i, member             ! counter
      real(pwp)                 :: invdim_ens            ! Inverse ensemble size

      ! Initialize numbers
      invdim_ens    = 1.0_pwp / real(dim_ens, pwp)

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
   end subroutine get_meanstate


   ! *** Compute sampled variances ***
   subroutine get_variance(dim_ens, dim_p, ens_p, state_p, state_tmp)
      integer,   intent(in)     :: dim_ens               !< Size of state ensemble
      integer,   intent(in)     :: dim_p                 !< PE-local state dimension
      real(pwp), intent(inout)  :: state_p(dim_p)        !< PE-local forecast/analysis state
      real(pwp), intent(inout)  :: state_tmp(dim_p)      !< PE-local forecast/analysis state
      real(pwp), intent(inout)  :: ens_p(dim_p, dim_ens) !< PE-local state ensemble
      ! *** local variables ***
      integer                   :: j, member             ! counter
      real(pwp)                 :: invdim_ensm1          ! Inverse of ensemble size minus 1

      ! Initialize numbers
      invdim_ensm1  = 1.0_pwp / real(dim_ens - 1, pwp)

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
   end subroutine get_variance


   ! ************************************************************
   ! *** Compute RMS errors according to sampled covar matrix ***
   ! ************************************************************
   subroutine get_rmse(state_tmp, forana, rmse_est_p, dimfield_p, rmse_est, dimfield)
      use mpi
      use mod_parallel_pdaf, &
           only: mype=>mype_filter, comm_filter, MPIerr
      use mod_statevector_pdaf, &
           only: n_fields, sfields
      ! arguments
      real(pwp), intent(in)    :: state_tmp(:)   ! temporary state vector to hold model state variances or increment
      character(*), intent(in) :: forana         ! String indicating forecast or analysis
      real(pwp), intent(out)   :: rmse_est_p(:)  ! PE-local estimated RMS errors (ensemble standard deviations)
      integer,   intent(out)   :: dimfield_p(:)  ! Local field dimensions

      real(pwp), intent(out)   :: rmse_est(:)    ! Global estimated RMS errors (ensmeble standard deviations)
      integer,   intent(out)   :: dimfield(:)    ! Global field dimensions

      ! local variables
      integer :: i, j ! counter

      rmse_est_p  = 0.0_pwp

      dimfield_p(:) = sfields(:)%dim
      call MPI_Allreduce(dimfield_p, dimfield, n_fields, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)

      ! total estimated mean-square error per field per process
      do j = 1, n_fields
         do i = 1+sfields(j)%off, sfields(j)%dim+sfields(j)%off
            rmse_est_p(j) = rmse_est_p(j) + state_tmp(i)
         enddo
         rmse_est_p(j) = rmse_est_p(j) / real(dimfield(j), pwp)
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
   end subroutine get_rmse


   ! *** Write variance into nc file ***
   subroutine write_variance(ndastp_str, forana, dim_p, state_tmp)
      use mod_parallel_pdaf, &
           only: mype=>mype_filter
      use mod_io_pdaf, &
            only: file_PDAF_variance, &
                  save_var_time, &
                  write_field_mv
      ! arguments
      character(len=8), intent(in)    :: ndastp_str       ! String for model date
      character(len=3), intent(in)    :: forana           ! String indicating forecast or analysis
      integer,          intent(in)    :: dim_p            !< PE-local state dimension
      real(pwp),        intent(inout) :: state_tmp(dim_p) ! temporary state vector to hold model state variances or increment
      ! local variables
      character(len=200) :: titleVar   ! Strings for file titles
      integer            :: nsteps     ! Number of steps written into file

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

      !!!! todo: attime argument is currently always 1.0

         if (forana/='ini') then
            call write_field_mv(state_tmp, trim(file_PDAF_variance)//'_'//trim(ndastp_str)//'.nc', &
                 titleVar, 1.0_pwp, nsteps, writestep_var)
         else
            call write_field_mv(state_tmp, trim(file_PDAF_variance)//'_'//trim(ndastp_str)//'_ini.nc', &
                 titleVar, 1.0_pwp, nsteps, writestep_var)
         end if
         if (forana/='ini') writestep_var = writestep_var + 1

      elseif (writestep_var>1 .and. (trim(save_var_time)=='ana' .or. trim(save_var_time)=='both')) then

         if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write variance after analysis step'

         call write_field_mv(state_tmp, trim(file_PDAF_variance)//'_'//trim(ndastp_str)//'.nc', &
              titleVar, 1.0_pwp, 2, writestep_var)

         writestep_var = 1
      end if
   end subroutine write_variance


   ! *** Write state into nc file ***
   subroutine write_state(ndastp_str, forana, dim_p, state_p, state_tmp)
      use mod_parallel_pdaf, &
           only: mype=>mype_filter
      use mod_io_pdaf, &
            only: file_PDAF_state, &
                  write_field_mv
      ! arguments
      character(len=8), intent(in)    :: ndastp_str       ! String for model date
      character(len=3), intent(in)    :: forana           ! String indicating forecast or analysis
      integer,          intent(in)    :: dim_p            !< PE-local state dimension
      real(pwp),        intent(in) :: state_p(dim_p) ! PE-local state vector
      real(pwp),        intent(inout) :: state_tmp(dim_p) ! temporary state vector to hold model state variances or increment
      ! local variables
      character(len=200)     :: titleState  ! Strings for file titles

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

      !!!! todo: attime argument is currently always 1.0

      ! Write separate files for forecast and analysis
      ! call write_field_mv(state_tmp, dim_p, trim(file_PDAF_state)//'_'//forana//'.nc', titleState, 1.0, 1, 1)
      ! Write forecast and analysis into the same file
      if (forana/='ini') then
         call write_field_mv(state_tmp, &
                             trim(file_PDAF_state)//'_'//trim(ndastp_str)//'.nc', &
                             titleState, 1.0_pwp, 2, writestep_state)
      else
         call write_field_mv(state_tmp, &
                             trim(file_PDAF_state)//'_'//trim(ndastp_str)//'_ini.nc', &
                             titleState, 1.0_pwp, 2, writestep_state)
      end if
   end subroutine write_state


   ! *** Write ensemble and increments into nc file ***
   subroutine write_increments(ndastp_str, dim_p, dim_ens, ens_p)
      use mod_io_pdaf, &
          only: path_restart, &
                update_restart_mv, &
                write_increment, &
                save_incr, &
                save_restart, &
                file_PDAF_incr
      ! arguments
      character(len=8), intent(in)    :: ndastp_str            ! String for model date
      integer,          intent(in)    :: dim_p                 !< PE-local state dimension
      integer,          intent(in)    :: dim_ens               !< ensemble size
      real(pwp),        intent(inout) :: ens_p(dim_p, dim_ens) ! PE-local ensemble
      ! local variables
      integer             :: i     ! counter
      integer             :: ios   ! file status
      character(len=200)  :: path  ! Path of restart file

      ! get increments
      ens_p = ens_p - ens_p_bkg


      if (save_restart) then
         ! read directory for each ensemble member
         open ( unit=10, file=trim(path_restart), iostat=ios)
         if (ios /= 0) write(*,*) 'Could not open file ', trim(path_restart)
         ! write to restart file
         do i = 1, dim_ens
            read(10, '(A)', iostat=ios) path
            if (ios/=0) write(*,*) 'Not sufficient directories for ensemble ', trim(path_restart)
            call update_restart_mv(trim(path), ens_p_bkg(:, i), ens_p(:, i))
         end do
         close (10)
      end if

      if (save_incr) then
         ! write increments
         ! read directory for each ensemble member
         open ( unit=10, file=trim(path_restart), iostat=ios)
         if (ios /= 0) write(*,*) 'Could not open file ', trim(path_restart)
         do i = 1, dim_ens
            read(10, '(A)', iostat=ios) path
            if (ios/=0) write(*,*) 'Not sufficient directories for ensemble ', trim(path_restart)
            call write_increment(ens_p(:, i), &
                                 trim(path)//trim(file_PDAF_incr)//'_'//trim(ndastp_str)//'.nc')
         end do
         close (10)
      end if
   end subroutine write_increments


   subroutine get_rank_histogram(dim_ens, dim_p, state_p, ens_p)
      use pdaf_interfaces_module, &
            only: PDAF_diag_ensstats, &
                  PDAF_diag_histogram
      ! arguments
      integer,   intent(in)     :: dim_ens               !< Size of state ensemble
      integer,   intent(in)     :: dim_p                 !< PE-local state dimension
      real(pwp), intent(inout)  :: state_p(dim_p)        !< PE-local forecast/analysis state
      real(pwp), intent(inout)  :: ens_p(dim_p, dim_ens) !< PE-local state ensemble
      ! local variables
      integer                :: j, k
      integer                :: ios
      integer,   allocatable :: hist_ens(:)     ! ensemble rank histogram
      real(pwp)              :: delta, skewness, kurtosis
      integer                :: status_hist, status_ensstats
      integer                :: elemDiagn
      character(len=20)      :: rankfile

      write(*,*) "*** Ensemble diagnostics post***"
      allocate(hist_ens(dim_ens+1))
      ! if (firsttime) call memcount(3,'r',dim_ens+1)

      ! Read file with elem number for diagnostics
      write(*,*) "*** Reading elem list for diagnostics:"
      write(*,*) 'ElemForEnsStats.txt'

      open (unit=154,file='ElemForEnsStats.txt',iostat=ios)
      if (ios /= 0) write(*,*) 'Could not open file ElemForEnsStats.txt'

      open (157, file='PostEnsDiagnostics.dat',iostat=ios)
       write(157,*)'elem  skewness  kurtosis'

      j=1
      hist_ens=0.0D0
      do
         read (154,*,iostat=ios) elemDiagn
         if (ios/=0) exit
         call PDAF_diag_ensstats(dim_p, dim_ens,elemDiagn,state_p,ens_p,skewness,kurtosis,status_ensstats)
         write(157,*) elemDiagn,skewness,kurtosis
         call PDAF_diag_histogram(j,dim_p,dim_ens,elemDiagn,state_p,ens_p,hist_ens,delta,status_hist)
         j=j+1
      enddo
      close(154)

      rankfile=trim('PostRankHis_Diag.dat')
      print*,'rankfile',rankfile
      open (159, file=rankfile, iostat=ios)
         do k = 1, (dim_ens+1)
            write(159,'(I2,I10)') k, hist_ens(k)
         enddo
      close(159)

      hist_ens = 0.0D0
      call diag_histogram(1,dim_p,dim_ens,0,state_p,ens_p,hist_ens,delta,status_hist)
      call diag_ensstats(dim_p, dim_ens,0,state_p,ens_p,skewness,kurtosis,status_ensstats)

      !Write .dat file with j hist(j) delta skewness kurtosis
      open (158, file='PostRankHistogr.dat', iostat=ios)
         do j = 1, (dim_ens+1)
            write(158,'(I2,I10)') j, hist_ens(j)
         enddo
      close(158)

      write(157,*)'integral',skewness,kurtosis
      close(157)
      deallocate(hist_ens)
   end subroutine get_rank_histogram
end module mod_diagnostics_pdaf