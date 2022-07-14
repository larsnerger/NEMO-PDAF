!$Id$
!BOP
!
! !ROUTINE: init_ens_offline --- Initialize ensemble for SEIK in offline mode
!
! !INTERFACE:
SUBROUTINE init_ens_offline(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.  It has
! to initialize an ensemble of dim\_ens states.
! For the offline mode, the ensemble will be
! typically read-in from files.
!
! The routine is called by all filter processes and 
! initializes the ensemble for the PE-local domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code based on offline_1D
! Later revisions - see svn log
!
! !USES:
  use netcdf
  USE mod_assimilation_pdaf, &
       ONLY: ensfile_type, timeDA, &
       flate, genEnsMeanYearly, nyears, GaussTransf, trafoConst, &
       flateZ, flateTOP, flateBOT, nLevFB, nLevFE
  use mod_io_pdaf, &
       only : read_state_mv, read_ens, gen_ens_mv, write_state_ens, write_ens_files, &
       write_ens_states, write_ens_fields, &
       gen_ensMeanYearly, &
       read_ens_dim_ens_files, read_ens_dims_dim_ens_files, &
       flate_depth, gen_ensFlateZ, nfiles, ntimec, &
       path_state, file_state_date1, file_state_date2, path_ens, file_ens, ens_filelist, &
       coupling_nemo
  USE mod_parallel_pdaf, &
       ONLY: mype=>mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  ! It is not necessary to initialize the array 'state_p' for SEIK. 
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init    (as U_ens_init)
!EOP

! *** local variables ***
  INTEGER :: i, member           ! Counters
  REAL :: invdim_ensm1                 ! Inverse of ensemble size minus 1
  CHARACTER(len=200):: titleVar
  REAL :: incrTime,startEnsTime,endEnsTime
  INTEGER, ALLOCATABLE :: hist_ens(:)
  REAL              :: delta,skewness,kurtosis
  INTEGER           :: status_hist,status_ensstats 
  INTEGER           :: ios,elemDiagn,k
  CHARACTER(len=20) :: rankfile
  REAL, ALLOCATABLE :: flate_z(:)

  integer(4) :: dim_state_ergom
  integer(4) :: dim_state_ensfile
  integer(4) :: dim_ens_ensfile


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype==0) THEN
     WRITE (*, '(/9x, a)') 'Initialize state ensemble'
     WRITE (*, '(9x, a)') '--- read ensemble from files'
     WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
  END IF


! ***************************
! *** Check ensemble size ***
! ***************************
  
  ntimec = 1

!LN TODO: Move this to the actual ensemble reading
if (1==2) then
  if (ensfile_type == 1) then
     ! Read ensemble dimension from ensemble file
  elseif (ensfile_type == 2) then
     ! Read ensemble dimension from a list of files
  elseif (ensfile_type == 3) then
     ! Determine ensemble size from number of files present
     call read_ens_dims_dim_ens_files(path_ens,dim_state_ensfile,dim_ens_ensfile)
  else
     write (*,*) 'ERROR: No valid choice of ensfile_type', ensfile_type
  endif
endif

!   if (dim_ens /= dim_ens_ensfile) then
!      write (*,'(/1x,a)') 'dim_ens in pdaf.nml and from ensemble input file are not the same:', &
!           dim_ens, dim_ens_ensfile
!      write (*,'(/1x,a)')  'ERGOM-PDAF is gooing to exit.' 
!      call exit(0)  
!   endif


! ********************************
! *** Read ensemble from files ***
! ********************************

  ! Read ensemble central state vector state_p
  write (*,'(1x,a)') 'Read central model state of the ensemble'

  IF (coupling_nemo=='incr') THEN
     CALL read_state_mv(path_state, file_state_date1, file_state_date2, dim_p, timeDA, &
          coupling_nemo, state_p)
  ELSE
     CALL read_state_mv(path_state, file_state_date1, file_state_date2, dim_p, 1, &
          coupling_nemo, state_p)
  END IF

  IF (flateZ) THEN
     ALLOCATE (flate_z(dim_p))
     CALL flate_depth(dim_p,nLevFB,nLevFE,flateTOP,flateBOT,flate_z)
  ENDIF

  ! Read ensemble perturbations
  if ( ensfile_type == 1 ) then
     write (*,'(1x,a)') 'Initialize ensemble from single ensemble file'

     call read_ens(trim(path_ens)//trim(file_ens), dim_p, dim_ens, ens_p)

  else if  ( ensfile_type == 2 ) then

     IF (.NOT. genEnsMeanYearly) THEN
        IF (.NOT. flateZ) THEN
           if (mype==0) write (*,'(1x,a)') 'Initialize ensemble from NEMO/ERGOM output files'

           CALL gen_ens_mv(flate, ens_filelist, path_state, dim_p, dim_ens, ens_p)
        ELSE
!            call gen_ensFlateZ(flate_z,varname,ens_filelist,path_ens,dim_p,dim_ens,GaussTransf,trafoConst,ens_p)
        ENDIF
     ELSE
!        call gen_ensMeanYearly(flate,varname,ens_filelist,path_ens,dim_p,dim_ens,nyears,GaussTransf,trafoConst,ens_p)!TO DO: domain decomp parallelisation (dim_p=dim_state at the moment)
     ENDIF
     IF (flateZ) DEALLOCATE (flate_z)
     IF (write_ens_states) call write_state_ens(path_ens, file_ens, dim_p, dim_ens, ens_p)

  else if ( ensfile_type == 3) then
!     call read_ens_dim_ens_files(path_ens, dim_p, dim_ens, ens_p)
  endif

  ! write ensemble of files holding ensemble of fields
  ! Careful: This routine applied the backward state transformation
  IF (write_ens_fields) call write_ens_files(path_ens,'ensembleField',dim_p,dim_ens,ens_p)


  ! Add ensemble central state and perturbations
  DO member = 1, dim_ens
!$OMP PARALLEL DO  
     DO i=1, dim_p
        ens_p(i, member) = ens_p(i, member) + state_p(i)
     END DO
!$OMP END PARALLEL DO
  END DO


! ****************
! *** clean up ***
! ****************

END SUBROUTINE init_ens_offline
