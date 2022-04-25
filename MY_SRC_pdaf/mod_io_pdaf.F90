!> Module holding IO operations for NEMO-PDAF
!!
module mod_io_pdaf

  use mpi

  ! Include dimension information for model grid
  use mod_nemo_pdaf, &
       only: nlvls=>jpk, nlats=>jpjglo, nlons=>jpiglo, &
       depths=>gdept_1d, lons=>glamt, lats=>gphit, &
       tmp_4d, ni_p, nj_p, nk_p, istart, jstart, &
       nimpp, njmpp, nlei, nlej

  ! Include information on state vector
  use mod_statevector_pdaf, &
       only: id, sfields, n_fields

  ! Include parallelization information
  use mod_parallel_pdaf, &
       only: mype=>mype_ens, npes=>npes_ens, comm_filter

  ! Include auxiliary routines
  use mod_aux_pdaf, &
       only: field2state, state2field, transform_field, transform_field_mv

  implicit none
  save

  integer :: verbose=1   ! Set verbosity of IO routines (0,1,2,3)

  ! Control of IO
  character(len=4) :: save_var_time='both'   ! Write variance at 'fcst', 'ana', 'both', or 'none'
  logical :: save_ens_states=.false.         ! Write a single file of ensmeble state vectors
  logical :: save_ens_fields=.false.         ! Write set of files holding ensemble fields
  logical :: save_state=.true.               ! Write analysis state to file
  logical :: save_incr                       ! Write increment to file
  logical :: do_deflate=.false.              ! Deflate variables in NC files (this seems to fail for parallel nc)
  character(len=4)   :: coupling_nemo = 'odir'   ! offline: 'rest', 'incr', online: 'oinc', 'odir'

  character(len=100) :: file_PDAF_state='state'       ! File name for outputs of ensemble mean state
  character(len=100) :: file_PDAF_incr='incr'         ! File name for increment
  character(len=100) :: file_PDAF_variance='variance' ! File name for ensemble variance
  character(len=200) :: path_inistate      ! Path to NEMO files
  character(len=8)   :: file_ens_date1     ! Date 1 in NEMO file name
  character(len=8)   :: file_ens_date2     ! Date 2 in NEMO file name
  character(len=8)   :: file_inistate_date1     ! Date 1 in NEMO file name
  character(len=8)   :: file_inistate_date2     ! Date 2 in NEMO file name
  character(len=200) :: path_ens           ! Path of ensemble file  
  character(len=80)  :: file_ens           ! File name of ensemble file
  character(len=200) :: path_restart       ! Path of restart file
  character(len=80)  :: file_restart       ! file name of restart dile
  character(len=80)  :: ens_datelist       ! Name of file holding dates to read in ensembles states
  character(len=1)   :: datestype='_'      ! Linking character in dates (`_` or `-`)

   ! Temporary - from offline code
  real :: startEnsTime=1.0, endEnsTime=1.0, incrTime=1.0

  ! NEMO output file
  integer(4)        :: ntimec=1

  ! Missing value in netcdf file
  real(8) :: missing_value
    
contains
! ===================================================================================

!> Read fields from NEMO file into a state vector
!!
  subroutine read_state_mv(path, date1, date2, dim_p, itime, coupling, state)

    use netcdf
    
    implicit none
    
    character(len = *), intent(in)    :: path          !< Path to file
    character(len = *), intent(in)    :: date1         !< Start date in file name
    character(len = *), intent(in)    :: date2         !< End date in file name
    integer(4),         intent(in)    :: dim_p         !< PE-local state dimension
    integer(4),         intent(in)    :: itime         !< Time to read in file
    character(len = *), intent(in)    :: coupling      !< Type of NEMO coupling
    real(8),            intent(inout) :: state(dim_p)  !< State vector
    
    ! Local variables 
    integer(4) :: i, cnt          ! Counters
    integer(4) :: varid           ! Variable ID
    integer(4) :: ncid            ! NC file id
    character(len=50) :: filename ! Full file name
    character(len=17) :: dates    ! String with initial and final date

    if (verbose>0 .and. mype==0) &
         write(*,'(a,4x,a,i8)') 'NEMO-PDAF', '*** Read model output at time step: ', itime

    if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))

    ! Initialize state
    state = 0.0

    if (datestype=='_') then
       dates = trim(date1)//'_'//trim(date2)
    else
       dates = trim(date1)//'-'//trim(date2)
    end if

    do i = 1, n_fields
    
       filename = trim(sfields(i)%file)//trim(dates)//trim(sfields(i)%file_post)//'.nc'
       if (verbose>1 .and. mype==0) then 
          write(*,'(a,2x,a)') 'NEMO-PDAF', trim(path)//trim(filename)
          write (*,'(a,i5,a,a,a,i10)') &
               'NEMO-PDAF', i, 'Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
       end if

       ! Open the file
       call check( nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid) )

       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), varid) )
       if (coupling/='rest') then
          call check( nf90_get_att(ncid, varid, 'missing_value', missing_value) )
       else
          missing_value=0.0
       endif

       ! Read variable
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), varid) )

       if (sfields(i)%ndims == 3) then
          call check( nf90_get_var(ncid, varid, tmp_4d, &
               start=(/istart, jstart, 1, itime/), count=(/ni_p, nj_p, nlvls, 1/)) )
       else
          call check( nf90_get_var(ncid, varid, tmp_4d(:,:,1,1), &
               start=(/istart, jstart, itime/), count=(/ni_p, nj_p, 1/)) )
       end if

       call check( nf90_close(ncid) )

       ! Convert field to state vector
       call field2state(tmp_4d, state, sfields(i)%off, sfields(i)%ndims, missing_value)

    end do
   
    ! Potentially transform fields
    call transform_field_mv(1, state)

    if (verbose>2) then
       do i = 1, n_fields
          write(*,'(a, 1x, a, a10, 1x, a,5x, 2f12.6)') &
               'NEMO-PDAF', 'Min and max for ',trim(sfields(i)%variable),' :     ',              &
               minval(state(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim)), &
               maxval(state(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim))
       enddo
    end if

  end subroutine read_state_mv


!=============================================================================== 

!> Read an ensemble of model fields into the ensemble array
!!  
  subroutine read_ens_mv_loop(path, date1, date2, dim_p, dim_ens, coupling, ens)

    use netcdf
    
    implicit none

! *** Arguments ***    
    character(len = *), intent(in)   :: path      !< Path of file
    character(len = *), intent(in)   :: date1     !< Start date in file name
    character(len = *), intent(in)   :: date2     !< End date in file name
    integer(4),         intent(in)   :: dim_p     !< State dimension
    integer(4),         intent(in)   :: dim_ens   !< Ensemble size
    character(len = *), intent(in)   :: coupling  !< Type of NEMO coupling
    real(8),            intent(inout):: ens(:,:)  !< Ensemble array
    
! *** Local variables ***
    integer(4) :: i, cnt, member   ! Counters
    integer(4) :: ncid             ! NC file ID
    integer(4) :: varid            ! Variable ID
    character(len=50) :: filename  ! Full file name
    character(len=17) :: dates     ! Combined date string of file

    if (verbose>0 .and. mype==0) &
         write(*,'(a,4x,a)') 'NEMO-PDAF','*** Ensemble: Read model snapshots'

    if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))

    ! Initialize ensemble
    ens = 0.0

    do i = 1, n_fields

       if (datestype=='_') then
          dates = trim(date1)//'_'//trim(date2)
       else
          dates = trim(date1)//'-'//trim(date2)
       end if
       filename = trim(sfields(i)%file)//trim(dates)//trim(sfields(i)%file_post)//'.nc'
       if (verbose>1 .and. mype==0) then
          write(*,'(a,2x,a)') 'NEMO-PDAF', trim(path)//trim(filename)
          write (*,'(a,i5,a,a,a,i10)') &
               'NEMO-PDAF', i, 'Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
       end if

       ! Open the file
       call check( nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid) )

       !  Read field
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), varid) )

       ! Read missing value
       if (coupling/='rest') then
          call check( nf90_get_att(ncid, varid, 'missing_value', missing_value) )
       else
          missing_value=0.0
       endif

       do member = 1, dim_ens

          if (verbose>0 .and. mype==0 .and. i==1) &
               write (*,'(a,4x,a,i6)') 'NEMO-PDAF','--- read member', member

          if (sfields(i)%ndims == 3) then
             call check( nf90_get_var(ncid, varid, tmp_4d, &
                  start=(/istart, jstart, 1, member/), count=(/ni_p, nj_p, nlvls, 1/)) )
          else
             call check( nf90_get_var(ncid, varid, tmp_4d(:,:,1,1), &
                  start=(/istart, jstart, member/), count=(/ni_p, nj_p, 1/)) )
          end if

          ! Convert field to state vector
          call field2state(tmp_4d, ens(:,member), sfields(i)%off, sfields(i)%ndims, missing_value)

       enddo

       call check( nf90_close(ncid) )

    end do

    do member = 1, dim_ens
       ! Potentially transform fields
       call transform_field_mv(1, ens(:,member))
    end do

    if (verbose>2) then
       do i = 1, n_fields
          write(*,'(a, 1x, a, a10, 1x, a,1x, 2es13.6)') &
               'NEMO-PDAF','Ensemble min and max for ',trim(sfields(i)%variable),' :     ', &
               minval(ens(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim,:)), &
               maxval(ens(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim,:))
       enddo
    end if

  end subroutine read_ens_mv_loop


!===============================================================================  

!!> Read ensemble as state vectors from ensemble file
!!  
  subroutine read_ens(ensfile_fullname, dim_state, dim_ens, ens)

    use netcdf
    
    implicit none
    
    character(len=*), intent(in)    :: ensfile_fullname        !< Name and path of ensemble file
    integer(4),       intent(in)    :: dim_state               !< PE-local state dimension
    integer(4),       intent(in)    :: dim_ens                 !< Ensemble size
    real(8),          intent(inout) :: ens(dim_state, dim_ens) !< Ensemble array
    
    ! Local variables  
    integer(4) :: i                     ! Counter
    integer(4) :: ncid                  ! NC file id
    integer(4) :: dim_state_file        ! state dimension in file
    integer(4) :: dim_ens_file          ! Ensemble size in file
    character(len=400) :: varstr        ! String describing variables in state vector
    character(len=400) :: varstr_file   ! String describing variables in state vector
    integer(4) :: dimstate_dimid, dimens_dimid, ens_varid


! *** Generate string describing the state vector ***
    varstr = ''
    do i = 1, n_fields
       if (i==1) then
          varstr = trim(sfields(i)%variable)
       else
          varstr = trim(varstr)//' '//trim(sfields(i)%variable)
       endif
    end do

! *** Read file
    
    if (verbose>0 .and. mype==0) then
       write(*,'(1x,a,a)') "--- Read ensemble file: ", trim(ensfile_fullname)
    end if
    
    ! Open the file
    call check( nf90_open(ensfile_fullname, nf90_nowrite, ncid) )

    ! Read the string describing the state vector
    call check( nf90_get_att(ncid, NF90_GLOBAL, "state_fields", varstr_file) )

    ! Check consistency of state vector setup
    if (trim(varstr) == trim(varstr_file)) then
    
       ! Get the dimensions
       call check( nf90_inq_dimid(ncid, 'dim_state', dimstate_dimid) )  
       call check(nf90_inquire_dimension(ncid,dimstate_dimid,len=dim_state_file))
  
       call check( nf90_inq_dimid(ncid, 'dim_ens', dimens_dimid) )  
       call check(nf90_inquire_dimension(ncid,dimens_dimid,len=dim_ens_file))

       ! Check consistency of state dimension
       if (dim_state_file == dim_state) then

          ! Check consistency of ensemble size
          if (dim_ens_file >= dim_ens) then

             !  Read ensemble
             call check( nf90_inq_varid(ncid, 'ensemble', ens_varid) )
 
             call check( nf90_get_var(ncid, ens_varid, ens, start=(/1,1/),count=(/dim_state,dim_ens/)) )

             if (dim_ens_file> dim_ens) &
                  write (*,*) 'Notice: Ensemble in file is larger than dim_ens'
          else
             write (*,'(1x,a)') 'ERROR: Ensemble in file is too small'
             write (*,'(1x,a)')  'Stopping program!' 
             stop 10
          end if

       else
          write (*,'(1x,a)') 'ERROR: inconsistent state dimension'
          write (*,'(1x,a)')  'Stopping program!' 
          stop 10
       end if

    else
       write (*,'(1x,a)') 'ERROR: inconsistent variables in state'
       write (*,'(1x,a)')  'Stopping program!' 
       stop 10
    end if

    call check( nf90_close(ncid) )

  end subroutine read_ens


!================================================================================

!> Initialie ensemble array from a list of NEMO output files
!!
  subroutine gen_ens_mv(flate, infilelist, inpath, dim_p, dim_ens, ens)

  implicit none
  
! *** Arguments ***
  real(8),            intent(in)    :: flate        !< inflation
  character(len=*),   intent(in)    :: infilelist   !< Name of file holding dates of input files
  character(len=*),   intent(in)    :: inpath       !< Path to input files
  integer(4),         intent(in)    :: dim_p        !< State dimension
  integer(4),         intent(in)    :: dim_ens      !< Ensemble size
  real(8),            intent(inout) :: ens(:, :)    !< Ensemble array

! *** Local variables ***
  integer(4)        :: i, k, iens, ifile  ! Counters
  integer           :: ios                ! Flag for file reading
  character(len=8)  :: indate             ! Date string
  real(8)           :: ens_mean           ! Ensemble mean
  real(8)           :: invsteps           ! Inverse of ensemble size


! *** Read ensemble from files ***

  if (verbose>0 .and. mype==0) &
       write(*,'(/1x,a)') "*** Generating ensemble from output files ***"

  open (unit=10,file=trim(infilelist),iostat=ios)
  if (ios /= 0) write(*,*) 'Could not open file ',infilelist 
  
  iens=0

  ensloop: do

     read (10,*,iostat=ios) indate
     if (ios/=0) exit ensloop

     do k =1, ntimec
        iens = iens + 1
        if (verbose>0 .and. mype==0) write (*,*) '--- Read ensemble member', iens
        call read_ens_mv(inpath, indate, indate, dim_p, k, ens(:,iens))
     enddo

     if (iens==dim_ens) exit ensloop

  enddo ensloop
  
  close(10)

  ! Check ensemble size
  if (iens<dim_ens) then
     write (*,'(/1x,a)') 'ERROR: Available files less than ensemble size!'
     write (*,'(1x,a)')  'Stopping program!' 
     stop 10
  end if

 
! *** Subtract ensemble mean and inflate ensemble perturbations ***

  invsteps = 1.0/real(dim_ens)


!$OMP PARALLEL DO private(k, ens_mean)
  do k=1,dim_p
     ens_mean = 0.0
     do i=1,dim_ens
        ens_mean = ens_mean + invsteps*ens(k,i)
     end do

     do i=1,dim_ens
        ens(k,i) = flate*(ens(k,i)-ens_mean)
     end do
  end do
!$OMP END PARALLEL DO


end subroutine gen_ens_mv


!=============================================================================== 

!> Read a model field into the state vector of an ensemble array
!!  
  subroutine read_ens_mv(path, date1, date2, dim_p, itime, state)

    use netcdf
    
    implicit none

! *** Arguments ***    
    character(len = *), intent(in)   :: path      !< Path of file
    character(len = *), intent(in)   :: date1     !< First date in file name
    character(len = *), intent(in)   :: date2     !< Second date in file name
    integer(4),         intent(in)   :: dim_p     !< State dimension
    integer(4),         intent(in)   :: itime     !< Time in file to read
    real(8),            intent(inout):: state(:)  !< State vector
    
! *** Local variables ***
    integer(4) :: i, cnt           ! Counters
    integer(4) :: ncid             ! NC file ID
    integer(4) :: varid            ! Variable ID
    character(len=50) :: filename  ! Full file name
    character(len=17) :: dates     ! Combined date string of file

    if (verbose>1) &
         write(*,*) "*** Ensemble: Read model output at time step: ", itime

    ! Initialize state
    state = 0.0

    do i = 1, n_fields

       if (datestype=='_') then
          dates = trim(date1)//'_'//trim(date2)
       else
          dates = trim(date1)//'-'//trim(date2)
       end if
       filename = trim(sfields(i)%file)//trim(dates)//trim(sfields(i)%file_post)//'.nc'

       if (verbose>1) then
          write(*,*) trim(path)//trim(filename)
          write (*,*) i, 'Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
       end if

       ! Open the file
       call check( nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid) )

       !  Read field
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), varid) )

       if (sfields(i)%ndims == 3) then
          call check( nf90_get_var(ncid, varid, tmp_4d, &
               start=(/istart, jstart, 1, itime/), count=(/ni_p, nj_p, nlvls, 1/)) )
       else
          call check( nf90_get_var(ncid, varid, tmp_4d(:,:,1,1), &
               start=(/istart, jstart, itime/), count=(/ni_p, nj_p, 1/)) )
       end if

       call check( nf90_close(ncid) )


       ! Convert field to state vector
       call field2state(tmp_4d, state, sfields(i)%off, sfields(i)%ndims, missing_value)

    end do
 
    ! Potentially transform fields
    call transform_field_mv(1, state)

    if (verbose>1) then
       do i = 1, n_fields
          write(*,*) 'Min and max for ',trim(sfields(i)%variable),' :     ',              &
               minval(state(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim)), &
               maxval(state(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim))
       enddo
    end if

  end subroutine read_ens_mv

!================================================================================

!> Write an ensemble file holding the state vectors
!!
  subroutine write_state_ens(path, file, dim_state, dim_ens, ens)

    use netcdf

    implicit none 

! *** Arguments *** 
    character(len=*), intent(in):: path          !< Path of file
    character(len=*), intent(in):: file          !< File name
    integer(4),       intent(in):: dim_state     !< state dimension
    integer(4),       intent(in):: dim_ens       !< Ensemble size
    real(8),          intent(in):: ens(:,:)      !< Ensemble array

! *** Local variables *** 
    integer(4) :: i          ! Counter
    integer(4) :: fileid     ! NC file id
    integer(4) :: dimids(2)  ! dimension ids
    integer(4) :: id_ens     ! variable id
    real(8)    :: fillval    ! fill value
    integer(4) :: startv(2),countv(2)  ! Arrays for writing
    character(len=400) :: varstr   ! String describing variables in state vector
    character(len=200) :: filestr  ! String for file name

! *** Generate string describing the state vector ***
    varstr = ''
    do i = 1, n_fields
       if (i==1) then
          varstr = trim(sfields(i)%variable)
       else
          varstr = trim(varstr)//' '//trim(sfields(i)%variable)
       endif
    end do

    if (npes==1) then
       filestr = trim(file)//'.nc'
    else
       filestr = trim(file)//'_'//trim(str(mype))//'.nc'
    end if

! *** Write ensemble of state vectors ***

    ! *** Open file and initialize dimensions and fields *** 
    call check( NF90_CREATE(trim(filestr),NF90_NETCDF4,fileid) )
    call check( NF90_PUT_ATT(fileid,NF90_GLOBAL,'title', &
         'Ensemble matrix for NEMO') )
    call check( nf90_put_att(fileid, NF90_GLOBAL, "state_fields", trim(varstr)) )
  
    ! define dimensions
    call check( NF90_DEF_DIM(fileid,'dim_state',dim_state,dimids(1)) )
    call check( NF90_DEF_DIM(fileid,'dim_ens',dim_ens,dimids(2)) )

    ! define variables
    call check( NF90_DEF_VAR(fileid,'ensemble',NF90_DOUBLE,dimids(1:2),id_ens) )
    fillval = 0.0
    call check( nf90_put_att(fileid, id_ens, "_FillValue", fillval) )
    call check( nf90_put_att(fileid, id_ens, "missing_value", fillval) )
    call check( NF90_def_var_deflate(fileid,id_ens,0,1,1) )

    ! End define mode
    call check( NF90_ENDDEF(fileid) )
 
    do i=1,dim_ens
       startv(1) = 1
       startv(2) = i
       countv(1) = dim_state
       countv(2) = 1
       call check( nf90_put_var(fileid,id_ens,ens(1:dim_state,i), startv,countv) )
    end do

    ! *** close file with state sequence ***
    call check( NF90_CLOSE(fileid) )

  end subroutine write_state_ens


!================================================================================

!> Write ensemble as single files holding model fields
!!
  subroutine write_ens_files(path,file_ens,dim_state,dim_ens,ens)

    use netcdf

    implicit none 
 
! *** Arguments ***
    character(len=*), intent(in) :: path        !< Path of file
    character(len=*), intent(in) :: file_ens    !< Name stub of file
    integer(4),       intent(in) :: dim_state   !< State dimension
    integer(4),       intent(in) :: dim_ens     !< Ensemble size
    real(8),       intent(inout) :: ens(:, :)   !< Ensemble array

! *** Local variables ***
    integer(4)          :: i                ! Counter
    character(len=200)  :: file_ensemble   ! Full name of an ensemble size
    character(len=200)  :: titleEns         ! NC title of file
    real(8)             :: time             ! Time in file


! *** Write ensemble perturbation files

    time=10.0 !TO DO: this is random, time has to be read in pdaf.nml and set here

    titleEns='Ensemble perturbation (ens-mean) for PDAF'
    
    do i=1,dim_ens

       file_ensemble=trim(path)//trim(file_ens)//'_'//trim(str(i))//'.nc'

       call write_field_mv(ens(:, i), dim_state, file_ensemble, titleEns,time, 1, 1)
    enddo

  end subroutine write_ens_files

!================================================================================

!> Write a state vector as model fields into a file
!!
  subroutine write_field_mv(state, dim, filename, title, &
       attime, nsteps, step)
    
    use netcdf
   
    implicit none

! *** Arguments ***
    real(8),          intent(inout) :: state(:)  ! State vector
    integer(4),       intent(in) :: dim          ! State dimension
    character(len=*), intent(in) :: filename     ! File name
    character(len=*), intent(in) :: title        ! File title
    real(8),          intent(in) :: attime       ! Time attribute
    integer(4),       intent(in) :: nsteps       ! Number of time steps stored in file
    integer(4),       intent(in) :: step         ! Time index to write at

! *** Local variables ***
    integer(4) :: ncid
    integer(4) :: dimids_field(4)
    integer(4) :: cnt,i,j,k
    integer(4) :: dimid_time, dimid_lvls, dimid_lat, dimid_lon, dimid_one
    integer(4) :: id_lat, id_lon, id_lev, id_time, id_field
    integer(4) :: dimids(4)
    integer(4) :: startC(2), countC(2)
    integer(4) :: startt(4), countt(4)
    integer(4) :: startz(1), countz(1)
    real(8)    :: fillval
    real(8)    :: timeField(1)

    timeField(1)=attime

    if (step==1) then

! *** Create file ***

       if (verbose>0 .and. mype==0) &
            write (*,'(a,1x,a,a)') 'NEMO-PDAF', 'Create file: ', trim(filename)

       if (npes==1) then
          call check( NF90_CREATE(trim(filename),NF90_NETCDF4,ncid))
       else
          call check( NF90_CREATE_PAR(trim(filename), NF90_NETCDF4, comm_filter, MPI_INFO_NULL, ncid))
       end if
       call check( NF90_PUT_ATT(ncid, NF90_GLOBAL, 'title', trim(title)))
     
       ! define dimensions for NEMO-input file
       call check( NF90_DEF_DIM(ncid,'t', nsteps, dimid_time))
       call check( NF90_DEF_DIM(ncid, 'z', nlvls, dimid_lvls))
       call check( NF90_DEF_DIM(ncid, 'y', nlats, dimid_lat) )
       call check( NF90_DEF_DIM(ncid, 'x', nlons, dimid_lon) )
       call check( NF90_DEF_DIM(ncid, 'one', 1, dimid_one) )
   
       dimids_field(4)=dimid_time
       dimids_field(3)=dimid_lvls
       dimids_field(2)=dimid_lat
       dimids_field(1)=dimid_lon
       
       ! define variables
       call check( NF90_DEF_VAR(ncid, 'time', NF90_DOUBLE, id_time))
       call check( NF90_DEF_VAR(ncid, 'nav_lat', NF90_FLOAT, dimids_field(1:2), id_lat))
       call check( NF90_DEF_VAR(ncid, 'nav_lon', NF90_FLOAT, dimids_field(1:2), id_lon))
       call check( NF90_DEF_VAR(ncid, 'nav_lev', NF90_FLOAT, dimids_field(3), id_lev))
       if (do_deflate) then
          call check( NF90_def_var_deflate(ncid, id_lat, 0, 1, 1) )
          call check( NF90_def_var_deflate(ncid, id_lon, 0, 1, 1) )
          call check( NF90_def_var_deflate(ncid, id_lev, 0, 1, 1) )
       end if

       do i = 1, n_fields
          if (sfields(i)%ndims==3) then
             dimids_field(3)=dimid_lvls
             call check( NF90_DEF_VAR(ncid, trim(sfields(i)%variable), NF90_DOUBLE, dimids_field(1:4), id_field) )
          else
             dimids_field(3)=dimid_time
             call check( NF90_DEF_VAR(ncid, trim(sfields(i)%variable), NF90_DOUBLE, dimids_field(1:3), id_field) )
          end if
          if (do_deflate) &
               call check( NF90_def_var_deflate(ncid, id_field, 0, 1, 1) )
          call check( nf90_put_att(ncid, id_field, "coordinates", "nav_lat nav_lon") )
          fillval = 1.0e20
          call check( nf90_put_att(ncid, id_field, "_FillValue", fillval) )
          call check( nf90_put_att(ncid, id_field, "missing_value", fillval) )
       end do
       
       ! End define mode
       call check( NF90_ENDDEF(ncid) )

       ! write coordinates
       startz(1)=1
       countz(1)=nlvls

       startC(1) = nimpp
       countC(1) = nlei
       startC(2) = njmpp
       countC(2) = nlej

       call check( nf90_put_var(ncid,id_lon,lons,startC,countC))
       call check( nf90_put_var(ncid,id_lat,lats,startC,countC))

       if (mype==0) then
          call check( nf90_put_var(ncid,id_lev,depths,startz,countz))
       end if

    else
       if (verbose>0 .and. mype==0) &
            write (*,'(a,1x,a,a)') 'NEMO-PDAF', 'Open file: ', trim(filename)

       if (npes==1) then
          call check( nf90_open(trim(filename), NF90_WRITE, ncid) )
       else
          call check( nf90_open_par(trim(filename), NF90_WRITE, comm_filter, MPI_INFO_NULL, ncid) )
       end if

    end if


    ! *** Write fields

    call check( nf90_inq_varid(ncid, 'time', id_time) )
    call check( nf90_VAR_PAR_ACCESS(NCID, id_time, NF90_COLLECTIVE) )
!    call check( nf90_put_vara(ncid, id_time, timeField, start=(/step/), count=(/1/)))
    startt(1) = step
    countt(1) = 1 
    call check( nf90_put_var(ncid, id_time, timeField, startt(1:1), countt(1:1)))

    ! Backwards transformation of state
    call transform_field_mv(2, state)

    do i = 1, n_fields

       ! Convert state vector to field
       tmp_4d = 1.0e20
       call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims)

       if (verbose>1) &
            write (*,'(5x,a,a)') '--- write variable: ', trim(sfields(i)%variable)
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), id_field) )
!       call check( nf90_VAR_PAR_ACCESS(NCID, id_field, NF90_COLLECTIVE) )

       ! Attention with coordinates, in Nemo Restart it is var(time,depth,y,x)
       startt(1) = istart
       countt(1) = ni_p
       startt(2) = jstart
       countt(2) = nj_p
       startt(3) = 1
       countt(3) = nlvls
       startt(4) = step
       countt(4) = 1 
    
       if (sfields(i)%ndims==3) then
          startt(3) = 1
          countt(3) = nlvls
          call check( nf90_put_var(ncid, id_field, tmp_4d, startt, countt))
       else
          startt(3) = step
          countt(3) = 1 

          call check( nf90_put_var(ncid, id_field, tmp_4d, startt(1:3), countt(1:3)))
       end if

    end do

    ! *** close file with state sequence ***
    call check( NF90_CLOSE(ncid) )

  end subroutine write_field_mv


!================================================================================

!> Write increment file
!!
  subroutine write_increment(state, dim, filename, id_field)
    
    use netcdf
   
    implicit none

! *** Arguments ***
    real(8),       intent(inout) :: state(:)    !< State vector
    integer(4),       intent(in) :: dim         !< State dimension
    character(len=*), intent(in) :: filename    !< File name
    integer(4),       intent(in) :: id_field    !< Id of field in stte vector

! *** Local variables ***
    integer(4) :: cnt, i, j, k
    integer(4) :: ncid
    integer(4) :: dimids_field(4)
    integer(4) :: dimid_time, dimid_lvls, dimid_lat, dimid_lon
    integer(4) :: id_dateb, id_datef
    integer(4) :: id_lat, id_lon, id_lev, id_time, id_incr
    integer(4) :: dimids(4)
    integer(4) :: startC(2), countC(2)
    integer(4) :: startt(4), countt(4)
    integer(4) :: startz(1), countz(1)
    real(8)    :: fillval
    real(8)    :: timeField(1)
    real(8)    :: bgnTimeInterv(1), finTimeInterv(1), timeInIncr(1)

    ! NOTE: This routine currently only writes a single field from the state vector
    !       which is specified by id_field. this must be a 3D field

    !The time values are read in via namelist at the moment
    ! Fixme
    timeInIncr(1)=incrTime !time for direct initialisation in Nemo (time of restart file which is used for adding to increment file)
    bgnTimeInterv(1)=startEnsTime !Start date of interval on which increment is valid (later for time ramp initialisation of  increment)
    finTimeInterv(1)=endEnsTime !End date of interval on which increment is valid (later for time rap init of increment)! 

    if (verbose>0 .and. mype==0) &
         write (*,'(8x,a)') '--- Write increment file'

    if (npes==1) then
       call check( NF90_CREATE(trim(filename), NF90_NETCDF4, ncid))
    else
       call check( NF90_CREATE_PAR(trim(filename), NF90_NETCDF4, comm_filter, MPI_INFO_NULL, ncid))
    end if

    call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'title', &
         'Increment matrix for NEMO-PDAF coupling'))
     
    ! define dimensions for NEMO-input file
    call check( NF90_DEF_DIM(ncid, 't', 1, dimid_time))
    call check( NF90_DEF_DIM(ncid, 'z', nlvls, dimid_lvls))
    call check( NF90_DEF_DIM(ncid, 'y', nlats, dimid_lat) )
    call check( NF90_DEF_DIM(ncid, 'x', nlons, dimid_lon) )
   
    dimids_field(4)=dimid_time
    dimids_field(3)=dimid_lvls
    dimids_field(2)=dimid_lat
    dimids_field(1)=dimid_lon
       
    ! define variables
    call check( NF90_DEF_VAR(ncid, 'nav_lat', NF90_FLOAT, dimids_field(1:2), id_lat))
    call check( NF90_DEF_VAR(ncid, 'nav_lon', NF90_FLOAT, dimids_field(1:2), id_lon))
    call check( NF90_DEF_VAR(ncid, 'nav_lev', NF90_FLOAT, dimids_field(3), id_lev))
    if (do_deflate) then
       call check( NF90_def_var_deflate(ncid, id_lat, 0, 1, 1) )
       call check( NF90_def_var_deflate(ncid, id_lon, 0, 1, 1) )
       call check( NF90_def_var_deflate(ncid, id_lev, 0, 1, 1) )
    end if

    call check( NF90_DEF_VAR(ncid, 'time', NF90_DOUBLE, id_time))
    call check( NF90_DEF_VAR(ncid, 'z_inc_dateb', NF90_DOUBLE, id_dateb))
    call check( NF90_DEF_VAR(ncid, 'z_inc_datef', NF90_DOUBLE, id_datef)) 

    do i = 1, n_fields

       if (sfields(i)%ndims==3) then
          dimids_field(3)=dimid_lvls
          call check( NF90_DEF_VAR(ncid, trim(sfields(i)%name_incr), NF90_DOUBLE, dimids_field(1:4), id_incr) )
       else
          dimids_field(3)=dimid_time
          call check( NF90_DEF_VAR(ncid, trim(sfields(i)%name_incr), NF90_DOUBLE, dimids_field(1:3), id_incr) )
       end if

       if (do_deflate) &
            call check( NF90_def_var_deflate(ncid, id_incr, 0, 1, 1) )

       fillval = 0.0
       call check( nf90_put_att(ncid, id_incr, "long_name", trim(sfields(id_field)%name_incr)//trim('Increment')) )
       call check( nf90_put_att(ncid, id_incr, "units", trim(sfields(id_field)%unit)) )
       call check( nf90_put_att(ncid, id_incr, "coordinates", "nav_lat nav_lon") )
       call check( nf90_put_att(ncid, id_incr, "_FillValue", fillval) )
       call check( nf90_put_att(ncid, id_incr, "missing_value", fillval) )
    end do
       
    ! End define mode
    call check( NF90_ENDDEF(ncid) )


    ! *** write coordinates ***
    startz(1)=1
    countz(1)=nlvls

    startC(1)=1
    startC(2)=1
    countC(1)=nlons
    countC(2)=nlats

    if (mype==0) then
       call check( nf90_put_var(ncid, id_lev, depths, startz, countz))
       call check( nf90_put_var(ncid, id_lon, lons, startC, countC))
       call check( nf90_put_var(ncid, id_lat, lats, startC, countC))

       call check( nf90_put_var(ncid, id_time, timeInIncr, start=(/1/), count=(/1/)))
       call check( nf90_put_var(ncid, id_dateb, bgnTimeInterv, start=(/1/), count=(/1/)))
       call check( nf90_put_var(ncid, id_datef, finTimeInterv, start=(/1/), count=(/1/)))
    end if

    ! *** Write fields ***

    ! backwards transformation of increment
    call transform_field_mv(2, state)

    do i = 1, n_fields

       ! Convert state vector to field
       tmp_4d = 0.0
       call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims)

       if (verbose>1) &
            write (*,'(5x,a,a)') '--- write variable: ', trim(sfields(i)%name_incr)
       call check( nf90_inq_varid(ncid, trim(sfields(i)%name_incr), id_incr) )

       ! Attention with coordinates, in Nemo Restart it is var(time,depth,y,x)
       startt(1) = istart
       countt(1) = ni_p
       startt(2) = jstart
       countt(2) = nj_p
       startt(3) = 1
       countt(3) = nlvls
       startt(4) = 1
       countt(4) = 1 
           
       if (sfields(i)%ndims==3) then
          countt(3) = nlvls
          call check( nf90_put_var(ncid, id_incr, tmp_4d, startt, countt))
       else
          countt(3) = 1 
          call check( nf90_put_var(ncid, id_incr, tmp_4d, startt(1:3), countt(1:3)))
       end if

    end do

    ! *** close file with state sequence ***
    call check( NF90_CLOSE(ncid) )

  end subroutine write_increment


!================================================================================

!> Overwrite the NEMO restart file
!!
  subroutine update_restart_mv(path, file, state, state_tmp)

    use netcdf
   
    implicit none

! *** Arguments ***
    character(len=*), intent(in) :: path         !< Path for restart file
    character(len=*), intent(in) :: file         !< Restart file name
    real(8),       intent(inout) :: state(:)     !< State vector
    real(8),       intent(inout) :: state_tmp(:) !< tmp state vector (used for storage)

! *** Local variables ***
    integer :: i                      ! Counter
    integer :: ncid                   ! NC file ID
    integer :: lid, uid               ! index range in state vector
    integer :: id_var_n, id_var_b     ! NC variable IDs
    integer  :: startt(4), countt(4)  ! arrays for file writing
    character(len=30) :: rst_file     ! Name of restart file


    ! Attention in run script copy restart file from time of DA to file 'restart_trc_in_befDA.nc' 

    !Write oxy to TRNOXY of restart file (now, time t) -> 
    !restart Nemo with nn_euler=0 (TRBOXY is oxy for t-Delta t)

    ! Store name of restart file
    rst_file = sfields(1)%rst_file

    if (verbose>0) &
         write (*,'(a,3x,a,1x,a)') 'NEMO-PDAF', '--- Overwrite restart file:',trim(path_restart)//trim(rst_file) 

    ! Open file and retrieve field ids
    if (npes==1) then
       call check( nf90_open(trim(path_restart)//trim(rst_file),NF90_WRITE, ncid))
    else
       call check( nf90_open_par(trim(path_restart)//trim(rst_file),NF90_WRITE,comm_filter, MPI_INFO_NULL, ncid))
    end if

    ! field transformation (if save_Incr=.true. this was already done)
    if (.not.save_Incr) call transform_field_mv(2, state)

    do i = 1, n_fields

       if (trim(sfields(i)%rst_file) /= trim(rst_file)) then
       ! Open other restart file and retrieve field ids
          if (verbose>0 .and. mype==0) &
               write (*,'(a, 3x,a,1x,a)') 'NEMO-PDAF', '--- Open restart file:',trim(path_restart)//trim(sfields(i)%rst_file) 
          if (npes==1) then
             call check( nf90_open(trim(path_restart)//trim(sfields(i)%rst_file),NF90_WRITE, ncid))
          else
             call check( nf90_open_par(trim(path_restart)//trim(sfields(i)%rst_file), &
                  NF90_WRITE, comm_filter, MPI_INFO_NULL, ncid))
          end if
          
          ! Store name of restart file
          rst_file = sfields(i)%rst_file

       end if

       ! Retrieve field IDs
       call check( nf90_inq_varid(ncid, trim(sfields(i)%name_rest_n), id_var_n))
       call check( nf90_inq_varid(ncid, trim(sfields(i)%name_rest_b), id_var_b)) 


       ! backwards transformation state - only if not done by write_increment before

       ! Convert state vector to field
       tmp_4d = 0.0
       call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims)

       ! *** write variable for current time ***
       startt(1) = istart
       countt(1) = ni_p
       startt(2) = jstart
       countt(2) = nj_p
       startt(3) = 1
       countt(3) = nlvls
       startt(4) = 1
       countt(4) = 1 
       
       if (sfields(i)%ndims==3) then
          call check( nf90_put_var(ncid, id_var_n, tmp_4d, startt, countt))
       else
          countt(3) = 1
          call check( nf90_put_var(ncid, id_var_n, tmp_4d, startt(1:3), countt(1:3)))
       end if

       ! *** For second (past) time use increment ***

       ! Read field, add increment, and write field
       if (sfields(i)%ndims==3) then
          countt(3) = nlvls
          call check( nf90_get_var(ncid, id_var_b, tmp_4d, startt, countt))
       else
          countt(3) = 1
          call check( nf90_get_var(ncid, id_var_b, tmp_4d, startt(1:3), countt(1:3)))
       end if

       call field2state(tmp_4d, state, sfields(i)%off, sfields(i)%ndims, missing_value)

       lid = sfields(i)%off+1
       uid = sfields(i)%off+sfields(i)%dim
       state(lid : uid) = state(lid : uid) + state_tmp(lid : uid)

       call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims)

       if (sfields(i)%ndims==3) then
          countt(3) = nlvls
          call check( nf90_put_var(ncid, id_var_b, tmp_4d, startt, countt))
       else
          countt(3) = 1
          call check( nf90_put_var(ncid, id_var_b, tmp_4d, startt(1:3), countt(1:3)))
       end if

    end do

    call check( nf90_close(ncid))

  end subroutine update_restart_mv


!================================================================================

!> Check status of NC operation
!!   
  subroutine check(status)

    use netcdf

! *** Aruments ***
    integer, intent ( in) :: status   ! Reading status

    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if

  end subroutine check
 

! ==============================================================================

!> Add a trailing slash to a path string
!!
!! This routine ensures that a string defining a path
!! has a trailing slash.
!!
  subroutine add_slash(path)

    implicit none

! *** Arguments ***
    character(len=100) :: path  !< String holding the path

! *** Local variables ***
    integer :: strlength

! *** Add trailing slash ***
    strlength = len_trim(path)

    if (path(strlength:strlength) /= '/') then
       path = trim(path) // '/'
    end if
    
  end subroutine add_slash


!===============================================================================

!> Convert an integer to a strong of length 4
!!
  character(len=4) function str(k)

    implicit none

    integer, intent(in) :: k   !< number

    write (str, '(i4.4)') k

  end function str


!===============================================================================

!> Check whether a file exists
!!
  function file_exists(filename) result(res)

    implicit none

    character(len=*),intent(in) :: filename   !< File name
    logical                     :: res        !< Status of file

    ! Check if the file exists
    inquire( file=trim(filename), exist=res )

  end function file_exists

end module mod_io_pdaf
