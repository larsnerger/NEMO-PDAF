!> Project observed model field onto observation grid
!! This variant is adapted for NEMO DA output       
!! 
program obsgrid

  use netcdf
  use timer

  implicit none

  integer :: i, j, i_index, j_index, counter, step, mcnt, cnt_days, cnt, cnt0, iday
  integer :: ilon_m, ilat_m
  integer :: ilon_o, ilat_o
  integer :: ndays_m(12)
  integer :: day, month, year
  integer :: oncid, ncid, dimid, lonid, latid, varid_o, varid_m ! nc file IDs
  integer :: fileid, maskncid
  integer :: dimid_day, dimid_mon, dimid_2, dimids(4)
  integer :: id_mon, id_days
  integer :: wfileid
  integer :: wdimids(3)
  integer :: id_mobs_f, id_mobs_a, lonid_mobs, latid_mobs
  integer :: dim_olat, dim_olon                 ! Observation grid dimensions read from file
  integer :: dim_mlat, dim_mlon                 ! Model grid dimensions read from file
  integer :: startv(4), countv(4)               ! Index arrays for reading from nc file
  integer :: startv3(3), countv3(3)               ! Index arrays for reading from nc file
  integer :: firstmonth, lastmonth, nmonths
  integer :: months(12), years(12)
  integer :: idays(12)
  integer :: jpiglo, jpjglo
  integer :: maxstep

  real(8), allocatable :: field_m(:,:)
  real(8), allocatable :: field_m_obsgrid(:,:,:)
  real(8), allocatable :: field_o(:,:)
  real(8), allocatable :: tmask(:,:,:,:)
  real(8), allocatable :: lon_o(:), lat_o(:)  ! Obs. coordinates read from file
  real(8), allocatable :: lon_m(:), lat_m(:)  ! Obs. coordinates read from file
  real(8), allocatable :: gphit(:,:), glamt(:,:)  ! Obs. coordinates read from file
  real(8) :: dlat, dlon, nlat, slat, wlon, elon
  real(8) :: missing_value                           ! missing value
  real(4) :: fillval                           ! missing value
  real(8) :: missing_value_obs                       ! missing value in obs file
  character(len=200) :: path_o, file_o
  character(len=200) :: path_m, file_m
  character(len=200) :: path_mask, file_mask
  character(len=200) :: path_mobs, file_mobs
  character(len=50) :: exp, mfile, varname_m, ofile, file_o_stub, varname_o
  character(len=2) :: mstr, dstr
  character(len=4) :: ystr
  character(len=3) :: obstype

  logical :: first = .true.

! Variables particular for SST observations
  integer, allocatable :: ifield_o(:,:)         ! Integer value read from observation file


! *********************
! *** Configuration ***
! *********************

  ! First and last month of experiment
  firstmonth = 1
  lastmonth = 12


  ! *** Model settings

  ! Name of experiment
  exp = 'free_N30'
!  exp = 'sst-chl_Tonly_N30'

  ! Path to data assimilation output files
  path_m = '/scratch/projects/hbk00095/exp/exp.'//trim(exp)//'/DA'
!  path_m = '/scratch/usr/hbknerge/SEAMLESS/run/DA-SST-CHL/exp.'//trim(exp)//'/DA'
  
  ! Name model variable
  varname_m = 'votemper'



  ! *** Observation settings

  ! Choose observation type: L4 or L3S
  obstype = 'L4'

  ! Path to and name stub of observation files, and name of variable
  if (trim(obstype)=='L4') then
     path_o = '/scratch/usr/hzfblner/SEAMLESS/observations/SST_2015'
     file_o_stub = 'sst_REP_L4_'
     varname_o = 'analysed_sst'
  else
     path_o = '/scratch/usr/hzfblner/SEAMLESS/observations/SST_L3S_2015'
     file_o_stub = 'sst_L3S_'
     varname_o = 'sea_surface_temperature'
  end if

  ! Set missing value for observations
  missing_value_obs = -10.0


  ! *** Output settings

  ! Name stub of file holding observed model state in observation grid
  path_mobs = '.'
  file_mobs = 'sst_'//trim(obstype)//'_mobs'


  ! *** Specification of mask file

  ! Path and file name of NEMO T-grid output file to read mask
  path_mask = '/scratch/projects/hbk00095/exp/NEMOout'
  file_mask = '001_NORDIC_1d_grid_T_20150101-20150101.nc'


  ! Number of steps to process - for free run set maxstep=1, otherwise =2
  if (trim(exp)=='free_N30') then
     maxstep = 1
  else
     maxstep = 2
  end if


  ! *** End of configuration part ******************************

  ! set number of timers
  call timeit(7,'ini')

  ! set first timer
  call timeit(1,'new')

  write (*,'(10x,a)') '************************************************************'
  write (*,'(10x,a)') '* Generate surface CHL on observation grid of CMEMS SST    *'
  write (*,'(10x,a)') '************************************************************'

  write (*,'(5x,a,1x,a)') 'write mapped field into file ', trim(file_mobs)//'_'//trim(exp)

  ! Definitions for months
  ndays_m = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
  months = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /)
  years = (/2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015 /)

  ! Initialize number of months
  nmonths = lastmonth - firstmonth + 1

  ! Count total number of days
  cnt_days = 0
  do mcnt = firstmonth, lastmonth
     month = months(mcnt)
     cnt_days = cnt_days + ndays_m(month)
  end do
  write (*,*) 'total number of days', cnt_days

  ! Get first day of year in time series
  iday = 1
  do mcnt = 1, firstmonth-1
     iday = iday + ndays_m(mcnt)
  end do
  write (*,*) 'Start processing at day of year', iday


! *** Read NEMO output to get mask information

  write (*,'(a,a)') 'Read mask information from:', trim(path_mask)//'/'//trim(file_mask)

  call check( nf90_open(trim(path_mask)//'/'//trim(file_mask), NF90_NOWRITE, ncid) )
  call check( nf90_inq_dimid(ncid, "y", dimid) )
  call check( nf90_Inquire_dimension(ncid, dimid, len=dim_mlat) )
  call check( nf90_inq_dimid(ncid, "x", dimid) )
  call check( nf90_Inquire_dimension(ncid, dimid, len=dim_mlon) )

  allocate(tmask(dim_mlon, dim_mlat, 1, 1))

  call check( nf90_inq_varid(ncid, 'votemper', varid_m) )
  call check( nf90_get_att(ncid, varid_m, 'missing_value', missing_value) )
  call check( nf90_get_var(ncid, varid_m, tmask, (/1, 1, 1, 1/), &
       (/dim_mlon, dim_mlat, 1, 1/) ) )

  call check( nf90_close(ncid) )


! *** Start loop computing the timeseries ***

  cnt = 1

  mloop: do mcnt = firstmonth, lastmonth

     month = months(mcnt)
     year = years(mcnt)

     write (ystr, '(i4.4)') year
     write (mstr, '(i2.2)') month

     write (*,*) 'Project onto observation grid for month, year', month, year

     ! Open monthly observation file
     ofile = trim(file_o_stub)//ystr//mstr//'.nc'
     file_o = trim(path_o)//'/'//ofile

     write (*, '(8x,a,i3,a)') 'Read observations from file:'
     write (*, '(10x,a)') trim(file_o)

     call check( nf90_open(file_o, NF90_NOWRITE, oncid) )

     if (first) then
        ! Read dimensions of observation grid
        call check( nf90_inq_dimid(oncid, "lat", dimid) )
        call check( nf90_Inquire_dimension(oncid, dimid, len=dim_olat) )
        call check( nf90_inq_dimid(oncid, "lon", dimid) )
        call check( nf90_Inquire_dimension(oncid, dimid, len=dim_olon) )

        ! Allocate arrays
        allocate(field_o(dim_olon, dim_olat))
        allocate(ifield_o(dim_olon, dim_olat))
        allocate(field_m_obsgrid(dim_olon, dim_olat, 2))
        allocate(lon_o(dim_olon), lat_o(dim_olat))

        ! Read coordinates
        call check( nf90_inq_varid(oncid, "lon", lonid) )
        call check( nf90_inq_varid(oncid, "lat", latid) )

        call check( nf90_get_var(oncid, lonid, lon_o) )
        call check( nf90_get_var(oncid, latid, lat_o) )

        write (*,'(a,2i8)') 'Size of observation grid (lat, lon)', dim_olat, dim_olon
        write (*,'(a,4f10.2)') 'Limits of observation grid (lat, lon)', lat_o(1), lat_o(dim_olat), lon_o(1), lon_o(dim_olon)
     end if

     ! Get variable IDs
     call check( nf90_inq_varid(oncid, trim(varname_o), varid_o) )
     fillval = REAL(missing_value_obs,4)


     ! *** Create file to store observed model state ***

     call check( NF90_CREATE(trim(path_mobs)//'/'//trim(file_mobs)//'_'//trim(exp)//'_'//ystr//mstr//'.nc', NF90_NETCDF4, wfileid) )
     call check( NF90_PUT_ATT(wfileid,NF90_GLOBAL,'title', &
          'Observed model state for chlorophyll') )

     ! Define dimensions
     call check( NF90_DEF_DIM(wfileid,'lon',dim_olon,wdimids(1)) )
     call check( NF90_DEF_DIM(wfileid,'lat',dim_olat,wdimids(2)) )
     call check( NF90_DEF_DIM(wfileid,'time',ndays_m(month),wdimids(3)) )

     ! define variables
     call check( NF90_DEF_VAR(wfileid,'lat',NF90_FLOAT,wdimids(2),latid_mobs) )
     call check( NF90_def_var_deflate(wfileid, latid_mobs, 0, 1, 1) )
     call check( NF90_DEF_VAR(wfileid,'lon',NF90_FLOAT,wdimids(1),lonid_mobs) )
     call check( NF90_def_var_deflate(wfileid, lonid_mobs, 0, 1, 1) )
     call check( NF90_DEF_VAR(wfileid,'SST_f',NF90_FLOAT,wdimids(1:3),id_mobs_f) )
     call check( NF90_def_var_deflate(wfileid, id_mobs_f, 0, 1, 1) )
     call check( nf90_put_att(wfileid, id_mobs_f, "_FillValue", fillval) )
     call check( nf90_put_att(wfileid, id_mobs_f, "missing_value", fillval) )
     call check( NF90_DEF_VAR(wfileid,'SST_a',NF90_FLOAT,wdimids(1:3),id_mobs_a) )
     call check( NF90_def_var_deflate(wfileid, id_mobs_a, 0, 1, 1) )
     call check( nf90_put_att(wfileid, id_mobs_a, "_FillValue", fillval) )
     call check( nf90_put_att(wfileid, id_mobs_a, "missing_value", fillval) )

     ! End define mode
     call check( NF90_ENDDEF(wfileid) )

     ! Write coordinates
     startv(1) = 1
     countv(1) = dim_olon
     call check( nf90_put_var(wfileid, lonid_mobs, lon_o))
     startv(1) = 1
     countv(1) = dim_olat
     call check( nf90_put_var(wfileid, latid_mobs, lat_o))


     ! *** Loop over days of the month

     dloop: do day = 1, ndays_m(month)

        write (dstr, '(i2.2)') day

        ! *** Read observations ***

        ! Read observed field
        ! SSTs are in deg C but have to be scaled by sst_scale (1/100).
        startv3(1) = 1 ! lon
        startv3(2) = 1 ! lat
        startv3(3) = day ! time
        countv3(1) = dim_olon
        countv3(2) = dim_olat
        countv3(3) = 1
        call check( nf90_get_var(oncid, varid_o, ifield_o, start=startv3, count=countv3) )

        ! Scale SST and store as real
        do j = 1, dim_olat
           do i = 1, dim_olon
              field_o(i,j) = real(ifield_o(i, j), 8) * 0.01_8
           end do
        end do

        ! *** Read model ***

        file_m = trim(path_m)//'/'//'state_'//ystr//mstr//dstr//'.nc'

        if (day==1) then
           write (*, '(8x,a,a,a)') 'Read model ', trim(varname_m), ' from file:'
           write (*, '(10x,a)') trim(file_m)
        end if

        call check( nf90_open(file_m, NF90_NOWRITE, ncid) )

        if (first) then
           ! Read dimensions of observation grid
           call check( nf90_inq_dimid(ncid, "y", dimid) )
           call check( nf90_Inquire_dimension(ncid, dimid, len=dim_mlat) )
           call check( nf90_inq_dimid(ncid, "x", dimid) )
           call check( nf90_Inquire_dimension(ncid, dimid, len=dim_mlon) )

           ! Allocate arrays
           allocate(field_m(dim_mlon, dim_mlat))
           allocate(gphit(dim_mlon, dim_mlat), glamt(dim_mlon, dim_mlat))

           call check( nf90_inq_varid(ncid, "nav_lon", lonid) )
           call check( nf90_inq_varid(ncid, "nav_lat", latid) )

           call check( nf90_get_var(ncid, lonid, glamt) )
           call check( nf90_get_var(ncid, latid, gphit) )

           ! Get coordinate vectors
           allocate(lon_m(dim_mlon), lat_m(dim_mlat))

           lat_m(:) = 0.0
           do j = 1, dim_mlat
              do i = 1, dim_mlon
                 if (gphit(i,j) >= -180.0 .and. gphit(i,j) <= 180.0) then
                    lat_m(j) = gphit(i,j)
                 endif
              enddo
           enddo

           lon_m(:) = 0.0
           do j = 1, dim_mlat
              do i = 1, dim_mlon
                 if (glamt(i,j) >= -180.0 .and. glamt(i,j) <= 180.0) then
                    lon_m(i) = glamt(i,j)
                 endif
              enddo
           enddo

           if (first) then
              write (*,'(a,2i8)') 'Size of model grid (lat, lon)', dim_mlat, dim_mlon
              write (*,'(a,4f10.2)') 'Limits of model grid (lat, lon)', lat_m(1), lat_m(dim_mlat), lon_m(1), lon_m(dim_mlon-1)
           end if

           dlat = lat_m(2) - lat_m(1)
           dlon = lon_m(2) - lon_m(1)

           nlat = lat_m(dim_mlat)
           slat = lat_m(1)
           wlon = lon_m(1)
           elon = lon_m(dim_mlon-1)

        end if

        ! Get variable IDs and read data
        call check( nf90_inq_varid(ncid, trim(varname_m), varid_m) )

        ! Loop over two steps in DA output file
        steploop: do step = 1, maxstep

           ! Read model field
           startv(1) = 1 ! lon
           startv(2) = 1 ! lat
           startv(3) = 1 ! layer
           startv(4) = step ! time
           countv(1) = dim_mlon
           countv(2) = dim_mlat
           countv(3) = 1
           countv(4) = 1
           call check( nf90_get_var(ncid, varid_m, field_m, start=startv, count=countv) )

           ! Project model field onto observation grid

           ! Initialize counters and sums
           counter = 0
           cnt0 = 0
           field_m_obsgrid(:,:,step) = missing_value_obs
           if (maxstep==1) field_m_obsgrid(:,:,2) = missing_value_obs

           ! Loop through observation grid

           do j = 1, dim_olat
           
              if (lat_o(j) <= nlat .and. lat_o(j) >= slat) then

                 j_index = floor((lat_o(j) - slat) / dlat) + 1
                 if (j_index > dim_mlat) j_index = dim_mlat

                 do i = 1, dim_olon

                    if (lon_o(i) >= wlon .and. lon_o(i) <= elon) then

                       i_index = floor((lon_o(i) - wlon) / dlon) + 1
                       if (i_index > dim_mlon) i_index = dim_mlon

                       cnt0 = cnt0+1
                       if (abs(tmask(i_index, j_index, 1, 1) - missing_value) > 0.1 &
                          .and. field_o(i,j) > missing_value_obs) then

                          field_m_obsgrid(i,j,step) = field_m(i_index, j_index)

                       end if

                    end if
                 end do ! i
              end if
           end do ! j

           ! Write observed model field on observation grid
           startv3(1) = 1 ! lon
           startv3(2) = 1 ! lat
           startv3(3) = day ! time
           countv3(1) = dim_olon
           countv3(2) = dim_olat
           countv3(3) = 1
           if (step==1) then
              call check( nf90_put_var(wfileid, id_mobs_f, field_m_obsgrid(:,:,1), startv3, countv3))
              if (maxstep==1) then
                 ! For free run write field onlo into step 2
                 call check( nf90_put_var(wfileid, id_mobs_a, field_m_obsgrid(:,:,1), startv3, countv3))
              end if
           else
              call check( nf90_put_var(wfileid, id_mobs_a, field_m_obsgrid(:,:,2), startv3, countv3))
           end if
        
           first = .false.

        end do steploop

        write (*,'(a, 1x, 2i3)') 'month, day', month, day

        cnt = cnt + 1
        iday = iday + 1

        ! Close the file
        call check( nf90_close(ncid) )

     end do dloop

     ! Close observation file
     call check( nf90_close(oncid) )

     ! close output file
     call check( NF90_CLOSE(wfileid) )

  end do mloop


! ********************
! *** Finishing up ***
! ********************

! *** set first timer ***
  call timeit(1,'old')

  ! *** Timings ***
  write (*,'(//23x,a)') 'Timing information'
  write (*,'(14x,a)') '----------------------------------'
  write (*,'(21x,a,F13.3,1x,a)') 'Total run time:',time_tot(1),'s'




contains
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      write(*,*) trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check
end program obsgrid
