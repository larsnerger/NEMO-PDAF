!> Compute RMS errors for SST
!! Programm to compute RMS errors for SST
!! for CMEMS SST data vs. the model.
!! This variant is adapted for NEMO DA output       
!! 
program rmse

  use netcdf
  use timer

  implicit none

  integer :: i, j, i_index, j_index, counter, iens
  integer :: mcnt, cnt_days, cnt, cnt0, cntmin
  integer :: dim_ens
  integer :: ilon_m, ilat_m
  integer :: ilon_o, ilat_o
  integer :: ndays_m(12)
  integer :: day, month, year
  integer :: oncid, ncid, dimid, lonid, latid, varid_o, varid_m ! nc file IDs
  integer :: fileid, maskncid, wfileid
  integer :: dimid_day, dimid_mon, dimid_2, dimids(4), wdimids(3)
  integer :: id_mon, id_days, id_rmse, id_mean
  integer :: id_mobs_f, id_mobs_a, lonid_mobs, latid_mobs
  integer :: dim_olat, dim_olon                 ! Observation grid dimensions read from file
  integer :: dim_mlat, dim_mlon                 ! Model grid dimensions read from file
  integer :: startv(4), countv(4)               ! Index arrays for reading from nc file
  integer :: startv3(3), countv3(3)               ! Index arrays for reading from nc file
  integer :: nmonths
  integer :: months(12), years(12)
  integer :: idays(12)
  integer :: jpiglo, jpjglo


  real(8), allocatable :: sst_m(:,:)
  real(8), allocatable :: sst_m_obsgrid(:,:)
  real(8), allocatable :: sst_o(:,:)
  real(8), allocatable :: tmask(:,:,:,:)
  real(8), allocatable :: lon_o(:), lat_o(:)  ! Obs. coordinates read from file
  real(8), allocatable :: lon_m(:), lat_m(:)  ! Obs. coordinates read from file
  real(8), allocatable :: gphit(:,:), glamt(:,:)  ! Obs. coordinates read from file
  real(8) :: dlat, dlon, nlat, slat, wlon, elon
  real(8) :: diff, diff_squared, ssum, diff_ssum
  real(8) :: rmse_val, mean_of_difference_squared, mean_of_difference
  real(8), allocatable :: rmses(:,:), means(:,:)
  real(8) :: missing_value                           ! missing value
  real(4) :: fillval                           ! missing value
  real(8) :: missing_value_obs                       ! missing value in obs file
  real(8) :: minchl      ! Minimum CHL limit for model state
  character(len=200) :: path_o, path_m, path_mobs
  character(len=200) :: file_o, file_m, file_rms, file_mobs, file_model
  character(len=200) :: path_mask, file_mask
  character(len=50) :: exp, mfile, ofile
  character(len=2) :: mstr, dstr
  character(len=4) :: ystr
  character(len=3) :: ensstr

  logical :: first = .true.

  logical :: async = .true.

! *** set number of timers ***
  call timeit(7,'ini')

! *** set first timer ***
  call timeit(1,'new')

  write (*,'(10x,a)') '************************************'
  write (*,'(10x,a)') '*     RMSE for CMEMS CHL Baltic    *'
  write (*,'(10x,a)') '************************************'

  minchl = 0.00001

  ! Number of days per month
  ndays_m = (/30, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

  !months = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /)
  months = (/4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3 /)
  years = (/2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2016, 2016 /)

  ! How many months to compute
  nmonths = 2

  ! Ensemble size
  dim_ens = 30

  ! TRUE for directory structure at HLRN / FALSE for Ollie
  async = .true.

  ! Path to observation files
  path_o = '/scratch/usr/hzfblner/SEAMLESS/observations/CHL_2015'

  ! Path to data assimilation output files
  path_m ='/scratch/usr/hzfblner/SEAMLESS/run/ens_NORDIC_ndays/output_N20_DA_Smago/DA'
  path_m ='/scratch/usr/hzfblner/SEAMLESS/run/ens_ERGOM_DA/001/output/data'
  file_m = 'NORDIC_1d_ERGOM_T'


  ! Path and file name of NEMO T-grid output file to read mask
  path_mask ='/scratch/usr/hzfblner/SEAMLESS/run/ens_ERGOM_DA/001/output/data'
  !path_mask = '/scratch/usr/hzfblner/SEAMLESS/outputs/out_ERGOM_free'
  file_mask = '001_NORDIC_1d_grid_T_20150401-20150401.nc'

  ! Experiment name (used for output file names)
  exp = 'hnk_ska'

  ! Name stub of file holding observed model state in observation grid
  path_mobs = '.' !exp
  file_mobs = 'chl_mobs'




  ! *** End of configuration part ***
  
  ! Count total number of days
  cnt_days = 0
  do mcnt = 1, nmonths

     month = months(mcnt)
     cnt_days = cnt_days + ndays_m(month)
  end do

  write (*,*) 'total number of days', cnt_days
  allocate(rmses(cnt_days, 2), means(cnt_days, 2))

  cnt = 1


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
write (*,*) 'missing_Value', missing_value
  call check( nf90_close(ncid) )

  
  ensloop: do iens = 1, dim_ens

     write (ensstr, '(i3.3)') iens

     write (*,*) 'Process ensemble member',iens

     mloop: do mcnt = 1, nmonths

        month = months(mcnt)
        year = years(mcnt)

        write (ystr, '(i4.4)') year
        write (mstr, '(i2.2)') month

        write (*,*) 'Compute RMSes for month, year', month, year

        ! Open monthly observation file
        ofile = 'chl_ba_MY_'//ystr//mstr//'.nc'
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
           allocate(sst_o(dim_olon, dim_olat))
           allocate(sst_m_obsgrid(dim_olon, dim_olat))
           allocate(lon_o(dim_olon), lat_o(dim_olat))

           ! Read coordinates
           call check( nf90_inq_varid(oncid, "lon", lonid) )
           call check( nf90_inq_varid(oncid, "lat", latid) )

           call check( nf90_get_var(oncid, lonid, lon_o) )
           call check( nf90_get_var(oncid, latid, lat_o) )

           write (*,'(a,2i8)') 'Size of observation grid (lat, lon)', dim_olat, dim_olon
           write (*,'(a,4f10.2)') 'Limits of observation grid (lat, lon)', lat_o(1), lat_o(dim_olat), lon_o(1), lon_o(dim_olon)

        endif

        ! Get variable IDs
        call check( nf90_inq_varid(oncid, "CHL", varid_o) )
        call check( nf90_get_att(oncid, varid_o, 'missing_value', missing_value_obs) )


        ! *** Create file to store observed model state ***

        call check( NF90_CREATE(trim(path_mobs)//'/'//ensstr//'_'//trim(file_mobs)//'_'//trim(exp)//'_'//ystr//mstr//'.nc', NF90_NETCDF4, wfileid) )
        call check( NF90_PUT_ATT(wfileid,NF90_GLOBAL,'title', &
             'Observed model state for chlorophyll') )

        ! Define dimensions
        call check( NF90_DEF_DIM(wfileid,'lon',dim_olon,wdimids(1)) )
        call check( NF90_DEF_DIM(wfileid,'lat',dim_olat,wdimids(2)) )
        call check( NF90_DEF_DIM(wfileid,'time',ndays_m(month),wdimids(3)) )

        ! define variables
        fillval = -999.0
        call check( NF90_DEF_VAR(wfileid,'lat',NF90_FLOAT,wdimids(2),latid_mobs) )
        call check( NF90_def_var_deflate(wfileid, latid_mobs, 0, 1, 1) )
        call check( NF90_DEF_VAR(wfileid,'lon',NF90_FLOAT,wdimids(1),lonid_mobs) )
        call check( NF90_def_var_deflate(wfileid, lonid_mobs, 0, 1, 1) )
        call check( NF90_DEF_VAR(wfileid,'CHL',NF90_FLOAT,wdimids(1:3),id_mobs_f) )
        call check( NF90_def_var_deflate(wfileid, id_mobs_f, 0, 1, 1) )
        call check( nf90_put_att(wfileid, id_mobs_f, "_FillValue", fillval) )
        call check( nf90_put_att(wfileid, id_mobs_f, "missing_value", fillval) )

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

           ! Read SST values
           ! They are in deg C but have to be scaled by sst_scale (1/100).
           startv3(1) = 1 ! lon
           startv3(2) = 1 ! lat
           startv3(3) = day ! time
           countv3(1) = dim_olon
           countv3(2) = dim_olat
           countv3(3) = 1
           call check( nf90_get_var(oncid, varid_o, sst_o, start=startv3, count=countv3) )

           ! *** Read model ***


           file_model = trim(path_m)//'/'//ensstr//'_'//trim(file_m)//'_'//ystr//mstr//dstr//'-'//ystr//mstr//dstr//'.nc'

           if (day==1) then
              write (*, '(8x,a,i4,1x,i2,1x,i2,a)') 'Read model CHL from file:'
              write (*, '(10x,a)') trim(file_model)
           end if

           call check( nf90_open(file_model, NF90_NOWRITE, ncid) )

           if (first) then
              ! Read dimensions of observation grid
              call check( nf90_inq_dimid(ncid, "y", dimid) )
              call check( nf90_Inquire_dimension(ncid, dimid, len=dim_mlat) )
              call check( nf90_inq_dimid(ncid, "x", dimid) )
              call check( nf90_Inquire_dimension(ncid, dimid, len=dim_mlon) )

              ! Allocate arrays
              allocate(sst_m(dim_mlon, dim_mlat))
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
                    if (gphit(i,j) >= -180.0 .and. gphit(i,j) <= 180.0 &
                         .and. gphit(i,j) /= 0.0) then
                       lat_m(j) = gphit(i,j)
                    endif
                 enddo
              enddo

              lon_m(:) = 0.0
              do j = 1, dim_mlat
                 do i = 1, dim_mlon
                    if (glamt(i,j) >= -180.0 .and. glamt(i,j) <= 180.0 &
                         .and. glamt(i,j) /= 0.0) then
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
           call check( nf90_inq_varid(ncid, "CHL", varid_m) )

           sst_m_obsgrid = fillval


           ! Read SST values
           startv(1) = 1 ! lon
           startv(2) = 1 ! lat
           startv(3) = 1 ! layer
           startv(4) = 1 ! time
           countv(1) = dim_mlon
           countv(2) = dim_mlat
           countv(3) = 1
           countv(4) = 1
           call check( nf90_get_var(ncid, varid_m, sst_m, start=startv, count=countv) )

           cntmin = 0
           do j = 1, dim_mlat
             do i = 1, dim_mlon
               if (abs(sst_m(i,j)-missing_value)>0.1 .and. sst_m(i,j)<minchl) then
                  sst_m(i,j)=minchl
                  cntmin=cntmin+1
               endif
             end do
           end do
!           write (*,*) 'number of corrected CHL values', cntmin

           ! Compute RMSE

           ! Initialize counters and sums
           ssum = 0.0
           diff_ssum = 0.0
           counter = 0
           cnt0 = 0
           sst_m_obsgrid(:,:) = fillval           

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
                          .and. sst_o(i,j) > missing_value_obs) then
                          sst_m_obsgrid(i,j) = sst_m(i_index, j_index)
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
           call check( nf90_put_var(wfileid, id_mobs_f, sst_m_obsgrid(:,:), startv3, countv3))
        
           first = .false.

           write (*,'(a, 1x, 2i3)') 'month, day', month, day

           cnt = cnt + 1

           ! Close the file
           call check( nf90_close(ncid) )

        end do dloop

        ! Close observation file
        call check( nf90_close(oncid) )

        ! close output file
        call check( NF90_CLOSE(wfileid) )

     end do mloop

  end do ensloop


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
end program rmse
