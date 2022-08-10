!> Module holding IO operations for NEMO-PDAF
!!
module mod_io_pdaf
use mod_kind_pdaf
use netcdf
implicit none

integer :: verbose=1   ! Set verbosity of IO routines (0,1,2,3)


contains
   !> get number of units for varname
   !!   
   function get_units(filename, varname) result(units)
      character(*), intent(in) :: filename
      character(*), intent(in) :: varname

      character(256) :: units

      integer :: ncid
      integer :: varid

      ! Open the file
      call check( nf90_open(filename, nf90_nowrite, ncid) )
      call check( nf90_inq_varid(ncid, varname, varid) )

      if (nf90_inquire_attribute(ncid, varid, 'units') == nf90_noerr) &
         call check( nf90_get_att(ncid, varid, 'units', units) )
      call check( nf90_close(ncid) )
   end function get_units
   ! ===================================================================================


   !> get number of dimensions for varname
   !!   
   function get_ndims(filename, varname) result(ndims)
      character(*), intent(in) :: filename
      character(*), intent(in) :: varname

      integer :: ndims

      integer :: ncid
      integer :: varid

      call check( nf90_open(filename, nf90_nowrite, ncid) )
      call check( nf90_inq_varid(ncid, varname, varid) )
      call check( nf90_inquire_variable(ncid, varid, ndims=ndims) )
      call check( nf90_close(ncid) )
   end function get_ndims
   ! ===================================================================================


   !> Read total time steps of NEMO files in a directory
   !> Here we assume all variables have the same length
   !!
   function get_steps(filename) result(steps)
      character(*), intent(in) :: filename
      integer :: steps

      integer :: ncid   ! NC file id
      integer :: dimid  ! dimension id

      call check( nf90_open(filename, nf90_nowrite, ncid) )
      ! Get the number of time steps
      call check( nf90_inq_dimid(ncid, 'time_counter', dimid) )
      call check( nf90_inquire_dimension(ncid, dimid, len=steps) )
      call check( nf90_close(ncid) )
   end function get_steps
   ! ===================================================================================
 

   !> Read fields from NEMO file into a state vector
   !!
   subroutine read_state_mv(filename, varname, ndims, field, nx, ny, nlvls, nt)
      character(len = *), intent(in) :: filename  !< filename
      character(len = *), intent(in) :: varname   !< variable name
      integer(4), intent(in) :: ndims             !< number of dimensions
      integer(4), intent(in) :: nx, ny, nlvls     !< size of spatial dimension
      integer(4), optional, intent(in) :: nt      !< number of time steps for read
      real(pdp), intent(out)   :: field(:, :, :, :) !< field

      ! ! Local variables 
      integer(4) :: it              ! Counters
      integer(4) :: nt_             ! number of time steps
      
      integer(4) :: varid           ! Variable ID
      integer(4) :: ncid            ! NC file id
      real(pdp) :: missing_value      ! missing value
      

      if (verbose>0) &
         write(*,'(a,4x,a,a)') 'NEMO-PDAF', '*** Read model output from: ', filename

      if (verbose>1) then 
         write (*,'(a,a,a)') &
         'NEMO-PDAF', 'Variable: ', varname
      end if

      ! Open the file
      call check( nf90_open(filename, nf90_nowrite, ncid) )

      call check( nf90_inq_varid(ncid, varname, varid) )

      missing_value=0.0_pdp
      if (nf90_inquire_attribute(ncid, varid, 'missing_value') == nf90_noerr) &
         call check( nf90_get_att(ncid, varid, 'missing_value', missing_value) )

      if (present(nt)) then
         nt_ = nt
      else
         nt_ = get_steps(filename)
      end if

      do it = 1, nt_
         ! Read variable
         if (ndims == 3) then
            call check( nf90_get_var(ncid, varid, field(:, :, :, it), &
            start=[1, 1, 1, it], count=[nx, ny, nlvls, 1]) )
         else
            call check( nf90_get_var(ncid, varid, field(:, :, 1, it), &
            start=[1, 1, it], count=[nx, ny, 1]) )
         end if
      end do

      call check( nf90_close(ncid) )

      ! if (verbose>2) then
      !    do i = 1, n_fields
      !       write(*,'(a, 1x, a, a10, 1x, a,5x, 2f12.6)') &
      !       'NEMO-PDAF', 'Min and max for ',trim(sfields(i)%variable),' :     ', &
      !       minval(state(off(i)+1:off(i)+sfields(i)%dim)), &
      !       maxval(state(off(i)+1:off(i)+sfields(i)%dim))
      !    enddo
      ! end if
   end subroutine read_state_mv
   !===============================================================================


   !> write singular value and vectors
   !!   
   subroutine write_state(filename, varnames, ndims, off, &
                          n2d, n3d, nx, ny, nlvls, dim_state, &
                          maxtimes, nfields, svals, &
                          svdu, rank, run_meanstate, do_mv)
      character(*), intent(in) :: filename
      character(*), dimension(:), intent(in) ::varnames
      integer, intent(in) :: ndims(:)
      integer, intent(in) :: off(:)
      real(pdp), intent(in) :: svals(:)
      real(pdp), intent(in) :: run_meanstate(:, :)
      real(pdp), intent(in) :: svdu(:, :)
      integer, intent(in) :: do_mv
      integer, intent(in) :: rank
      integer, intent(in) :: n2d, n3d
      integer, intent(in) :: nx, ny, nlvls
      integer, intent(in) :: dim_state
      integer, intent(in) :: maxtimes
      integer, intent(in) :: nfields

      integer :: i, j ! counter
      integer :: ncid ! netcdf file id
      integer :: dimid_one ! size one dimension
      integer :: dimid_rank ! dimension id of rank of svd
      integer :: dimid_x ! dimension id of x dimension
      integer :: dimid_y ! dimension id of y dimension
      integer :: dimid_z ! dimension id of vertical dimension
      integer :: dimid_state
      integer :: dimid_nfields
      integer :: varid_sigma
      integer :: varid_mean(nfields)
      integer :: varid_svd(nfields)
      real(pdp), allocatable :: field(:, :, :)
      character(150) :: attstr

      ! *********************************************************
      ! *** write mean state and decomposed covariance matrix ***
      ! ********************************************************

      write (*,'(/1x,a)') '------- write decomposed covariance matrix -------------'

      allocate(field(nx, ny, nlvls))

      ! *** initialize file
      call check( nf90_create(filename, nf90_netcdf4, ncid) )

      if (do_mv == 1) then
         attstr = 'running mean state, scaled singular vectors and values of multivariate decomposed covariance matrix for nemo'
      else
         attstr = 'running mean state, singular vectors and values of decomposed covariance matrix for nemo'
      end if

      ! write dim: rank. nodes_2d, nodes_3d, dim_state, one, nfields
      call check( nf90_put_att(ncid, nf90_global, 'title', trim(attstr)) )

      do i = 1, nfields
         attstr = trim(attstr)//trim(varnames(i))
      end do
      call check( nf90_put_att(ncid, nf90_global, 'state_fields', trim(attstr)) )

      ! define dimensions
      call check( nf90_def_dim(ncid, 'rank',  rank, dimid_rank) )
      call check( nf90_def_dim(ncid, 'x', nx, dimid_x) )
      call check( nf90_def_dim(ncid, 'y', ny, dimid_y) )
      call check( nf90_def_dim(ncid, 'z', nlvls, dimid_z) )
      call check( nf90_def_dim(ncid, 'dim_state', dim_state, dimid_state) )
      call check( nf90_def_dim(ncid, 'one',  1, dimid_one) )
      call check( nf90_def_dim(ncid, 'nfields',  nfields, dimid_nfields) )

      ! define variables
      ! sigular values
      call check( nf90_def_var(ncid, 'sigma', nf90_double, dimid_rank, varid_sigma) )

      ! mean state and singular vectors
      do i = 1, nfields
         if (ndims(i) == 3) then
            call check( nf90_def_var(ncid, trim(varnames(i))//'_mean', nf90_double, &
                                     [dimid_x, dimid_y, dimid_z, dimid_one], varid_mean(i)) )
            call check( nf90_def_var(ncid, trim(varnames(i))//'_svd', nf90_double, &
                                     [dimid_x, dimid_y, dimid_z, dimid_rank], varid_svd(i)) )
         else
            call check( nf90_def_var(ncid, trim(varnames(i))//'_mean', nf90_double, &
                                     [dimid_x, dimid_y, dimid_one], varid_mean(i)) )
            call check( nf90_def_var(ncid, trim(varnames(i))//'_svd', nf90_double, &
                                     [dimid_x, dimid_y, dimid_rank], varid_svd(i)) )
         endif
      end do
      ! end define mode
      call check( nf90_enddef(ncid) )

      ! write singular values
      call check( nf90_put_var(ncid, varid_sigma, svals(1:rank)) )


      ! write running mean (for the last snap shot)
      do i = 1, nfields
         if (ndims(i) == 3) then
            field(:, :, :) = reshape(run_meanstate(1+off(i) : n3d+off(i), maxtimes), [nx, ny, nlvls])
            call check( nf90_put_var(ncid, varid_mean(i), &
                                     field, &
                                     [1, 1, 1], [nx, ny, nlvls, 1]) )
         else
            field(:, :, 1) = reshape(run_meanstate(1+off(i) : n2d+off(i), maxtimes), [nx, ny])
            call check( nf90_put_var(ncid, varid_mean(i), &
                                     field(:, :, 1), &
                                     [1, 1], [nx, ny, 1]) )
         endif
      end do

      ! *** write singular vectors
      writevectors: do i = 1, rank
         do j = 1, nfields
            if (ndims(j) == 3) then
               field(:, :, :) = reshape(svdu(1+off(j) : n3d+off(j), i), [nx, ny, nlvls])
               call check( nf90_put_var(ncid, varid_svd(j), &
                                        field, &
                                        [1, 1, 1, i], [nx, ny, nlvls, 1]) )
            else
               field(:, :, 1) = reshape(svdu(1+off(j) : n3d+off(j), i), [nx, ny])
               call check( nf90_put_var(ncid, varid_svd(j), &
                                        field(:, :, 1), &
                                        [1, 1, i], [nx, ny, 1]) )
            endif
         end do
      end do writevectors

      ! close file
      call check( nf90_close(ncid))

      deallocate(field)
   end subroutine write_state   
   ! ==============================================================================


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


   !> Convert an integer to a string of length 4
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
