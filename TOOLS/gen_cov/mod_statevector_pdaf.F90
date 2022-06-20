!> Building the Statevector
!!
!! This module provides variables & routines for
!! building the state vector.
!!
!! The module contains three routines
!! - **init_id** - initialize the array `id`
!! - **init_sfields** - initialize the array `sfields`
!! - **setup_statevector** - generic routine controlling the initialization
!!
!! The declarations of **id** and **sfields** as well as the
!! routines ~~init_id** and **init_sfields** usually need to be
!! adapted to a particular modeling case.
!!
module mod_statevector_pdaf
use mod_kind_pdaf, only: pdp
implicit none
save


! Declare Fortran type holding the definitions for model fields
type state_field
   character(len=10) :: variable=''    ! Name of field
   character(len=50) :: file=''        ! File name stub to read field from
   logical :: read_from_list = .false. ! if we read actual filename from a list of files

   integer :: transform = 0            ! Type of variable transformation; not implemented
   real(pdp) :: trafo_shift = 0.0      ! Constant to shift value in transformation; not implemented

   character(len=20) :: unit           ! Unit of variable
   integer :: ndims = 0                ! Number of field dimensions (2 or 3)
   integer :: dim = 0                  ! Dimension of the field
   integer :: off  = 0                 ! Offset of field in state vector
end type state_field


!---- The next variables usually do not need editing -----

integer :: screen=1          ! Verbosity flag

! Type variable holding field IDs in state vector


contains
   ! ==============================================================
   subroutine init_sfields(sfields, nfields, n2d, n3d)
      use mod_io_pdaf, only: get_var_dims
      type(state_field), allocatable, intent(out) :: sfields(:)
      integer, intent(out) :: nfields ! number of model fields
      integer, intent(out) :: n2d     ! size of 2d field
      integer, intent(out) :: n3d     ! size of 3d field

      integer :: nx, ny ,nlvls

      namelist /n_statevector/ nfields
      read(unit=20, nml=n_statevector)
      rewind(20)

      if (allocated(sfields)) deallocate(sfields)
      allocate(sfields(nfields))

      ! initialise field
      call init_sfields_namelist(sfields)
      call init_sfields_unit(sfields, nfields)
      call init_sfields_ndim(sfields, nfields)
      ! get the spatial dimension size of the domain
      call get_var_dims(trim(get_filename(sfields(1))), nlvls, ny, nx)
      ! set variable sizes
      n2d = ny*nx
      n3d = n2d*nlvls
      call init_sfields_off(sfields, nfields, n2d, n3d)
   end subroutine init_sfields


   !> This initializes the sfields from namelist
   !!
   subroutine init_sfields_namelist(sfields)
      use mod_kind_pdaf
      implicit none
      ! Type variable holding the defintions of model fields
      type(state_field), intent(inout) :: sfields(:)

      namelist /state_vector/ sfields
      read(unit=20, nml=state_vector)
      rewind(20)
   end subroutine init_sfields_namelist
   ! ==============================


   subroutine init_sfields_unit(sfields,  nfields)
      use mod_io_pdaf, only: get_units
      integer, intent(in) :: nfields
      type(state_field), intent(inout) :: sfields(nfields)

      ! local variable
      integer :: i

      do i = 1, nfields
         sfields(i)%unit = trim(get_units(trim(get_filename(sfields(i))), &
                                          trim(sfields(i)%variable)  &
                                          ) &
                                )
      end do
   end subroutine init_sfields_unit


   subroutine init_sfields_ndim(sfields, nfields)
      use mod_io_pdaf, only: get_ndims
      integer, intent(in) :: nfields
      type(state_field), intent(inout) :: sfields(nfields)


      integer :: i

      do i = 1, nfields
         sfields(i)%ndims = get_ndims(trim(get_filename(sfields(i))), &
                                      trim(sfields(i)%variable) &
                                      ) - 1
      end do
   end subroutine init_sfields_ndim


   subroutine init_sfields_off(sfields, nfields, n2d, n3d)
      integer, intent(in) :: nfields
      type(state_field), intent(inout) :: sfields(nfields)
      integer, intent(in) :: n2d
      integer, intent(in) :: n3d

      integer :: count
      integer :: i

      do i = 1, nfields
         if (sfields(i)%ndims == 2) then
            sfields(i)%dim = n2d
         else if (sfields(i)%ndims == 3) then
            sfields(i)%dim = n3d
         else
            stop 'init_sfields_off: number of dimension in sfields are wrong'
         end if
      end do

      count = 0
      do i = 1, nfields
         sfields(i)%off = count
         count = count + sfields(i)%dim
      end do
   end subroutine init_sfields_off


   ! an iterator of filenames for each variable
   function get_filename(sfield, do_init, do_exit) result(filename)
      type(state_field), intent(in) :: sfield
      logical, intent(in), optional :: do_init
      logical, intent(out), optional :: do_exit

      character(100) :: filename

      integer, save :: counter = 1 ! number function calls
      logical :: do_init_ = .true.
      logical :: opened  ! if a file is opened
      integer :: ios     ! status of the read

      if (present(do_init)) do_init_ = do_init

      if (sfield%read_from_list) then
         ! close opened file
         if (do_init_) then
            counter = 1
            inquire(unit=10, opened=opened)
            if (opened) close(10)
         endif

         ! open filelist for reading
         inquire(unit=10, opened=opened)
         if (.not. opened) then
            open (unit=10, file=trim(sfield%file), iostat=ios)
            if (ios /= 0) write(*,*) 'Could not open file ', trim(sfield%file)
         end if

         ! read a filename
         read (10, *, iostat=ios) filename
         ! check if it is the end of filename
         if (ios/=0) then
            if (present(do_exit)) do_exit = .true.
            close(10)
         else
            if (present(do_exit)) do_exit = .false.
         end if
      else
         if (do_init_) counter = 1
         filename = trim(sfield%file)
         if ((counter > 1) .and. (present(do_exit))) do_exit = .true.
      end if
      counter = counter + 1
   end function get_filename

end module mod_statevector_pdaf
