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

  use mod_kind_pdaf
  implicit none
  save

  !---- `field_ids` and `state_field` need to be adapted for a DA case -----

  ! Declare Fortran type holding the indices of model fields in the state vector
  ! This can be extended to any number of fields - it serves to give each field a name
  type field_ids
     ! Ocean Physics
     integer :: ssh = 0
     integer :: temp = 0
     integer :: salt = 0
     integer :: uvel = 0
     integer :: vvel = 0
  end type field_ids

  ! Declare Fortran type holding the definitions for model fields
  type state_field
     integer :: ndims = 0                  ! Number of field dimensions (2 or 3)
     integer :: dim = 0                    ! Dimension of the field
     integer :: off = 0                    ! Offset of field in state vector
     character(len=10) :: variable = ''    ! Name of field
     character(len=20) :: name_incr = ''   ! Name of field in increment file
     character(len=20) :: name_rest_n = '' ! Name of field in restart file (n-field)
     character(len=20) :: name_rest_b = '' ! Name of field in restart file (b-field)
     character(len=50) :: file = ''        ! File name stub to read field from
     character(len=50) :: file_state = ''  ! File name stub to read field from
     character(len=30) :: rst_file = ''    ! Name of restart file
     character(len=20) :: unit = ''        ! Unit of variable
     integer :: transform = 0              ! Type of variable transformation
     real(pwp) :: trafo_shift = 0.0_pwp    ! Constant to shift value in transformation
     integer :: limit = 0                  ! Whether to limit the value of the variable
                     ! 0: no limits, 1: lower limit, 2: upper limit, 3: both limits
     real(pwp) :: max_limit = 0.0_pwp      ! Upper limit of variable
     real(pwp) :: min_limit = 0.0_pwp      ! Lower limit of variable
  end type state_field


  !---- The next variables usually do not need editing -----

  integer :: screen=1          ! Verbosity flag

  ! Type variable holding field IDs in state vector
  type(field_ids) :: id

  ! Type variable holding the defintions of model fields
  type(state_field), allocatable :: sfields(:)

  ! Variables to handle multiple fields in the state vector
  integer :: n_fields      !< number of fields in state vector

contains

!> This routine initializes the array id
!!
  subroutine init_id(nfields)

    implicit none

! *** Arguments ***
    integer, intent(out) :: nfields

! *** Local variables *** 
    integer :: cnt               ! Counter

    ! Variables to activate a field from the namelist
    logical :: sv_temp = .false. ! Whether to include temperature in state vector
    logical :: sv_salt = .false. ! Whether to include salinity in state vector
    logical :: sv_ssh = .false.  ! Whether to include SSH in state vector
    logical :: sv_uvel = .false. ! Whether to include u-velocity in state vector
    logical :: sv_vvel = .false. ! Whether to include v-velocity in state vector

    ! Namelist to define active parts of state vector
    namelist /state_vector/ screen, &
         sv_temp, sv_salt, sv_ssh, sv_uvel, sv_vvel


! **********************
! *** Initialization ***
! **********************

! *** Read namelist file for state vector setup

    open (500,file='namelist_cfg.pdaf')
    read (500,NML=state_vector)
    close (500)

! *** Now setup field indices in state vector

    cnt = 0
    if (sv_ssh) then
       cnt = cnt + 1
       id%ssh = cnt
    end if

    if (sv_temp) then
       cnt = cnt + 1
       id%temp = cnt
    end if

    if (sv_salt) then
       cnt = cnt + 1
       id%salt = cnt
    end if

    if (sv_uvel) then
       cnt = cnt + 1
       id%uvel = cnt
    end if

    if (sv_vvel) then
       cnt = cnt + 1
       id%vvel = cnt
    end if

    ! Set number of fields in state vector
    nfields = cnt

  end subroutine init_id
! ===================================================================================

!> This initializes the array sfields
!!
!! This routine initializes the sfields array with specifications
!! of the fields in the state vector.
  subroutine init_sfields()

    use mod_kind_pdaf
    use mod_nemo_pdaf, &
         only: sdim2d, sdim3d

    implicit none

! *** Local variables *** 
    integer :: id_var            ! Index of a variable in state vector

    namelist /sfields_nml/ sfields
! *** Specifications for each model field in state vector ***

    ! SSH
    id_var = id%ssh
    if (id_var>0) then
       sfields(id_var)%ndims = 2
       sfields(id_var)%dim = sdim2d
       sfields(id_var)%variable = 'zos'
       sfields(id_var)%name_incr = 'bckineta'
       sfields(id_var)%name_rest_n = 'sshn'
       sfields(id_var)%name_rest_b = 'sshb'
       sfields(id_var)%file = 'ORCA2_1d_'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'm'
    endif

    ! Temperature
    id_var = id%temp
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%dim = sdim3d
       sfields(id_var)%variable = 'thetao'
       sfields(id_var)%name_incr = 'bckint'
       sfields(id_var)%name_rest_n = 'tn'
       sfields(id_var)%name_rest_b = 'tb'
       sfields(id_var)%file = 'ORCA2_1d_'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'degC'
    endif

    ! Salinity
    id_var = id%salt
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%dim = sdim3d
       sfields(id_var)%variable = 'so'
       sfields(id_var)%name_incr = 'bckins'
       sfields(id_var)%name_rest_n = 'sn'
       sfields(id_var)%name_rest_b = 'sb'
       sfields(id_var)%file = 'ORCA2_1d_'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = '1e-3'
       sfields(id_var)%limit = 1             ! Apply lower limit
       sfields(id_var)%min_limit = 0.0_pwp   ! Salinity is never negative
    endif

    ! U-velocity
    id_var = id%uvel
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%dim = sdim3d
       sfields(id_var)%variable = 'uo'
       sfields(id_var)%name_incr = 'bckinu'
       sfields(id_var)%name_rest_n = 'un'
       sfields(id_var)%name_rest_b = 'ub'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'm/s'
    endif

    ! V-velocity
    id_var = id%vvel
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%dim = sdim3d
       sfields(id_var)%variable = 'vo'
       sfields(id_var)%name_incr = 'bckinv'
       sfields(id_var)%name_rest_n = 'vn'
       sfields(id_var)%name_rest_b = 'vb'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'm/s'
    endif

    open (500,file='namelist_cfg.pdaf')
    read (500,NML=sfields_nml)
    close (500)

    do id_var = 1, n_fields
      if (sfields(id_var)%ndims == 2) then
        sfields(id_var)%dim = sdim2d
      else if (sfields(id_var)%ndims == 3) then
        sfields(id_var)%dim = sdim3d
      else
        write (*, '(a,i2,a)') 'NEMO-PDAF: cannot handle', sfields(id_var)%ndims, ' number of dimensions.'
      end if
    end do
  end subroutine init_sfields
! ===================================================================================

!> Calculate the dimension of the process-local statevector.
!!
!! This routine is generic. case-specific adaptions should only
!! by done in the routines init_id and init_sfields.
!!
  subroutine setup_statevector(dim_state, dim_state_p)

    use mod_kind_pdaf
    use mod_parallel_pdaf, &
         only: mype=>mype_ens, npes=>npes_ens, task_id, comm_ensemble, &
         comm_model, MPI_SUM, MPI_INTEGER, MPIerr

    implicit none

! *** Arguments ***
    integer, intent(out) :: dim_state    !< Global dimension of state vector
    integer, intent(out) :: dim_state_p  !< Local dimension of state vector

! *** Local variables *** 
    integer :: i                 ! Counters


! ***********************************
! *** Initialize the state vector ***
! ***********************************

! *** Initialize array `id` ***

    call init_id(n_fields)

! *** Initialize array `sfields` ***

    allocate(sfields(n_fields))

    call init_sfields()

! *** Compute offsets ***

    ! Define offsets in state vector
    sfields(1)%off = 0
    do i = 2, n_fields
       sfields(i)%off = sfields(i-1)%off + sfields(i-1)%dim
    end do

! *** Set state vector dimension ***

    dim_state_p = sum(sfields(:)%dim)

! *** Write information about the state vector ***

    if (mype==0) then 
       write (*,'(/a,2x,a)') 'NEMO-PDAF', '*** Setup of state vector ***' 
       write (*,'(a,5x,a,i5)') 'NEMO-PDAF', '--- Number of fields in state vector:', n_fields
       if (npes==1) then
          write (*,'(a,5x,a2,3x,a8,6x,a5,7x,a3,7x,a6)') 'NEMO-PDAF','ID', 'variable', 'ndims', 'dim', 'offset'
          write (*,'(a,5x,49a)') 'NEMO-PDAF', ('-',i=1,49)
          do i = 1, n_fields
             write (*,'(a, 2x,i5,3x,a10,2x,i5,3x,i10,3x,i10)') &
                  'NEMO-PDAF',i, sfields(i)%variable, sfields(i)%ndims, sfields(i)%dim, sfields(i)%off
          end do
       else
          if (task_id==1) then
             write (*,'(a,a4,7x,a2,3x,a8,6x,a5,7x,a3,7x,a6)') 'NEMO-PDAF','pe','ID', 'variable', 'ndims', 'dim', 'offset'
          end if
       end if
    end if

    if (npes>1 .and. (mype==0 .or. (task_id==1 .and. screen>2))) then
       do i = 1, n_fields
          write (*,'(a,i3,5x,i5,3x,a10,2x,i5,3x,i10,3x,i10)') &
               'NEMO-PDAF',mype, i, sfields(i)%variable, sfields(i)%ndims, sfields(i)%dim, sfields(i)%off
       end do
    end if

    if (npes==1) then
       write (*,'(a,2x,a,1x,i10)') 'NEMO-PDAF', 'Full state dimension: ',dim_state_p
    else
       if (task_id==1) then
          if (screen>1 .or. mype==0) &
               write (*,'(a,2x,a,1x,i4,2x,a,1x,i10)') &
               'NEMO-PDAF', 'PE', mype, 'PE-local full state dimension: ',dim_state_p

          call MPI_Reduce(dim_state_p, dim_state, 1, MPI_INTEGER, MPI_SUM, 0, COMM_model, MPIerr)
          if (mype==0) then
             write (*,'(a,2x,a,1x,i10)') 'NEMO-PDAF', 'Global state dimension: ',dim_state
          end if
       end if
    end if
    call MPI_Barrier(comm_ensemble, MPIerr)

  end subroutine setup_statevector

end module mod_statevector_pdaf
