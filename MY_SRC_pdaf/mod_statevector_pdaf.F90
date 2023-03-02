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
#if defined key_top
  use par_trc, &
       only: jptra
#endif
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
#if defined key_top
     ! Biogeochemistry
     integer, allocatable  :: bgc1(:)
     integer, allocatable  :: bgc2(:)
#endif
  end type field_ids

  ! Declare Fortran type holding the definitions for model fields
  type state_field
     integer :: ndims = 0                  ! Number of field dimensions (2 or 3)
     integer :: dim = 0                    ! Dimension of the field
     integer :: off = 0                    ! Offset of field in state vector
     integer :: jptrc = 0                  ! index of the tracer in nemo tracer variable
     logical :: update = .true.            ! Whether to update this variable 
     character(len=10) :: variable = ''    ! Name of field
     character(len=20) :: name_incr = ''   ! Name of field in increment file
     character(len=20) :: name_rest_n = '' ! Name of field in restart file (n-field)
     character(len=20) :: name_rest_b = '' ! Name of field in restart file (b-field)
     character(len=50) :: file = ''        ! File name stub to read field from
     character(len=50) :: file_state = ''  ! File name to read model state
     character(len=30) :: rst_file = ''    ! Name of restart file
     character(len=20) :: unit = ''        ! Unit of variable
     integer :: transform = 0              ! Type of variable transformation
     real(pwp) :: trafo_shift = 0.0_pwp    ! Constant to shift value in transformation
     integer :: limit = 0                  ! Whether to limit the value of the variable
                     ! 0: no limits, 1: lower limit, 2: upper limit, 3: both limits
     real(pwp) :: max_limit = 0.0_pwp      ! Upper limit of variable
     real(pwp) :: min_limit = 0.0_pwp      ! Lower limit of variable
     real(pwp) :: ensscale = 1.0           ! Scale factor for initial ensemble perturbations
  end type state_field

#if defined key_top
  integer :: n_trc = 0                     !< number of tracer fields
  integer :: n_bgc1 = 0                    !< number of prognostic tracer fields
  integer :: n_bgc2 = 0                    !< number of diagnostic tracer fields
  integer, parameter :: jptra2 = 4         !< number of total diagnostic tracer fields
#endif

  ! Variables to activate a field from the namelist
  logical :: sv_temp = .false. !< Whether to include temperature in state vector
  logical :: sv_salt = .false. !< Whether to include salinity in state vector
  logical :: sv_ssh = .false.  !< Whether to include SSH in state vector
  logical :: sv_uvel = .false. !< Whether to include u-velocity in state vector
  logical :: sv_vvel = .false. !< Whether to include v-velocity in state vector
#if defined key_top
  logical, allocatable :: sv_bgc1(:) !< Whether to include ERGOM in state vector
  logical, allocatable :: sv_bgc2(:) !< Whether to include diagnosed ERGOM variables
#endif

  integer :: id_chl=0          ! Index of chlorophyll field in state vector
  integer :: id_dia=0          ! Index of diatom field in state vector
  integer :: id_fla=0          ! Index of flagellate field in state vector
  integer :: id_cya=0          ! Index of cyanobacteria field in state vector
  integer :: id_netpp=0        ! Index of net primary production in state vector

  !---- The next variables usually do not need editing -----

  integer :: screen=1          ! Verbosity flag

  ! Type variable holding field IDs in state vector
  type(field_ids) :: id

  ! Type variable holding the defintions of model fields
  type(state_field), allocatable :: sfields(:)

  ! Variables to handle multiple fields in the state vector
  integer :: n_fields          !< number of fields in state vector
  integer :: n_fields_covar=0  !< number of fields to read from covariance matrix file

  logical :: update_phys = .true.  !< Whether to update NEMO physics after analysis step
  logical :: update_phyto = .true.  !< Whether to update phytoplankton variables of ERGOM
  logical :: update_nophyto = .true.  !< Whether to update non-phytoplankton variables of ERGOM

contains

!> This routine initializes the array id
!!
  subroutine init_id(nfields)

    implicit none

! *** Arguments ***
    integer, intent(out) :: nfields

! *** Local variables ***
    integer :: cnt               ! Counter
#if defined key_top
    integer :: id_bgc1           ! Counter
    integer :: id_bgc2           ! Counter

    allocate(id%bgc1(jptra))
    allocate(id%bgc2(jptra2))
    id%bgc1(:)=0
    id%bgc2(:)=0

    allocate(sv_bgc1(jptra))
    allocate(sv_bgc2(jptra2))
    sv_bgc1(:) = .false.
    sv_bgc2(:) = .false.
#endif

    ! Namelist to define active parts of state vector
#if defined key_top
    namelist /state_vector/ screen, n_fields_covar, &
         sv_temp, sv_salt, sv_ssh, sv_uvel, sv_vvel, &
         sv_bgc1, sv_bgc2, update_phys, update_phyto, update_nophyto
#else
    namelist /state_vector/ screen, n_fields_covar, &
         sv_temp, sv_salt, sv_ssh, sv_uvel, sv_vvel
#endif

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

#if defined key_top
    do id_bgc1 = 1, jptra
      if (sv_bgc1(id_bgc1)) then
        cnt = cnt + 1
        id%bgc1(id_bgc1) = cnt
        n_bgc1=n_bgc1+1
      end if
    end do

    do id_bgc2 = 1, jptra2
      if (sv_bgc2(id_bgc2)) then
        cnt = cnt + 1
        id%bgc2(id_bgc2) = cnt
        n_bgc2=n_bgc2+1
      end if
    end do
#endif

    ! Set number of fields in state vector
    nfields = cnt
#if defined key_top
    n_trc=n_bgc1+n_bgc2
#endif
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
#if defined key_top
    integer :: id_bgc1           ! Counter
    integer :: id_bgc2           ! Counter
#endif

    namelist /sfields_nml/ sfields
! *** Specifications for each model field in state vector ***

    ! SSH
    id_var = id%ssh
    if (id_var>0) then
       sfields(id_var)%ndims = 2
       sfields(id_var)%dim = sdim2d
       sfields(id_var)%variable = 'SSH_inst'
       sfields(id_var)%name_incr = 'bckineta'
       sfields(id_var)%name_rest_n = 'sshn'
       sfields(id_var)%name_rest_b = 'sshb'
       sfields(id_var)%file = 'files_surf_T.txt'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'm'
       sfields(id_var)%transform = 0
       sfields(id_var)%trafo_shift = 0.0
       sfields(id_var)%limit = 0
    endif

    ! Temperature
    id_var = id%temp
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%dim = sdim3d
       sfields(id_var)%variable = 'votemper'
       sfields(id_var)%name_incr = 'bckint'
       sfields(id_var)%name_rest_n = 'tn'
       sfields(id_var)%name_rest_b = 'tb'
       sfields(id_var)%file = 'files_T.txt'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'degC'
       sfields(id_var)%transform = 0
       sfields(id_var)%trafo_shift = 0.0
    endif

    ! Salinity
    id_var = id%salt
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%dim = sdim3d
       sfields(id_var)%variable = 'vosaline'
       sfields(id_var)%name_incr = 'bckins'
       sfields(id_var)%name_rest_n = 'sn'
       sfields(id_var)%name_rest_b = 'sb'
       sfields(id_var)%file = 'files_T.txt'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'psu'
       sfields(id_var)%transform = 0
       sfields(id_var)%trafo_shift = 0.0
    endif

    ! U-velocity
    id_var = id%uvel
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%dim = sdim3d
       sfields(id_var)%variable = 'uos'
       sfields(id_var)%name_incr = 'bckinu'
       sfields(id_var)%name_rest_n = 'un'
       sfields(id_var)%name_rest_b = 'ub'
       sfields(id_var)%file = 'files_U.txt'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'm/s'
       sfields(id_var)%transform = 0
       sfields(id_var)%trafo_shift = 0.0
    endif

    ! V-velocity
    id_var = id%vvel
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%dim = sdim3d
       sfields(id_var)%variable = 'vos'
       sfields(id_var)%name_incr = 'bckinv'
       sfields(id_var)%name_rest_n = 'vn'
       sfields(id_var)%name_rest_b = 'vb'
       sfields(id_var)%file = 'files_V.txt'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'm/s'
       sfields(id_var)%transform = 0
       sfields(id_var)%trafo_shift = 0.0
    endif

#if defined key_top
    ! BGC
    do id_bgc1 = 1, jptra
      if (sv_bgc1(id_bgc1)) then
        id_var=id%bgc1(id_bgc1)
        sfields(id_var)%ndims = 3
        sfields(id_var)%dim = sdim3d
        sfields(id_var)%jptrc = id_bgc1
        sfields(id_var)%file = 'NORDIC_1d_ERGOM_T_'
        sfields(id_var)%rst_file = 'restart_trc_in.nc'
        sfields(id_var)%transform = 0   ! log-transform
        sfields(id_var)%limit = 1
        sfields(id_var)%min_limit = 0.00001_pwp
        sfields(id_var)%trafo_shift = 0.0
        sfields(id_var)%ensscale = 0.5

        select case (id_bgc1)
        case (1)
          sfields(id_var)%variable = 'NH4'
          sfields(id_var)%name_incr = 'bckinnh4'
          sfields(id_var)%name_rest_n = 'TRNNH4'
          sfields(id_var)%name_rest_b = 'TRBNH4'
          sfields(id_var)%unit = 'mmol m-3'
!        sfields(id_var)%transform = 2   ! log-transform
        case (2)
          sfields(id_var)%variable = 'NO3'
          sfields(id_var)%name_incr = 'bckinno3'
          sfields(id_var)%name_rest_n = 'TRNNO3'
          sfields(id_var)%name_rest_b = 'TRBNO3'
          sfields(id_var)%unit = 'mmol m-3'
        case (3)
          sfields(id_var)%variable = 'PO4'
          !sfields(id_var)%name_incr = 'bckinpo3'
          sfields(id_var)%name_rest_n = 'TRNPO4'
          sfields(id_var)%name_rest_b = 'TRBPO4'
          sfields(id_var)%unit = 'mmol m-3'
        case (4)
          sfields(id_var)%variable = 'SIL'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNSIL'
          sfields(id_var)%name_rest_b = 'TRBSIL'
          sfields(id_var)%unit = 'mmol m-3'
        case (5)
          sfields(id_var)%variable = 'DIA'
          sfields(id_var)%name_incr = 'bckindia'
          sfields(id_var)%name_rest_n = 'TRNDIA'
          sfields(id_var)%name_rest_b = 'TRBDIA'
          sfields(id_var)%unit = 'mmol m-3'
          id_dia = id_var
!          sfields(id_var)%limit = 1
        case (6)
          sfields(id_var)%variable = 'FLA'
          sfields(id_var)%name_incr = 'bckinfla'
          sfields(id_var)%name_rest_n = 'TRNFLA'
          sfields(id_var)%name_rest_b = 'TRBFLA'
          sfields(id_var)%unit = 'mmol m-3'
          id_fla = id_var
        case (7)
          sfields(id_var)%variable = 'CYA'
          sfields(id_var)%name_incr = 'bckincya'
          sfields(id_var)%name_rest_n = 'TRNCYA'
          sfields(id_var)%name_rest_b = 'TRBCYA'
          sfields(id_var)%unit = 'mmol m-3'
          id_cya = id_var
        case (8)
          sfields(id_var)%variable = 'MEZ'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNMEZ'
          sfields(id_var)%name_rest_b = 'TRBMEZ'
          sfields(id_var)%unit = 'mmol m-3'
        case (9)
          sfields(id_var)%variable = 'MIZ'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNMIZ'
          sfields(id_var)%name_rest_b = 'TRBMIZ'
          sfields(id_var)%unit = 'mmol m-3'
        case (10)
          sfields(id_var)%variable = 'DET'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNDET'
          sfields(id_var)%name_rest_b = 'TRBDET'
          sfields(id_var)%unit = 'mmol m-3'
        case (11)
          sfields(id_var)%variable = 'DETs'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNDETs'
          sfields(id_var)%name_rest_b = 'TRBDETs'
          sfields(id_var)%unit = 'mmol m-3'
        case (12)
          sfields(id_var)%variable = 'FE'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNFE'
          sfields(id_var)%name_rest_b = 'TRBFE'
          sfields(id_var)%unit = 'mmol m-3'
        case (13)
          sfields(id_var)%variable = 'LDON'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNLDON'
          sfields(id_var)%name_rest_b = 'TRBLDON'
          sfields(id_var)%unit = 'mmol m-3'
        case (14)
          sfields(id_var)%variable = 'DIC'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNDIC'
          sfields(id_var)%name_rest_b = 'TRBDIC'
          sfields(id_var)%unit = 'mmol m-3'
        case (15)
          sfields(id_var)%variable = 'ALK'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNALK'
          sfields(id_var)%name_rest_b = 'TRBALK'
          sfields(id_var)%unit = 'mmol m-3'
        case (16)
          sfields(id_var)%variable = 'OXY'
          sfields(id_var)%name_incr = 'bckinoxy'
          sfields(id_var)%name_rest_n = 'TRNOXY'
          sfields(id_var)%name_rest_b = 'TRBOXY'
          sfields(id_var)%unit = 'mmol m-3'
          sfields(id_var)%transform = 0     ! Oxy will not be transformed
          sfields(id_var)%limit = 0         ! No limits for Oxy
        end select
      end if
    end do

    ! Diagnostic BGC fields - not part of restart files
    do id_bgc2 = 1, jptra2
      if (sv_bgc2(id_bgc2)) then
        id_var=id%bgc2(id_bgc2)
        sfields(id_var)%ndims = 3
        sfields(id_var)%dim = sdim3d
        sfields(id_var)%jptrc = id_bgc2
        sfields(id_var)%file = 'NORDIC_1d_ERGOM_T_'
        sfields(id_var)%transform = 0
        sfields(id_var)%trafo_shift = 0.0

        select case (id_bgc2)
        case (1)
          sfields(id_var)%variable = 'PCO2'
          sfields(id_var)%unit = 'micro atm'
        case (2)
          sfields(id_var)%variable = 'PH'
          sfields(id_var)%unit = '-'
        case (3)
          id_chl = id_var        ! Store ID of chlorophyll to be used in observation module
          sfields(id_var)%variable = 'CHL'
          sfields(id_var)%unit = 'mg m-3'
          ! Set log-transform if prognostic ERGOM variables are transformed
          if (sfields(id%bgc1(1))%transform==2) then
             sfields(id_var)%transform = 2   ! log-transform
          end if
        case (4)
          id_netpp = id_var      ! Store ID of NETPP to be used in observation module
          sfields(id_var)%variable = 'PP'
          sfields(id_var)%unit = 'microgC m-3* -d'
        end select
      end if
    end do
#endif

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
if (mype==11 .and. task_id==2) write (*,*) 'dims', sfields(:)%dim
if (mype==11 .and. task_id==2) write (*,*) 'offs', sfields(:)%off

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
