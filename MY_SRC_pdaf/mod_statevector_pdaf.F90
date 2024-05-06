!> Building the state vector
!!
!! This module provides variables & routines for
!! defining the state vector.
!!
!! The module contains three routines
!! - **init_id** - initialize the array `id`
!! - **init_sfields** - initialize the array `sfields`
!! - **setup_statevector** - generic routine controlling the initialization
!!
!! The declarations of **id** and **sfields** as well as the
!! routines **init_id** and **init_sfields** might need to be
!! adapted to a particular modeling case. However, for most
!! parts also the configruation using the namelist is possible.
!!
module mod_statevector_pdaf

  use mod_kind_pdaf

  implicit none
  save

  !---- `field_ids` and `state_field` can be adapted for a DA case -----

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
     integer, allocatable  :: bgc_prog(:)
     integer, allocatable  :: bgc_diag(:)
#endif
  end type field_ids

  !< Declare Fortran type holding the definitions for model fields
  type state_field
     integer :: ndims = 0                  !< Number of field dimensions (2 or 3)
     integer :: dim = 0                    !< Dimension of the field
     integer :: off = 0                    !< Offset of field in state vector
     integer :: jptrc = 0                  !< index of the tracer in nemo tracer variable
     logical :: update = .false.           !< Whether to update this variable in the analysis step
     logical :: vloc = .false.             !< Whether to apply vertical localization
     real(pwp) :: vloc_limit = 100.0       !< Limit depth for vertical localization
     character(len=3) :: type = ''         !< Type of field (phy/bio)
     character(len=10) :: variable = ''    !< Name of field
     character(len=20) :: name_incr = ''   !< Name of field in increment file
     character(len=20) :: name_rest_n = '' !< Name of field in restart file (n-field)
     character(len=20) :: name_rest_b = '' !< Name of field in restart file (b-field)
     character(len=50) :: file = ''        !< File name stub to read field from
     character(len=50) :: file_state = ''  !< File name to read model state
     character(len=30) :: rst_file = ''    !< Name of restart file
     character(len=20) :: unit = ''        !< Unit of variable
     integer :: transform = 0              !< Type of variable transformation
     real(pwp) :: trafo_shift = 0.0_pwp    !< Constant to shift value in transformation
     integer :: limit = 0                  !< Whether to limit the value of the variable
                     !< 0: no limits, 1: lower limit, 2: upper limit, 3: both limits
     real(pwp) :: max_limit = 0.0_pwp      !< Upper limit of variable
     real(pwp) :: min_limit = 0.0_pwp      !< Lower limit of variable
     real(pwp) :: ensscale = 1.0           !< Scale factor for initial ensemble perturbations
  end type state_field

  ! Declare Fortran type holding the definitions for local model fields
  ! This is separate from state_field to support OpenMP
  type state_field_l
     integer :: dim = 0                    !< Dimension of the field
     integer :: off = 0                    !< Offset of field in state vector
  end type state_field_l

  ! Variables to activate a field from the namelist
  logical :: sv_temp = .false.             !< Whether to include temperature in state vector
  logical :: sv_salt = .false.             !< Whether to include salinity in state vector
  logical :: sv_ssh  = .false.             !< Whether to include SSH in state vector
  logical :: sv_uvel = .false.             !< Whether to include u-velocity in state vector
  logical :: sv_vvel = .false.             !< Whether to include v-velocity in state vector

  ! Control whether to apply assimilation increment to a variable
  logical :: update_temp  = .false.     !< Whether to update NEMO physics after analysis step
  logical :: update_salt  = .false.     !< Whether to update NEMO physics after analysis step
  logical :: update_vel   = .false.     !< Whether to update NEMO physics after analysis step
  logical :: update_ssh   = .false.     !< Whether to update NEMO physics after analysis step

  ! Variables for biogeochemistry
  integer :: n_trc = 0                     !< number of tracer fields
  integer :: n_bgc_prog = 0                !< number of active prognostic tracer fields
  integer :: n_bgc_diag = 0                !< number of active diagnostic tracer fields
  integer :: jpbgc_prog = 0                !< number of total prognistic tracer fields
  integer :: jpbgc_diag = 4                !< number of total diagnostic tracer fields

  ! Variables to activate a field from the namelist
  logical, allocatable :: sv_bgc_prog(:) !< Whether to include BGC in state vector
  logical, allocatable :: sv_bgc_diag(:) !< Whether to include diagnostic BGC variables

  ! Control whether to apply assimilation increment to a variable group
  logical :: update_phyto = .false.     !< Whether to update phytoplankton variables of ERGOM (DIA, FLA, CYA)
  logical :: update_zoo   = .false.     !< Whether to update zooplankton variables of ERGOM (MIZ, MEZ)
  logical :: update_det   = .false.     !< Whether to update detritus variables of ERGOM (DET, DETs)
  logical :: update_nut   = .false.     !< Whether to update nutrient variables of ERGOM (NH4, NO3, PO4, FE)
  logical :: update_other = .false.     !< Whether to update non-phytoplankton variables of ERGOM
  logical :: update_diag  = .false.     !< Whether to update diagnostic variables of ERGOM

  ! Control whether a single BGC variable is updated
  logical :: update_NH4  = .false.     !< Whether to update variable NH4
  logical :: update_NO3  = .false.     !< Whether to update variable NO3
  logical :: update_PO4  = .false.     !< Whether to update variable PO4
  logical :: update_SIL  = .false.     !< Whether to update variable SIL
  logical :: update_DIA  = .false.     !< Whether to update variable DIA
  logical :: update_FLA  = .false.     !< Whether to update variable FLA
  logical :: update_CYA  = .false.     !< Whether to update variable CYA
  logical :: update_MEZ  = .false.     !< Whether to update variable MEZ
  logical :: update_MIZ  = .false.     !< Whether to update variable MIZ
  logical :: update_DET1 = .false.     !< Whether to update variable DET
  logical :: update_DETs = .false.     !< Whether to update variable DETs
  logical :: update_FE   = .false.     !< Whether to update variable FE
  logical :: update_LDON = .false.     !< Whether to update variable LDON
  logical :: update_DIC  = .false.     !< Whether to update variable DIC
  logical :: update_ALK  = .false.     !< Whether to update variable ALK
  logical :: update_OXY  = .false.     !< Whether to update variable OXY
  logical :: update_CHL  = .false.     !< Whether to update diagnostic variable CHL
  logical :: update_PCO2 = .false.     !< Whether to update diagnostic variable PCO2
  logical :: update_PH   = .false.     !< Whether to update diagnostic variable PH
  logical :: update_PP   = .false.     !< Whether to update diagnostic variable PP

  ! Control vertical localization in groups
  logical :: vloc_phys = .false.       !< Whether physics variables use vertical localization
  logical :: vloc_bgc  = .false.       !< Whether ERGOM variables use vertical localization
  real(pwp) :: vloc_depth_phys         !< vertcial localization depth of physics variables
  real(pwp) :: vloc_depth_bgc          !< vertcial localization depth of ERGOM variables


  ! Helper variables to point to particular fields
  integer :: id_chl=0          !< Index of chlorophyll field in state vector
  integer :: id_dia=0          !< Index of diatom field in state vector
  integer :: id_fla=0          !< Index of flagellate field in state vector
  integer :: id_cya=0          !< Index of cyanobacteria field in state vector
  integer :: id_netpp=0        !< Index of net primary production in state vector


  !---- The next variables usually do not need editing -----

  integer :: screen=1          ! Verbosity flag

  ! Type variable holding field IDs in state vector
  type(field_ids) :: id

  ! Type variable holding the definitions of model fields
  type(state_field), allocatable :: sfields(:)

  ! Type variable holding the definitions of local model fields
  ! This is separate from sfields to support OpenMP
  type(state_field_l), allocatable :: sfields_l(:)

!$OMP THREADPRIVATE(sfields_l)

  ! Variables to handle multiple fields in the state vector
  integer :: n_fields          !< number of fields in state vector
  integer :: n_fields_covar=0  !< number of fields to read from covariance matrix file

contains

!> This routine initializes the array id
!!
  subroutine init_id(nfields)

#if defined key_top
  use mod_nemo_pdaf, &
       only: jptra
#endif

    implicit none

! *** Arguments ***
    integer, intent(out) :: nfields

! *** Local variables ***
    integer :: cnt               ! Counter
    integer :: id_bgc_prog       ! Counter for prognostic BGC variables
    integer :: id_bgc_diag       ! Counter for diagnostic BGC variables


! **********************
! *** Initialization ***
! **********************

#if defined key_top
    ! Set total number of prognostic and diagnostic BGC fields
    jpbgc_prog = jptra           ! Number of prognostic BGC fields
    jpbgc_diag = 4               ! Number of diagnostic BGC fields

    ! Prepare arrays for indices and switches for BGC fields
    allocate(id%bgc_prog(jpbgc_prog))
    allocate(id%bgc_diag(jpbgc_diag))
    id%bgc_prog(:)=0
    id%bgc_diag(:)=0

    allocate(sv_bgc_prog(jpbgc_prog))
    allocate(sv_bgc_diag(jpbgc_diag))
    sv_bgc_prog(:) = .false.
    sv_bgc_diag(:) = .false.
#endif

! *** Read namelist file for state vector setup

#if defined key_top
    namelist /state_vector/ screen, n_fields_covar, &
         sv_temp, sv_salt, sv_ssh, sv_uvel, sv_vvel, &
         sv_bgc_prog, sv_bgc_diag
#else
    namelist /state_vector/ screen, n_fields_covar, &
         sv_temp, sv_salt, sv_ssh, sv_uvel, sv_vvel
#endif

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
    do id_bgc_prog = 1, jpbgc_prog
       if (sv_bgc_prog(id_bgc_prog)) then
          cnt = cnt + 1
          id%bgc_prog(id_bgc_prog) = cnt
          n_bgc_prog=n_bgc_prog+1
       end if
    end do

    do id_bgc_diag = 1, jpbgc_diag
       if (sv_bgc_diag(id_bgc_diag)) then
          cnt = cnt + 1
          id%bgc_diag(id_bgc_diag) = cnt
          n_bgc_diag=n_bgc_diag+1
       end if
    end do

    ! Set number of BGC fields in state vector
    n_trc=n_bgc_prog+n_bgc_diag
#endif

    ! Set total number of fields in state vector
    nfields = cnt

  end subroutine init_id
! ===================================================================================

!> This initializes the array sfields
!!
!! This routine initializes the sfields array with specifications
!! of the fields in the state vector.
!!
  subroutine init_sfields()

    use mod_kind_pdaf
    use mod_nemo_pdaf, &
         only: sdim2d, sdim3d

    implicit none

! *** Local variables ***
    integer :: id_var            ! Index of a variable in state vector
#if defined key_top
    integer :: id_bgc_prog       ! Counter
    integer :: id_bgc_diag       ! Counter
#endif

    namelist /sfields_nml/ sfields
! *** Specifications for each model field in state vector ***

    ! SSH
    id_var = id%ssh
    if (id_var>0) then
       sfields(id_var)%ndims = 2
       sfields(id_var)%variable = 'SSH_inst'
       sfields(id_var)%name_incr = 'bckineta'
       sfields(id_var)%name_rest_n = 'sshn'
       sfields(id_var)%name_rest_b = 'sshb'
       sfields(id_var)%file = 'files_surf_T.txt'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'm'
       sfields(id_var)%type = 'phy'
       sfields(id_var)%transform = 0
       sfields(id_var)%trafo_shift = 0.0
       sfields(id_var)%limit = 0
       if (update_ssh) sfields(id_var)%update = .true.
    endif

    ! Temperature
    id_var = id%temp
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%variable = 'votemper'
       sfields(id_var)%name_incr = 'bckint'
       sfields(id_var)%name_rest_n = 'tn'
       sfields(id_var)%name_rest_b = 'tb'
       sfields(id_var)%file = 'files_T.txt'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'degC'
       sfields(id_var)%type = 'phy'
       sfields(id_var)%transform = 0
       sfields(id_var)%trafo_shift = 0.0
       if (update_temp) sfields(id_var)%update = .true.
       if (vloc_phys) then
          sfields(id_var)%vloc = .true.
          sfields(id_var)%vloc_limit = vloc_depth_phys
       end if
    endif

    ! Salinity
    id_var = id%salt
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%variable = 'vosaline'
       sfields(id_var)%name_incr = 'bckins'
       sfields(id_var)%name_rest_n = 'sn'
       sfields(id_var)%name_rest_b = 'sb'
       sfields(id_var)%file = 'files_T.txt'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'psu'
       sfields(id_var)%type = 'phy'
       sfields(id_var)%transform = 0
       sfields(id_var)%trafo_shift = 0.0
       sfields(id_var)%limit = 0
       sfields(id_var)%min_limit = 0.000001
       sfields(id_var)%max_limit = 36.0
       if (update_salt) sfields(id_var)%update = .true.
       if (vloc_phys) then
          sfields(id_var)%vloc = .true.
          sfields(id_var)%vloc_limit = vloc_depth_phys
       end if
    endif

    ! U-velocity
    id_var = id%uvel
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%variable = 'uos'
       sfields(id_var)%name_incr = 'bckinu'
       sfields(id_var)%name_rest_n = 'un'
       sfields(id_var)%name_rest_b = 'ub'
       sfields(id_var)%file = 'files_U.txt'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'm/s'
       sfields(id_var)%type = 'phy'
       sfields(id_var)%transform = 0
       sfields(id_var)%trafo_shift = 0.0
       if (update_vel) sfields(id_var)%update = .true.
       if (vloc_phys) then
          sfields(id_var)%vloc = .true.
          sfields(id_var)%vloc_limit = vloc_depth_phys
       end if
    endif

    ! V-velocity
    id_var = id%vvel
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%variable = 'vos'
       sfields(id_var)%name_incr = 'bckinv'
       sfields(id_var)%name_rest_n = 'vn'
       sfields(id_var)%name_rest_b = 'vb'
       sfields(id_var)%file = 'files_V.txt'
       sfields(id_var)%rst_file = 'restart_in.nc'
       sfields(id_var)%unit = 'm/s'
       sfields(id_var)%type = 'phy'
       sfields(id_var)%transform = 0
       sfields(id_var)%trafo_shift = 0.0
       if (update_vel) sfields(id_var)%update = .true.
       if (vloc_phys) then
          sfields(id_var)%vloc = .true.
          sfields(id_var)%vloc_limit = vloc_depth_phys
       end if
    endif

#if defined key_top

    ! Activate variable groups 
    if (update_phyto) then
       update_DIA = .true.
       update_FLA = .true.
       update_CYA = .true.
    end if

    if (update_zoo) then
       update_MEZ = .true.
       update_MIZ = .true.
    end if

    if (update_det) then
       update_DET1 = .true.
       update_DETs = .true.
    end if

    if (update_nut) then
       update_NH4 = .true.
       update_NO3 = .true.
       update_PO4 = .true.
       update_SIL = .true.
       update_FE  = .true.
    end if

    if (update_other) then
       update_LDON = .true.
       update_DIC  = .true.
       update_ALK  = .true.
    end if

    if (update_diag) then
       update_CHL = .true.
       update_PCO2 = .true.
       update_PH = .true.
       update_PP = .true.
    end if

    ! BGC
    do id_bgc_prog = 1, jpbgc_prog
      if (sv_bgc_prog(id_bgc_prog)) then
        id_var=id%bgc_prog(id_bgc_prog)
        sfields(id_var)%ndims = 3
        sfields(id_var)%jptrc = id_bgc_prog
        sfields(id_var)%file = 'NORDIC_1d_ERGOM_T_'
        sfields(id_var)%rst_file = 'restart_trc_in.nc'
        sfields(id_var)%type = 'bio'
        sfields(id_var)%transform = 0   ! log-transform
        sfields(id_var)%limit = 1
        sfields(id_var)%min_limit = 0.00000001_pwp
        sfields(id_var)%trafo_shift = 0.0
        sfields(id_var)%ensscale = 0.5
        if (vloc_bgc) then
           sfields(id_var)%vloc = .true.
           sfields(id_var)%vloc_limit = vloc_depth_bgc
        end if

        select case (id_bgc_prog)
        case (1)
          sfields(id_var)%variable = 'NH4'
          sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNNH4'
          sfields(id_var)%name_rest_b = 'TRBNH4'
          sfields(id_var)%unit = 'mmol m-3'
          if (update_NH4) sfields(id_var)%update = .true.
        case (2)
          sfields(id_var)%variable = 'NO3'
          sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNNO3'
          sfields(id_var)%name_rest_b = 'TRBNO3'
          sfields(id_var)%unit = 'mmol m-3'
          if (update_NO3) sfields(id_var)%update = .true.
        case (3)
          sfields(id_var)%variable = 'PO4'
          sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNPO4'
          sfields(id_var)%name_rest_b = 'TRBPO4'
          sfields(id_var)%unit = 'mmol m-3'
          if (update_PO4) sfields(id_var)%update = .true.
        case (4)
          sfields(id_var)%variable = 'SIL'
          sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNSIL'
          sfields(id_var)%name_rest_b = 'TRBSIL'
          sfields(id_var)%unit = 'mmol m-3'
          if (update_SIL) sfields(id_var)%update = .true.
        case (5)
          sfields(id_var)%variable = 'DIA'
          sfields(id_var)%name_incr = 'bckindia'
          sfields(id_var)%name_rest_n = 'TRNDIA'
          sfields(id_var)%name_rest_b = 'TRBDIA'
          sfields(id_var)%unit = 'mmol m-3'
          id_dia = id_var
          if (update_DIA) sfields(id_var)%update = .true.
        case (6)
          sfields(id_var)%variable = 'FLA'
          sfields(id_var)%name_incr = 'bckinfla'
          sfields(id_var)%name_rest_n = 'TRNFLA'
          sfields(id_var)%name_rest_b = 'TRBFLA'
          sfields(id_var)%unit = 'mmol m-3'
          id_fla = id_var
          if (update_FLA) sfields(id_var)%update = .true.
        case (7)
          sfields(id_var)%variable = 'CYA'
          sfields(id_var)%name_incr = 'bckincya'
          sfields(id_var)%name_rest_n = 'TRNCYA'
          sfields(id_var)%name_rest_b = 'TRBCYA'
          sfields(id_var)%unit = 'mmol m-3'
          id_cya = id_var
          if (update_CYA) sfields(id_var)%update = .true.
        case (8)
          sfields(id_var)%variable = 'MEZ'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNMEZ'
          sfields(id_var)%name_rest_b = 'TRBMEZ'
          sfields(id_var)%unit = 'mmol m-3'
          if (update_MEZ) sfields(id_var)%update = .true.
        case (9)
          sfields(id_var)%variable = 'MIZ'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNMIZ'
          sfields(id_var)%name_rest_b = 'TRBMIZ'
          sfields(id_var)%unit = 'mmol m-3'
          if (update_MIZ) sfields(id_var)%update = .true.
        case (10)
          sfields(id_var)%variable = 'DET'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNDET'
          sfields(id_var)%name_rest_b = 'TRBDET'
          sfields(id_var)%unit = 'mmol m-3'
          if (update_DET1) sfields(id_var)%update = .true.
        case (11)
          sfields(id_var)%variable = 'DETs'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNDETs'
          sfields(id_var)%name_rest_b = 'TRBDETs'
          sfields(id_var)%unit = 'mmol m-3'
          if (update_DETs) sfields(id_var)%update = .true.
        case (12)
          sfields(id_var)%variable = 'FE'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNFE'
          sfields(id_var)%name_rest_b = 'TRBFE'
          sfields(id_var)%unit = 'mmol m-3'
          if (update_FE) sfields(id_var)%update = .true.
        case (13)
          sfields(id_var)%variable = 'LDON'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNLDON'
          sfields(id_var)%name_rest_b = 'TRBLDON'
          sfields(id_var)%unit = 'mmol m-3'
          if (update_LDON) sfields(id_var)%update = .true.
        case (14)
          sfields(id_var)%variable = 'DIC'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNDIC'
          sfields(id_var)%name_rest_b = 'TRBDIC'
          sfields(id_var)%unit = 'mmol m-3'
          if (update_DIC) sfields(id_var)%update = .true.
        case (15)
          sfields(id_var)%variable = 'ALK'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNALK'
          sfields(id_var)%name_rest_b = 'TRBALK'
          sfields(id_var)%unit = 'mmol m-3'
          if (update_ALK) sfields(id_var)%update = .true.
        case (16)
          sfields(id_var)%variable = 'OXY'
          sfields(id_var)%name_incr = 'bckinoxy'
          sfields(id_var)%name_rest_n = 'TRNOXY'
          sfields(id_var)%name_rest_b = 'TRBOXY'
          sfields(id_var)%unit = 'mmol m-3'
          sfields(id_var)%transform = 0     ! Oxy will not be transformed
          sfields(id_var)%limit = 0         ! No limits for Oxy
          if (update_OXY) sfields(id_var)%update = .true.
        end select
      end if
    end do

    ! Diagnostic BGC fields - not part of restart files
    do id_bgc_diag = 1, jpbgc_diag

      if (sv_bgc_diag(id_bgc_diag)) then
        id_var=id%bgc_diag(id_bgc_diag)
        sfields(id_var)%ndims = 3
        sfields(id_var)%jptrc = id_bgc_diag
        sfields(id_var)%file = 'NORDIC_1d_ERGOM_T_'
        sfields(id_var)%type = 'bio'
        sfields(id_var)%transform = 0
        sfields(id_var)%trafo_shift = 0.0
        if (vloc_bgc) then
           sfields(id_var)%vloc = .true.
           sfields(id_var)%vloc_limit = vloc_depth_bgc
        end if

        select case (id_bgc_diag)
        case (1)
           sfields(id_var)%variable = 'PCO2'
           sfields(id_var)%unit = 'micro atm'
           if (update_PCO2) sfields(id_var)%update = .true.
        case (2)
           sfields(id_var)%variable = 'PH'
           sfields(id_var)%unit = '-'
           if (update_PH) sfields(id_var)%update = .true.
        case (3)
           id_chl = id_var        ! Store ID of chlorophyll to be used in observation module
           sfields(id_var)%variable = 'CHL'
           sfields(id_var)%unit = 'mg m-3'
           ! Set log-transform if prognostic chlorophyll are transformed
           if (sfields(id%bgc_prog(1))%transform==2) sfields(id_var)%transform = 2   ! log-transform
           if (update_CHL) sfields(id_var)%update = .true.
        case (4)
           id_netpp = id_var      ! Store ID of NETPP to be used in observation module
           sfields(id_var)%variable = 'PP'
           sfields(id_var)%unit = 'microgC m-3* -d'
           if (update_PP) sfields(id_var)%update = .true.
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

! *** Set state vector dimension ***

    dim_state_p = sum(sfields(:)%dim)

! *** Write information about the state vector ***

    if (mype==0) then
       write (*,'(/a,2x,a)') 'NEMO-PDAF', '*** Setup of state vector ***'
       write (*,'(a,5x,a,i5)') 'NEMO-PDAF', '--- Number of fields in state vector:', n_fields
       write (*,'(a,a4,3x,a2,3x,a8,6x,a5,7x,a3,7x,a6,4x,a6,a6)') &
            'NEMO-PDAF','pe','ID', 'variable', 'ndims', 'dim', 'offset', 'update', 'vloc'
    end if

    if (mype==0 .or. (task_id==1 .and. screen>2)) then
       do i = 1, n_fields
          write (*,'(a, i4, i5,3x,a10,2x,i5,3x,i10,3x,i10,4x,l,6x,l)') 'NEMO-PDAF', &
               mype, i, sfields(i)%variable, sfields(i)%ndims, sfields(i)%dim, sfields(i)%off, &
               sfields(i)%update, sfields(i)%vloc
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
