!> Initialise PDAF
!!
!! This modules contains the initialisation routine for PDAF
!! `init_pdaf`. Here the ensemble is initialised and distributed
!! and the statevector and state variable information are computed.
!! 
module mod_init_pdaf

   implicit none

contains

  !> This routine collects the initialization of variables for PDAF.
  !!
  !! The initialization routine `PDAF_init` is called
  !! such that the internal initialization of PDAF is performed.
  !! The initialization is used to set-up local domain and filter options
  !! such as the filter type, inflation, and localization radius.
  !! This variant is for the online mode of PDAF.
  !!
  !! The ensemble is initialised in `init_ens_pdaf`, and is then
  !! distributed to the model in `distribute_state_pdaf`. The arrays
  !! for the incremental analysis update (IAU) are initialised in
  !! `asm_inc_init_pdaf`.
  !!
  !! The NEMO grid information, e.g. on wet surface points is 
  !! initialized in `set_nemo_grid`.
  !!
  !! The statevector dimension, and the offsets and dimensions of the
  !! statevector variables are calculated in `setup_state`.
  !!
  !! Much of the initialisation is read from a PDAF-specific namelist.
  !! This is performed in `read_config_pdaf`.
  !!
  !! **Calling Sequence**
  !!
  !! - Called from: `nemogcm.F90`
  !!
  !! - Calls: `setup_state`
  !!          `set_nemo_grid`
  !!          `read_config_pdaf`
  !!          `init_pdaf_info`
  !!          `PDAF_set_comm_pdaf`
  !!          `PDAF_init`
  !!          `asm_inc_init_pdaf`
  !!          `PDAF_get_state`
  !!
  subroutine init_pdaf()

    use mod_kind_pdaf
    use mod_parallel_pdaf, &
         only: n_modeltasks, task_id, COMM_model, COMM_filter, &
         COMM_couple, COMM_ensemble, mype_ens, filterpe, abort_parallel
    use mod_assimilation_pdaf, &
         only: dim_state, dim_state_p, screen, filtertype, subtype, dim_ens, &
         incremental, type_forget, forget, rank_analysis_enkf, &
         type_trans, type_sqrt, delt_obs, locweight
    use mod_iau_pdaf, &
         only: asm_inc_init_pdaf
    use mod_nemo_pdaf, &
         only: set_nemo_grid
    use mod_statevector_pdaf, &
         only: setup_statevector
    use mod_util_pdaf, &
         only: init_info_pdaf, read_config_pdaf
    use mod_obs_sst_cmems_pdafomi, &
         only: assim_sst_cmems, rms_obs_sst_cmems, &
         lradius_sst_cmems, sradius_sst_cmems, mode_sst_cmems, dist_sst_cmems
    use mod_obs_ssh_mgrid_pdafomi, &
         only: assim_ssh_mgrid, rms_ssh_mgrid, &
         lradius_ssh_mgrid, sradius_ssh_mgrid

    implicit none

! *** Local variables
    integer :: filter_param_i(7) ! Integer parameter array for filter
    real    :: filter_param_r(2) ! Real parameter array for filter
    integer :: status_pdaf       ! PDAF status flag
    integer :: doexit, steps     ! Not used in this implementation
    real(pwp) :: timenow         ! Not used in this implementation
      
    external :: init_ens_pdaf         ! Ensemble initialization
    external :: next_observation_pdaf ! Determine how long until next observation
    external :: distribute_state_init_pdaf ! Distribute a state vector to model fields
    external :: prepoststep_ens_pdaf  ! User supplied pre/poststep routine


! ***************************
! ***   Initialize PDAF   ***
! ***************************

    if (mype_ens == 0) then
       write (*, '(/a,1x,a)') 'NEMO-PDAF', 'INITIALIZE PDAF - ONLINE MODE'
    end if

    ! **********************************************************
    ! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
    ! **********************************************************

    ! *** IO options ***
    screen = 2  ! Write screen output (1) for output, (2) add timings

    ! *** Filter specific variables
    filtertype = 7    ! Type of filter
      !   (1) SEIK
      !   (2) EnKF
      !   (3) LSEIK
      !   (4) ETKF
      !   (5) LETKF
      !   (6) ESTKF
      !   (7) LESTKF
    dim_ens = n_modeltasks  ! Size of ensemble for all ensemble filters
      !   We use n_modeltasks here, initialized in init_parallel_pdaf
    subtype = 0       ! subtype of filter:
      !   ESTKF:
      !     (0) Standard form of ESTKF
      !   LESTKF:
      !     (0) Standard form of LESTKF
    type_trans = 0    ! Type of ensemble transformation
      !   SEIK/LSEIK and ESTKF/LESTKF:
      !     (0) use deterministic omega
      !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
      !     (2) use product of (0) with random orthonormal matrix with
      !         eigenvector (1,...,1)^T
      !   ETKF/LETKF:
      !     (0) use deterministic symmetric transformation
      !     (2) use product of (0) with random orthonormal matrix with
      !         eigenvector (1,...,1)^T
    type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
      !   (0) fixed
      !   (1) global adaptive
      !   (2) local adaptive for LSEIK/LETKF/LESTKF
    forget = 1.0     ! Forgetting factor
    type_sqrt = 0     ! Type of transform matrix square-root
      !   (0) symmetric square root, (1) Cholesky decomposition
    incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)


    ! *********************************************************************
    ! ***   Settings for analysis steps  - used in call-back routines   ***
    ! *********************************************************************

    ! *** Forecast length (time interval between analysis steps) ***
    delt_obs = 496      ! Number of time steps between analysis/assimilation steps

    ! Which observations to assimilate
    assim_sst_cmems = .false.   ! Whether to assimilate SST data from CMEMS
    assim_ssh_mgrid = .false.   ! Whether to assimilate SSH data on model grid

    ! *** Localization settings
    locweight = 0     ! Type of localizating weighting
      !   (0) constant weight of 1
      !   (1) exponentially decreasing with SRADIUS
      !   (2) use 5th-order polynomial
      !   (3) regulated localization of R with mean error variance
      !   (4) regulated localization of R with single-point error variance

    ! Settings for CMEMS satellite SST
    rms_ssh_mgrid     = 0.8_8     ! Observation error stddev for SSH data on mopdel grid
    lradius_ssh_mgrid = 10000.0_8 ! Radius in km for lon/lat (or in grid points)
    sradius_ssh_mgrid = lradius_ssh_mgrid  ! Support radius for 5th-order polynomial
      ! or distance for 1/e for exponential weighting

    ! Settings for CMEMS satellite SST
    rms_obs_sst_cmems = 0.8_8     ! Observation error stddev for SST data from CMEMS
    mode_sst_cmems = 1            ! Observation mode for SST_CMEMS: 
                                  !  (0) linear interpolation to observation grid
                                  !  (1) super-obbing: average 4 observation values
    lradius_sst_cmems = 10000.0_8 ! Radius in km for lon/lat (or in grid points)
    sradius_sst_cmems = lradius_sst_cmems  ! Support radius for 5th-order polynomial
      ! or distance for 1/e for exponential weighting


! **********************************
! *** Namelist and screen output ***
! **********************************

    ! Read namelist file for PDAF
    call read_config_pdaf()

    ! Screen output for PDAF parameters
    if (mype_ens == 0) call init_info_pdaf()


! ************************************************
! *** Specify state vector and state dimension ***
! ************************************************

    ! Initialize dimension information for NEMO grid
    call set_nemo_grid()

    ! Setup state vector
    call setup_statevector(dim_state, dim_state_p)


! *****************************************************
! *** Set communicator within which PDAF operates.  ***
! *****************************************************

    call PDAF_set_comm_pdaf(COMM_ensemble)


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! *****************************************************

    whichinit: if (filtertype == 2) then
       ! *** EnKF with Monte Carlo init ***
       filter_param_i(1) = dim_state_p ! State dimension
       filter_param_i(2) = dim_ens     ! Size of ensemble
       filter_param_i(3) = rank_analysis_enkf ! Rank of speudo-inverse in analysis
       filter_param_i(4) = incremental ! Whether to perform incremental analysis
       filter_param_i(5) = 0           ! Smoother lag (not implemented here)
       filter_param_r(1) = forget      ! Forgetting factor

       call PDAF_init(filtertype, subtype, 0, &
            filter_param_i, 6, &
            filter_param_r, 2, &
            COMM_model, COMM_filter, COMM_couple, &
            task_id, n_modeltasks, filterpe, init_ens_pdaf, &
            screen, status_pdaf)
    else
       ! *** All other filters                       ***
       ! *** SEIK, LSEIK, ETKF, LETKF, ESTKF, LESTKF ***
       filter_param_i(1) = dim_state_p ! State dimension
       filter_param_i(2) = dim_ens     ! Size of ensemble
       filter_param_i(3) = 0           ! Smoother lag (not implemented here)
       filter_param_i(4) = incremental ! Whether to perform incremental analysis
       filter_param_i(5) = type_forget ! Type of forgetting factor
       filter_param_i(6) = type_trans  ! Type of ensemble transformation
       filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
       filter_param_r(1) = forget      ! Forgetting factor

       call PDAF_init(filtertype, subtype, 0, &
            filter_param_i, 7, &
            filter_param_r, 2, &
            COMM_model, COMM_filter, COMM_couple, &
            task_id, n_modeltasks, filterpe, init_ens_pdaf, &
            screen, status_pdaf)
    end if whichinit

    ! *** Check whether initialization of PDAF was successful ***
    if (status_pdaf /= 0) then
       write (*, '(/1x,a6,i3,a43,i4,a1/)') &
            'ERROR ', status_pdaf, &
            ' in initialization of PDAF - stopping! (PE ', mype_ens, ')'
       call abort_parallel()
    end if
      

! **************************************
! *** Initialise PDAF arrays for IAU ***
! **************************************

    call asm_inc_init_PDAF()


! **********************************
! *** Prepare ensemble forecasts ***
! **********************************

    call PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
         distribute_state_init_pdaf, prepoststep_ens_pdaf, status_pdaf)

  end subroutine init_pdaf

end module mod_init_pdaf
