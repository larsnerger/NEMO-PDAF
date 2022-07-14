!$Id$
!BOP
!
! !MODULE:
module mod_assimilation_pdaf

! !DESCRIPTION:
! This module provides variables needed for the 
! assimilation within the routines of the dummy model.
! For simplicity, all assimilation-related variables
! are stored here, even if they are only used in
! the main program for the filter initialization.
! Most variables can be specified as a command line 
! argument.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  implicit none
  save
!EOP

! *** specific variables for NEMO-PDAF ***

  integer :: iday                ! Day of the month
  integer :: dim_state           ! Global model state dimension
  integer :: dim_state_p         ! Model state dimension for PE-local domain

  real                 :: domain_coords(2)  ! Grid point coordinates for local analysis domain
  integer(4), allocatable :: id_lstate_in_pstate(:) ! Indices of local state vector in global vector


! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF_offline                 ***

! !PUBLIC MEMBER FUNCTIONS:
! ! Settings for time stepping - available as command line options
  logical :: model_error   ! Control application of model error
  real    :: model_err_amp ! Amplitude for model error

! ! Settings for observations - available as command line options
  integer :: delt_obs      ! time step interval between assimilation steps

! ! General control of PDAF - available as command line options
  integer :: screen       ! Control verbosity of PDAF
                          ! (0) no outputs, (1) progess info, (2) add timings
                          ! (3) debugging output
  integer :: dim_ens      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
                          ! Number of EOFs to be used for SEEK
  integer :: filtertype   ! Select filter algorithm:
                          !   SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4)
                          !   LETKF (5), ESTKF (6), LESTKF (7), LEnKF (8), NETF (9), LNETF (10)
  integer :: subtype      ! Subtype of filter algorithm
                          !   SEEK: 
                          !     (0) evolve normalized modes
                          !     (1) evolve scaled modes with unit U
                          !     (2) fixed basis (V); variable U matrix
                          !     (3) fixed covar matrix (V,U kept static)
                          !   SEIK:
                          !     (0) ensemble forecast; new formulation
                          !     (1) ensemble forecast; old formulation
                          !     (2) fixed error space basis
                          !     (3) fixed state covariance matrix
                          !     (4) SEIK with ensemble transformation
                          !   EnKF:
                          !     (0) analysis for large observation dimension
                          !     (1) analysis for small observation dimension
                          !   LSEIK:
                          !     (0) ensemble forecast;
                          !     (2) fixed error space basis
                          !     (3) fixed state covariance matrix
                          !     (4) LSEIK with ensemble transformation
                          !   ETKF:
                          !     (0) ETKF using T-matrix like SEIK
                          !     (1) ETKF following Hunt et al. (2007)
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to SEIK subtypes 2/3
                          !   LETKF:
                          !     (0) LETKF using T-matrix like SEIK
                          !     (1) LETKF following Hunt et al. (2007)
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to LSEIK subtypes 2/3
                          !   ESTKF:
                          !     (0) standard ESTKF 
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to SEIK subtypes 2/3
                          !   LESTKF:
                          !     (0) standard LESTKF 
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to LSEIK subtypes 2/3
                          !   LEnKF:
                          !     (0) Standard form of EnKF with covariance localization
                          !   NETF:
                          !     (0) standard NETF 
                          !   LNETF:
                          !     (0) standard LNETF 
  integer :: incremental  ! Perform incremental updating in LSEIK
  integer :: dim_lag      ! Number of time instances for smoother

! ! Filter settings - available as command line options
!    ! General
  integer :: type_forget  ! Type of forgetting factor
  real    :: forget       ! Forgetting factor for filter analysis
                          !   (0) fixed
                          !   (1) global adaptive
                          !   (2) local adaptive for LSEIK/LETKF/LESTKF
  integer :: dim_bias     ! dimension of bias vector
!    ! SEEK
  integer :: int_rediag   ! Interval to perform re-diagonalization in SEEK
  real    :: epsilon      ! Epsilon for gradient approx. in SEEK forecast
!    ! ENKF
  integer :: rank_analysis_enkf  ! Rank to be considered for inversion of HPH
!    ! SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF
  integer :: type_trans    ! Type of ensemble transformation
                           ! SEIK/LSEIK:
                           ! (0) use deterministic omega
                           ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! ETKF/LETKF with subtype=4:
                           ! (0) use deterministic symmetric transformation
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! ESTKF/LESTKF:
                           ! (0) use deterministic omega
                           ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! NETF/LNETF:
                           ! (0) use random orthonormal transformation orthogonal to (1,...,1)^T
                           ! (1) use identity transformation
!    ! LSEIK/LETKF/LESTKF
  integer :: locweight     ! Type of localizing weighting of observations
                    ! For LESTKF, LETKF, and LSEIK
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
                    ! For LEnKF
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) 5th-order polynomial weight function
!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  integer :: type_sqrt     ! Type of the transform matrix square-root 
                    !   (0) symmetric square root, (1) Cholesky decomposition

! Constants for coordinate calculations
  real(8), parameter  :: pi     = 3.14159265358979323846_8
  real :: deg2rad = pi / 180.0_8


!    !  File output  
  integer(4)           :: timeDA 
  logical             :: shiftObsInWet
  real                :: flate
  logical             :: genEnsMeanYearly
  integer(4)          :: nyears
  integer             :: GaussTransf
  real                :: trafoConst
  logical             :: EnsDiagnos 
  logical             :: flateZ
  real                :: flateTOP, flateBOT
  integer             :: nLevFB,nLevFE 

!    ! ERGOM specific parameters
  integer(4)          :: ensfile_type           ! (1) for using ens file, (2) for using ergom output files



!    ! Other variables - _NOT_ available as command line options!
  integer :: covartype     ! For SEIK: Definition of ensemble covar matrix
                           ! (0): Factor (r+1)^-1 (or N^-1)
                           ! (1): Factor r^-1 (or (N-1)^-1) - real ensemble covar.
                           ! This setting is only for the model part; The definition
                           ! of P has also to be specified in PDAF_filter_init.
                           ! Only for upward-compatibility of PDAF!
  real    :: time          ! model time


!   ! Variables to handle multiple fields in the state vector
!   integer :: n_fields      !< number of fields in state vector
! 
!   ! Declare Fortran type holding the indices of model fields in the state vector
!   ! This can be extended to any number of fields - it serves to give each field a name
!   type field_ids
!      ! Ocean Physics
!      integer :: ssh
!      integer :: temp
!      integer :: salt
!      integer :: uvel
!      integer :: vvel
! 
!      ! ERGOM
!      integer :: oxy
!   end type field_ids
! 
!   ! Type variable holding field IDs in state vector
!   type(field_ids) :: id
! 
!   type state_field
!      integer :: ndims                  ! Number of field dimensions (2 or 3)
!      integer :: dim                    ! Dimension of the field
!      integer :: off                    ! Offset of field in state vector
!      character(len=10) :: variable     ! Name of field
!      character(len=20) :: name_incr    ! Name of field in increment file
!      character(len=20) :: name_rest_n  ! Name of field in restart file (n-field)
!      character(len=20) :: name_rest_b  ! Name of field in restart file (b-field)
!      character(len=30) :: file         ! File name stub to read field from
!      character(len=30) :: rst_file     ! Name of restart file
!      character(len=20) :: unit         ! Unit of variable
!      integer :: transform = 0          ! Type of variable transformation
!      real :: trafoConst = 0.0          ! Constant to shift value in transformation
!      integer :: limit = 0              ! Whether to limit the value of the variable
!                                        ! 0: no limits, 1: lower limit, 2: upper limit, 3: both limits
!      real :: max_limit = 0.0           ! Upper limit of variable
!      real :: min_limit = 0.0           ! Lower limit of variable
!   end type state_field
! 
!   type(state_field), allocatable :: sfields(:)

!$OMP threadprivate(domain_coords, id_lstate_in_pstate)

end module mod_assimilation_pdaf
