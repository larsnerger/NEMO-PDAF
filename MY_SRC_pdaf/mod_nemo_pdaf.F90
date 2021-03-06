module mod_nemo_pdaf

  use mod_kind_pdaf

  ! Include variables from NEMO
  use par_oce, &
       only: jpk, jpiglo, jpjglo
  use dom_oce, &
       only: nldi, nldj, nlei, nlej, glamt, gphit, &
       nimpp, njmpp, tmask, gdept_1d, ndastp
  use par_oce, &
       only: jp_tem, jp_sal

  ! *** NEMO model variables
!  integer :: jpiglo, jpjglo, jpk        ! Global NEMO grid dimensions
!  integer :: nldi, nldj                 ! first inner index in i/j direction of sub-domain
!  integer :: nlei, nlej                 ! last inner index in i/j direction of sub-domain
!  integer :: nimpp, njmpp               ! start i,j of subdomain including halo

!  real, allocatable :: glamt(:,:)       ! Longitudes
!  real, allocatable :: gphit(:,:)       ! Latitudes
!  real, allocatable :: gdept_1d(:)      ! Depths
  
  ! *** Other grid variables
  real(pwp), allocatable :: tmp_4d(:,:,:,:)     ! 4D array used to represent full NEMO grid box
  real(pwp), allocatable :: lat1(:), lon1(:)    ! Vectors holding latitude and latitude

  integer :: dim_2d                        ! Dimension of 2d grid box    
  integer :: nwet                          ! Number of surface wet grid points
  integer :: nwet3d                        ! Number of 3d wet grid points
  integer, allocatable :: wet_pts(:,:)     ! Index array for wet grid points
                           ! (1) latitude, (2) langitude, (3) number wet layers, (4) index in 2d grid box
  integer, allocatable :: idx_wet_2d(:,:)  ! Index array for wet_pts row index in 2d box
  integer, allocatable :: idx_nwet(:,:)    ! Index array for wet_pts row index in wet surface grid points
  integer, allocatable :: nlev_wet_2d(:,:) ! Number of wet layers for ij position in 2d box

  integer :: use_wet_state=0               ! 1: State vector contains full columns when surface grid point is wet
                                           ! 2: State vector only contains wet grid points 
                                           ! other: State vector contains 2d/3d grid box

  integer :: i0, j0                         ! PE-local halo offsets
  integer :: ni_p, nj_p, nk_p               ! Size of decomposed grid
  integer :: istart, jstart                 ! Start indices for internal local domain
  integer :: dim_2d_p, dim_3d_p             ! Dimension of 2d/3d grid box of sub-domain   
  real(pwp), allocatable :: lat1_p(:), lon1_p(:) ! Vectors holding latitude and latitude for decomposition
  integer :: sdim2d, sdim3d                 ! 2D/3D dimension of field in state vector

  ! *** File name and path to read grid information
  character(len=200)  :: path_dims         ! Path for NEMO file holding dimensions
  character(len=80)   :: file_dims         ! File name NEMO file holding dimensions

contains

  !> Initialize grid information for the DA
  !!
  !! This routine initializes grid information for the data assimilation
  !! In particular index information is initialize to map in between
  !! the model grid and the state vector
  !!
  subroutine set_nemo_grid()

    use mod_kind_pdaf 
    use mod_parallel_pdaf, &
         only: mype_model
    use PDAFomi, &
         only: PDAFomi_set_domain_limits
    use mod_assimilation_pdaf, &
         only: deg2rad

    implicit none

    integer :: i, j, k
    integer :: cnt, cnt_all, cnt_layers
    real(pwp) :: lim_coords(2,2)      ! Limiting coordinates of sub-domain
    
! *** set dimension of 2d and 3d fields in state vector ***

    ! Local dimensions
    ni_p = nlei - nldi + 1
    nj_p = nlej - nldj + 1
    nk_p = jpk

    ! Compute halo offset
    i0 = nldi - 1
    j0 = nldj - 1

    ! Size of 2d/3d boxes without halo
    dim_2d_p = ni_p * nj_p
    dim_3d_p = ni_p * nj_p * nk_p

    ! Start indices for sub-domain without halo
    istart = nimpp+nldi-1
    jstart = njmpp+nldj-1

    ! Set coordinate vectors
    allocate(lat1_p(nj_p), lon1_p(ni_p))
    lat1_p = gphit(1, j0 + 1 : j0 + nj_p)
    lon1_p = glamt(i0 + 1 : i0 + ni_p, 1)

    ! Count number of wet surface points
    nwet = 0
    do j = 1, nj_p
       do i = 1, ni_p
          if (tmask(i + i0, j + j0, 1) == 1.0_pwp) then
             nwet = nwet + 1
          end if
       end do
    end do

    ! Initialize index arrays
    ! - for mapping from nx*ny grid to vector of wet points
    ! - mask for wet points

    if (nwet > 0) then

       allocate(wet_pts(7, nwet))

       cnt = 0
       cnt_all = 0
       do j = 1, nj_p
          do i = 1, ni_p
             cnt_all = cnt_all + 1
             if (tmask(i + i0, j + j0, 1) == 1.0_pwp) then
                cnt = cnt + 1

                wet_pts(1,cnt) = i + nimpp + nldi - 2     ! Global longitude index 
                wet_pts(2,cnt) = j + njmpp + nldj - 2     ! Global latitude index
                wet_pts(6,cnt) = i                ! Longitude index in subdomain
                wet_pts(7,cnt) = j                ! Latitidue index in subdomain

                ! Determine number of wet layers
                cnt_layers = 0
                do k = 1, nk_p
                   if (tmask(i + i0, j + j0, k) == 1.0_pwp) cnt_layers = cnt_layers + 1
                end do
                wet_pts(3,cnt) = cnt_layers
                wet_pts(4,cnt) = cnt_all
             end if
          end do
       end do

       ! row 5 stores wet_pts index for vertical column storage in state vector
       wet_pts(5,1) = 1
       do i = 2 , nwet
          wet_pts(5,i) = wet_pts(5,i-1) + wet_pts(3,i-1)
       end do

    else

       write (*, '(8x,a,i3)') 'WARNING: No valid local domains, PE=', mype_model
       nwet = -1

       allocate(wet_pts(3, 1))
       wet_pts(:, 1) = 1
       wet_pts(3, 1) = 0    ! Set number of layers=0
    end if

  ! Initialize index arrays 
  ! - for mapping from vector of wet points to 2d box
  ! - for mapping from vector model grid point coordinats to wet point index
  ! - for retrieving number of wet layers at a grid point coordinate
  ! these arrays also serve as mask arrays (0 indicates land point)

  allocate(idx_wet_2d(ni_p, nj_p))
  allocate(idx_nwet(ni_p, nj_p))
  allocate(nlev_wet_2d(ni_p, nj_p))

  idx_wet_2d = 0
  idx_nwet = 0
  nlev_wet_2d = 0
  do i = 1 , nwet
     idx_wet_2d(wet_pts(6,i), wet_pts(7,i)) = wet_pts(4,i)
     nlev_wet_2d(wet_pts(6,i), wet_pts(7,i)) = wet_pts(3,i)
  end do
  if (use_wet_state/=2) then
     do i = 1 , nwet
        idx_nwet(wet_pts(6,i), wet_pts(7,i)) = i
     end do
  else
     do i = 1 , nwet
        idx_nwet(wet_pts(6,i), wet_pts(7,i)) = wet_pts(5,i)
     end do
  end if


! ********************************************************
! *** Set dimension of 2D and 3D field in state vector ***
! ********************************************************
  
    if (use_wet_state==1) then
       ! State vector contains full columns when surface grid point is wet
       sdim3d = abs(nwet)*jpk
       sdim2d = abs(nwet)
    elseif (use_wet_state==2) then
       ! State vector only contains wet grid points
       sdim3d = abs(nwet3d)
       sdim2d = abs(nwet)
    else
       ! State vector contains 2d/3d grid box
       sdim3d = dim_3d_p
       sdim2d = dim_2d_p
    end if


! ******************************************************************
! *** Specify domain limits to limit observations to sub-domains ***
! ******************************************************************

    lim_coords(1,1) = glamt(i0 + 1, j0 + 1) * deg2rad
    lim_coords(1,2) = glamt(i0 + ni_p, j0 + 1) * deg2rad
    lim_coords(2,1) = gphit(i0 + ni_p, j0 + nj_p) * deg2rad
    lim_coords(2,2) = gphit(i0 + 1, j0 + 1) * deg2rad

    call PDAFomi_set_domain_limits(lim_coords)

  end subroutine set_nemo_grid

end module mod_nemo_pdaf
