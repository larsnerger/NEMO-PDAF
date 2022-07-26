module mod_nemo_pdaf
#ifndef key_PDAF_offline
  ! Include variables from NEMO
  use par_oce, &
       only: jpk, jpiglo, jpjglo
  use dom_oce, &
       only: nldi, nldj, nlei, nlej, glamt, gphit, &
       nimpp, njmpp, tmask, gdept_1d, ndastp
  use par_oce, &
       only: jp_tem, jp_sal
  use trc, &
       only: trb
  use oce, &
       only: sshb, tsb, ub, vb
  use lbclnk, &
       only: lbc_lnk, lbc_lnk_multi
#endif
  use mod_kind_pdaf
  implicit none

  ! *** NEMO model variables
#ifdef key_PDAF_offline
  integer                :: jpiglo, jpjglo, jpk        ! Global NEMO grid dimensions
  integer                :: nldi, nldj                 ! first inner index in i/j direction of sub-domain
  integer                :: nlei, nlej                 ! last inner index in i/j direction of sub-domain
  integer                :: nimpp, njmpp               ! start i,j of subdomain including halo
  logical                :: read_decomp=.false.        ! Whether to read domain decomposition from file
  integer                :: ndastp                     ! current date in integer format YYYYMMDD
  real(pwp), allocatable :: glamt(:,:)                 ! Longitudes
  real(pwp), allocatable :: gphit(:,:)                 ! Latitudes
  real(pwp), allocatable :: gdept_1d(:)                ! Depths
  real(pwp), allocatable :: tmask(:, :, :)             ! land maks at T-points
#endif

  ! *** Other grid variables
  real(pwp), allocatable :: tmp_4d(:,:,:,:)            ! 4D array used to represent full NEMO grid box

  integer                :: dim_2d                     ! Dimension of 2d grid box    
  integer                :: nwet                       ! Number of surface wet grid points
  integer                :: nwet3d                     ! Number of 3d wet grid points
  integer, allocatable   :: wet_pts(:,:)               ! Index array for wet grid points
                                                       !    (1) latitude, (2) langitude
                                                       !    (3) number wet layers, (4) index in 2d grid box
  integer, allocatable   :: idx_wet_2d(:,:)            ! Index array for wet_pts row index in 2d box
  integer, allocatable   :: idx_nwet(:,:)              ! Index array for wet_pts row index in wet surface grid points
  integer, allocatable   :: nlev_wet_2d(:,:)           ! Number of wet layers for ij position in 2d box

  integer                :: use_wet_state=0            ! 1: State vector contains full columns when surface grid point is wet
                                                       ! 2: State vector only contains wet grid points 
                                                       ! other: State vector contains 2d/3d grid box
  integer                :: i0, j0                     ! PE-local halo offsets
  integer                :: ni_p, nj_p, nk_p           ! Size of decomposed grid
  integer                :: istart, jstart             ! Start indices for internal local domain
  integer                :: dim_2d_p, dim_3d_p         ! Dimension of 2d/3d grid box of sub-domain   
  integer                :: sdim2d, sdim3d             ! 2D/3D dimension of field in state vector

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

    integer   :: i, j, k
    integer   :: cnt, cnt_all, cnt_layers
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
    dim_2d = jpiglo*jpjglo
    dim_2d_p = ni_p * nj_p
    dim_3d_p = ni_p * nj_p * nk_p

    ! Start indices for sub-domain without halo
    istart = nimpp+nldi-1
    jstart = njmpp+nldj-1

    ! Count number of wet surface points
    nwet = 0
    nwet3d = 0
    cnt = 0
    do k = 1, nk_p
       do j = 1, nj_p
          do i = 1, ni_p 
             cnt = cnt + 1
             if (tmask(i + i0, j + j0, k) == 1.0_pwp) then
                if (k==1) nwet = nwet + 1
                nwet3d = nwet3d + 1
             endif
          enddo
       enddo
    enddo

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

    !> Initialize index arrays 
    !! - for mapping from vector of wet points to 2d box
    !! - for mapping from vector model grid point coordinats to wet point index
    !! - for retrieving number of wet layers at a grid point coordinate
    !! these arrays also serve as mask arrays (0 indicates land point)
    allocate(idx_wet_2d(ni_p, nj_p))
    allocate(idx_nwet(ni_p, nj_p))
    allocate(nlev_wet_2d(ni_p, nj_p))
    ! call memcount(1, 'i', 3*nx*ny)

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

#ifdef key_PDAF_offline
  subroutine domain_decomposition(screen, mype, npes, file_decomp)
    use mod_parallel_pdaf, only: abort_parallel
    use mod_mppini_pdaf_offline, only: get_decomposition
    integer,      intent(in) :: screen
    integer,      intent(in) :: mype
    integer,      intent(in) :: npes
    character(*), intent(in) :: file_decomp
    
    integer                  :: n_domains_lon, n_domains_lat
    integer                  :: k
    integer, allocatable     :: dims_lat(:), dims_lon(:)
    integer, allocatable     :: nldi_all(:), nldj_all(:) ! first inner index in i/j direction for all PEs
    integer, allocatable     :: nlei_all(:), nlej_all(:) ! last inner index in i/j direction for all PEs
    ! ***********************************
    ! *** Define domain-decomposition ***
    ! ***********************************

    nldi = 1
    nldj = 1

    if (npes==1 .or. .not.read_decomp) then
      call get_decomposition(mype, jpiglo, jpjglo, npes, nlei, nlej, nimpp, njmpp)
    else
      ! Read decomposition file
      allocate(nldi_all(0:npes-1), nldj_all(0:npes-1))
      allocate(nlei_all(0:npes-1), nlej_all(0:npes-1))

      allocate(dims_lon(0:npes-1))
      allocate(dims_lat(0:npes-1))
      if (mype==0) write(*,'(/a,2x,a,a)') &
           'NEMO-PDAF','*** Read domain decomposition: ', file_decomp

      open(11,FILE=file_decomp)
      read(11,*) n_domains_lon, n_domains_lat
      if (mype==0) write (*,'(a,3x,a,2i5)') 'NEMO-PDAF', 'Number of tasks lon/lat: ', n_domains_lon, n_domains_lat

      do k = 0, n_domains_lon*n_domains_lat-1
         read(11, *)  nldj_all(k), nlej_all(k), nldi_all(k), nlei_all(k), dims_lon(k), dims_lat(k)
      end do
      close(11)
      nlei = nlei_all(mype) - nldi_all(mype) + 1
      nlej = nlej_all(mype) - nldj_all(mype) + 1
      nimpp = nldi_all(mype)
      njmpp = nldj_all(mype)
      if (n_domains_lon*n_domains_lat /= npes) then
         write (*,*) 'ERROR: number of domains in file inconsistent with number of processes'
         call abort_parallel()
      end if

      deallocate(dims_lon, dims_lat)
      deallocate(nldi_all, nldj_all)
      deallocate(nlei_all, nlej_all)
    end if

    ! Screen output
    if (screen > 0) then
       write (*,'(/a,3x,a)') 'NEMO-PDAF','Grid decomposition:' 
       write (*,'(a, 8x,a,2x,a,a,2x,a,a,1x,a,6(1x,a))') &
            'NEMO-PDAF','rank ', 'istart', '  iend', 'jstart', '  jend', '  idim', '  jdim'
       write (*,'(a,2x, a,i6,1x,2i7,2i7,2i7)') 'NEMO-PDAF', 'RANK',mype, nimpp+nldi-1, &
            nimpp+nlei-1, &
            njmpp+nldj-1, njmpp+nlej-1, nlei, nlej
    end if
  end subroutine domain_decomposition
#endif

  subroutine print_nemo_info(screen)
      use mod_parallel_pdaf, &
           only: mype => mype_ens, npes => npes_ens
      ! arguments
      integer, intent(in) :: screen
      if (mype==0 .and. screen>0) then
         write (*,'(/a,5x,a)') 'NEMO-PDAF', '*** NEMO: grid dimensions ***' 
         write(*,'(a,3x,2(6x,a),9x,a)') 'NEMO-PDAF', 'jpiglo','jpjglo','jpk' 
         write(*,'(a,3x,3i12)') 'NEMO-PDAF', jpiglo, jpjglo, jpk
         write(*,'(a,5x,a,i12)') 'NEMO-PDAF', 'Dimension of global 3D grid box', jpiglo*jpjglo*jpk
         write(*,'(a,5x,a,i12)') 'NEMO-PDAF', 'Number of global surface points', dim_2d
      end if

      if (npes>1 .and. screen>1) then
         write(*,'(a,2x,a,1x,i4,1x,a,i12)') &
              'NEMO-PDAF', 'PE', mype, 'Dimension of local 3D grid box', ni_p*nj_p*jpk
         write(*,'(a,2x,a,1x,i4,1x,a,i12)') &
              'NEMO-PDAF', 'PE', mype, 'Number of local surface points', dim_2d_p
      end if

      if (npes==1) then
         write(*,'(a,5x,a,3x,i11)') 'NEMO-PDAF', 'Number of wet surface points', nwet
         write(*,'(a,5x,a,8x,i11)') 'NEMO-PDAF', 'Number of 3D wet points', nwet3d
         write(*,'(a,5x,a,8x,i11)') 'NEMO-PDAF', '2D wet points * nlayers', nwet*jpk
      else 
         if (screen>1) then
            write(*,'(a,2x,a,1x,i4,2x,a,3x,i11)') &
                 'NEMO-PDAF', 'PE', mype, 'Number of wet surface points', nwet
            write(*,'(a,2x,a,1x,i4,2x,a,8x,i11)') &
                 'NEMO-PDAF', 'PE', mype, 'Number of 3D wet points', nwet3d
            write(*,'(a,2x,a,1x,i4,2x,a,8x,i11)') &
                 'NEMO-PDAF', 'PE', mype, '2D wet points * nlayers', nwet*jpk
         end if
      end if
  end subroutine print_nemo_info

  subroutine finalize_nemo_pdaf()
#ifdef key_PDAF_offline
    if (allocated(glamt)) deallocate(glamt)
    if (allocated(gphit)) deallocate(gphit)
    if (allocated(gdept_1d)) deallocate(gdept_1d)
    if (allocated(tmask)) deallocate(tmask)
#endif
    if (allocated(tmp_4d)) deallocate(tmp_4d)
    if (allocated(wet_pts)) deallocate(wet_pts)
    if (allocated(idx_wet_2d)) deallocate(idx_wet_2d)
    if (allocated(idx_nwet)) deallocate(idx_nwet)
    if (allocated(nlev_wet_2d)) deallocate(nlev_wet_2d)
  end subroutine finalize_nemo_pdaf
end module mod_nemo_pdaf
