!> Set dimension of local model state
!!
!! The routine is called during analysis step
!! in `PDAF_X_update` in the loop over all local
!! analysis domains. It has to set the dimension
!! of the local model state on the current analysis
!! domain.
!!
!! This code is for NEMO-PDAF
!! 
!! - Called from: `PDAFomi_assimilate_local`/`mod_assimilation_pdaf`
!!
subroutine init_dim_l_pdaf(step, domain_p, dim_l)

  use mod_kind_pdaf
  use mod_assimilation_pdaf, &
       only: domain_coords, id_lstate_in_pstate, deg2rad, dim_state_p
  use mod_statevector_pdaf, &
       only: n_fields, sfields
  use mod_nemo_pdaf, &
       only: glamt, gphit, i0, j0, &
             wet_pts, use_wet_state, nwet, sdim2d

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step     !< Current time step
  integer, intent(in)  :: domain_p !< Current local analysis domain
  integer, intent(out) :: dim_l    !< Local state dimension

! *** Local variables
  integer :: i, cnt, ifield        ! Counters
  integer :: id_surf               ! state vector index of surface grid point
  integer :: id_i, id_j            ! Grid coordinates for local analysis domain
  integer :: domain_all

  ! ****************************************
  ! *** Initialize local state dimension ***
  ! ****************************************

  ! *************************************************************
  ! dimension = (number of 2D state variables) 
  !     + (number of 3D variables * number of ocean vertical points).
  !
  ! The number of ocean vertical points is stored in wet_pts(3,:)
  ! (we do not include land points in our local state vector).
  ! *************************************************************

  ! Determine local state dimension - distinguish 2D and 3D fields
  dim_l = 0
  do i = 1, n_fields
     if (sfields(i)%ndims==2) then
        dim_l = dim_l + 1
     else
        dim_l = dim_l + wet_pts(3, domain_p)
     end if
  end do


  ! **********************************************
  ! *** Initialize coordinates of local domain ***
  ! **********************************************
  domain_all = wet_pts(4, domain_p)

  ! get grid point indices
  id_i = wet_pts(6, domain_p)
  id_j = wet_pts(7, domain_p)

  ! Use T-values to get local coordinates
  ! the coordinates are stored in radians (as required by PDAFOMI)
  domain_coords(1) = glamt(id_i + i0, id_j + j0) * deg2rad
  domain_coords(2) = gphit(id_i + i0, id_j + j0) * deg2rad

! ******************************************************
! *** Initialize array of indices of the local state ***
! ***  vector elements in the global state vector.   ***
! ******************************************************

  ! **********************************************************
  ! A local domain consists of all ocean points in a vertical
  ! column. Such a domain will have indices (x,y,:).
  ! Each 2d state variable in the global statevector will be
  ! located at ( (y-1)*dim_longitude ) + x + offset(field).
  ! Each 3d state variables in the global statevector will be
  ! located at ( (z-1)*dim_longitude*dim_latitude ) +
  ! ( (y-1)*dim_longitude ) + x + offset(field), where z can vary
  ! over all *ocean* points in the vertical column.
  ! offset(field) is the variable-specific offset.
  ! **********************************************************

  ! Allocate array
  if (allocated(id_lstate_in_pstate)) deallocate(id_lstate_in_pstate)
  allocate(id_lstate_in_pstate(dim_l))

  cnt = 1
  do ifield = 1, n_fields
     ! Set indices
     if (use_wet_state==1) then
        ! state vector contains full columns of surface wet points
        id_surf = domain_p + sfields(ifield)%off

        id_lstate_in_pstate(cnt) = id_surf
        cnt = cnt + 1
        if (sfields(ifield)%ndims==3) then
           do i = 2, wet_pts(3,domain_p)
              id_lstate_in_pstate(cnt) = id_surf + (i-1)*nwet
              cnt = cnt + 1
           enddo
        end if
     else if (use_wet_state == 2) then
        ! state vector only contains wet points - stored in leading vertical order
        if (sfields(ifield)%ndims==3) then
           id_surf = wet_pts(5, domain_p) + sfields(ifield)%off
           do i = 1, wet_pts(3,domain_p)
              id_lstate_in_pstate(cnt) = id_surf + i - 1
              cnt = cnt  + 1
           enddo

        else

           ! 2D field
           id_lstate_in_pstate(cnt) = domain_p + sfields(ifield)%off
           cnt = cnt  + 1

        end if
     else
        ! state vector contains all grid points
        id_surf = domain_all + sfields(ifield)%off
        id_lstate_in_pstate(cnt) = id_surf
        cnt = cnt + 1

        if (sfields(ifield)%ndims==3) then

           do i = 2, wet_pts(3,domain_p)
              id_lstate_in_pstate(cnt) = id_surf + (i-1)*sdim2d
              cnt = cnt + 1
           enddo

        end if
     end if
  end do

  if (dim_l > n_fields) then
     if (id_lstate_in_pstate(wet_pts(3,domain_p)) > dim_state_p) then
        write(*,*) 'Error: please check the global indices for local state vector'
     endif
  end if
end subroutine init_dim_l_pdaf
