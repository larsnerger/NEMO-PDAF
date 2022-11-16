!> Collecting the statevector variables
!!
!! The routine has to initialize the statevector of PDAF
!! from the fields of the model.
!! 
!! The routine is executed by each process that is
!! participating in the model integrations.
!! 
!! **Calling Sequence**
!! 
!!  - Called from:* `PDAFomi_assimilate_local`/`mod_assimilation_pdaf` (as U_coll_state)
!! 
subroutine collect_state_pdaf(dim_p, state_p)

  use mod_kind_pdaf
  use mod_parallel_pdaf, &
       only: mype=>mype_ens
  use mod_statevector_pdaf, &
       only: sfields, id
  use mod_nemo_pdaf, &
       only: ni_p, nj_p, nk_p, i0, j0, &
       jp_tem, jp_sal, ndastp
  use oce, &
       only: sshb, tsb, ub, vb

  implicit none

! *** Arguments ***
  integer, intent(in) :: dim_p               !< PE-local state dimension
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector

! *** Local variables ***
  integer :: i, j, k, cnt       ! Counters


  ! *********************************
  ! Collect state vector 2d variables
  ! *********************************


  ! Note: The loop limits account for the halo offsets i0 and j0

  ! SSH
  if (id%ssh > 0) then
     cnt = sfields(id%ssh)%off + 1
     do j = 1 + j0, nj_p + j0
        do i = 1 + i0, ni_p + i0
           state_p(cnt) = sshn(i, j)
           cnt = cnt + 1
        end do
     end do
  end if


  ! *********************************
  ! Collect state vector 3d variables
  ! *********************************

  ! T
  if (id%temp > 0) then
     cnt = sfields(id%temp)%off + 1
     do k = 1, nk_p
        do j = 1 + j0, nj_p + j0
           do i = 1 + i0, ni_p + i0
              state_p(cnt) = tsn(i, j, k, jp_tem)
              cnt = cnt + 1
           end do
        end do
     end do
  end if

  ! S
  if (id%salt > 0) then
     cnt = sfields(id%salt)%off + 1
     do k = 1, nk_p
        do j = 1 + j0, nj_p + j0
           do i = 1 + i0, ni_p + i0
              state_p(cnt) = tsn(i, j, k, jp_sal)
              cnt = cnt + 1
           end do
        end do
     end do
  end if

  ! U
  if (id%uvel > 0) then
     cnt = sfields(id%uvel)%off + 1
     do k = 1, nk_p
        do j = 1 + j0, nj_p + j0
           do i = 1 + i0, ni_p + i0
              state_p(cnt) = un(i, j, k)
              cnt = cnt + 1
           end do
        end do
     end do
  end if

  ! V
  if (id%vvel > 0) then
     cnt = sfields(id%vvel)%off + 1
     do k = 1, nk_p
        do j = 1 + j0, nj_p + j0
           do i = 1 + i0, ni_p + i0
              state_p(cnt) = vn(i, j, k)
              cnt = cnt + 1
           end do
        end do
     end do
  end if

end subroutine collect_state_pdaf
