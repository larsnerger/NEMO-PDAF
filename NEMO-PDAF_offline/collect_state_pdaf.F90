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
       only: sfields, id, n_ocean, n_bgc
#ifndef key_PDAF_offline
  use mod_nemo_pdaf, &
       only: ni_p, nj_p, nk_p, i0, j0, &
       jp_tem, jp_sal, ndastp, &
       trb, sshb, tsb, ub, vb
#endif

  implicit none

! *** Arguments ***
  integer, intent(in) :: dim_p               !< PE-local state dimension
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector
#ifndef key_PDAF_offline
! *** Local variables ***
  integer :: i, j, k, cnt, i_var       ! Counters


  ! *********************************
  ! Collect state vector 2d variables
  ! *********************************


  ! Note: The loop limits account for the halo offsets i0 and j0

  ! SSH
  if (id%ssh > 0) then
     cnt = sfields(id%ssh)%off + 1
     do j = 1 + j0, nj_p + j0
        do i = 1 + i0, ni_p + i0
           state_p(cnt) = sshb(i, j)
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
              state_p(cnt) = tsb(i, j, k, jp_tem)
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
              state_p(cnt) = tsb(i, j, k, jp_sal)
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
              state_p(cnt) = ub(i, j, k)
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
              state_p(cnt) = vb(i, j, k)
              cnt = cnt + 1
           end do
        end do
     end do
  end if

  do i_var = n_ocean + 1, n_bgc
     cnt = sfields(i_var)%off + 1
     do k = 1, nk_p
       do j = 1 + j0, nj_p + j0
         do i = 1 + i0, ni_p + i0
           state_p(cnt) = trb(i, j, k, sfields(i_var)%jptrc)
           cnt = cnt + 1
         end do
       end do
     end do
  end do
#endif
end subroutine collect_state_pdaf
