!> Distributing the statevector variables, computing the
!> statevector increments
!!
!! This routine either initializes the full fields of
!! the model from the statevector of PDAF (first timestep),
!! or computes the statevector increments (all other timesteps).
!! For all other timesteps, the increments are added to the
!! model during the NEMO timestepping routine. See `mod_iau_pdaf` for details.
!! 
!! The routine is executed by each process that is
!! participating in the model integrations.
!! 
!! **Calling Sequence**
!! 
!!  - Called from: `PDAF_get_state` (as U_dist_state)
!! 
!!  - Called from: `PDAFomi_assimilate_local` (as U_dist_state)
!!
subroutine distribute_state_init_pdaf(dim_p, state_p)

  use mod_kind_pdaf
  use mod_parallel_pdaf, &
       only: mype=>mype_ens
  use mod_iau_pdaf, &
       only: ssh_iau_pdaf, u_iau_pdaf, v_iau_pdaf, t_iau_pdaf, &
       s_iau_pdaf
  use mod_statevector_pdaf, &
       only: sfields, id
  use mod_nemo_pdaf, &
       only: ni_p, nj_p, nk_p, i0, j0, jp_tem, jp_sal, &
       tmask
  use oce, &
       only: sshn, tsn, un, vn, sshb, tsb, ub, vb
  use lbclnk, &
       only: lbc_lnk, lbc_lnk_multi
use dom_oce, only: neuler
USE in_out_manager, only: nit000
  implicit none

! *** Arguments ***
  integer, intent(in) :: dim_p               !< PE-local state dimension
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector

! *** Local variables ***
  integer :: i, j, k, cnt       ! Counters
  logical :: firststep = .true. ! Flag for first timestep



! ************************************
! Distribute state vector 2d variables
! ************************************

  ! Note: The loop limits account for the halo offsets i0 and j0

  if (mype==0) write (*,'(a,4x,a)') 'NEMO-PDAF', 'distribute state at initial time'

  ! SSH

  if (id%ssh > 0) then
     cnt = sfields(id%ssh)%off + 1
     do j = 1 + j0, nj_p + j0
        do i = 1 + i0, ni_p + i0
           if (tmask(i, j, 1) == 1.0_pwp) then
              sshn(i, j) = state_p(cnt)
           end if
           cnt = cnt + 1
        end do
     end do

     ! Fill halo regions
     call lbc_lnk('distribute_state_pdaf', sshn, 'T', 1.)

     ! Update before field 
     sshb = sshn
  endif


! ************************************
! Distribute state vector 3d variables
! ************************************

  ! T
  if (id%temp > 0) then
     cnt = sfields(id%temp)%off + 1
     do k = 1, nk_p
        do j = 1 + j0, nj_p + j0
           do i = 1 + i0, ni_p + i0
              if (tmask(i, j, 1) == 1.0_pwp) then
                 tsn(i, j, k, jp_tem) = state_p(cnt)
              end if
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
              if (tmask(i, j, 1) == 1.0_pwp) then
                 tsn(i, j, k, jp_sal) = state_p(cnt)
              end if
              cnt = cnt + 1
           end do
        end do
     end do
  end if

  if (id%temp>0 .or. id%salt>0) then
     ! Fill halo regions
     call lbc_lnk_multi('distribute_state_pdaf', tsn(:, :, :, jp_tem), 'T', &
          1., tsn(:, :, :, jp_sal), 'T', 1.)

     ! Update before fields
     tsb(:,:,:,jp_tem) = tsb(:,:,:,jp_tem)
     tsb(:,:,:,jp_sal) = tsb(:,:,:,jp_sal)
  end if

  ! U
  if (id%uvel > 0) then
     cnt = sfields(id%uvel)%off + 1
     do k = 1, nk_p
        do j = 1 + j0, nj_p + j0
           do i = 1 + i0, ni_p + i0
              if (tmask(i, j, 1) == 1.0_pwp) then
                 un(i, j, k) = state_p(cnt)
              end if
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
              if (tmask(i, j, 1) == 1.0_pwp) then
                 vn(i, j, k) = state_p(cnt)
              end if
              cnt = cnt + 1
           end do
        end do
     end do
  end if

  if (id%uvel>0 .or. id%vvel>0) then
     ! Fill halo regions
     call lbc_lnk_multi('distribute_state_pdaf', un, 'U', -1., vn, 'V', -1.)     

     ! Update before fields
     ub = un
     vb = vn
  end if

  ! Set Euler step
  neuler = 0

end subroutine distribute_state_init_pdaf
