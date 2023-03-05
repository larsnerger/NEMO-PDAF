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
subroutine distribute_state_pdaf(dim_p, state_p)

  use mod_kind_pdaf
  use mod_parallel_pdaf, &
       only: mype=>mype_ens
#if defined key_top
  use mod_iau_pdaf, &
       only: ssh_iau_pdaf, u_iau_pdaf, v_iau_pdaf, t_iau_pdaf, &
       s_iau_pdaf, bgc_iau_pdaf, div_damping_filter
  use mod_statevector_pdaf, &
       only: sfields, id, n_trc, jptra, jptra2, sv_bgc1, sv_bgc2, &
       id_dia, id_fla, id_cya, update_phys
  use mod_nemo_pdaf, &
       only: ni_p, nj_p, nk_p, i0, j0, jp_tem, jp_sal, trb, &
       sshb, tsb, ub, vb, tmask, lbc_lnk, lbc_lnk_multi, &
       trn, sshn, tsn, un, vn, &
       xph, xpco2, xchl, xnetpp
#else
  use mod_iau_pdaf, &
       only: ssh_iau_pdaf, u_iau_pdaf, v_iau_pdaf, t_iau_pdaf, &
       s_iau_pdaf, div_damping_filter
  use mod_statevector_pdaf, &
       only: sfields, id
  use mod_nemo_pdaf, &
       only: ni_p, nj_p, nk_p, i0, j0, jp_tem, jp_sal, trb, &
             sshb, tsb, ub, vb, tmask, lbc_lnk, lbc_lnk_multi, &
             sshn, tsn, un, vn
#endif
  use mod_aux_pdaf, &
       only: state2field, transform_field_mv
  implicit none

! *** Arguments ***
  integer, intent(in) :: dim_p               !< PE-local state dimension
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector

! *** Local variables ***
  logical :: dist_direct = .true. ! Flag for direct update of model fields
  integer :: i, j, k, cnt ! Counters
#if defined key_top
  integer :: id_var ! Counters
#endif
  integer :: verbose      ! Control verbosity


  ! **********************************************
  ! Only distribute full state on first time step.
  ! Otherwise compute increment.
  ! **********************************************

  ! Aply field transformations
  if (mype==0) then
     verbose = 1
  else
     verbose = 0
  end if

  call transform_field_mv(2, state_p, 21, verbose)  !21

!  direct: if (dist_direct) then
     ! ************************************
     ! Distribute state vector 2d variables
     ! ************************************

     if (mype==0) write (*,'(a,4x,a)') 'NEMO-PDAF', 'distribute state to model fields'

     ! SSH

     if (id%ssh > 0) then
        if (update_phys) &
             call state2field(state_p, sshn(1+i0:ni_p+i0, 1+j0:nj_p+j0), &
             sfields(id%ssh)%off, sfields(id%ssh)%ndims)

        ! Fill halo regions
        call lbc_lnk('distribute_state_pdaf', sshn, 'T', 1.)

        ! Update before field
        sshb = sshn
     endif


     ! ************************************
     ! Distribute state vector 3d variables
     ! ************************************

     if ((id%temp>0 .or. id%salt>0 .or. id%uvel>0 .or. id%vvel>0) .and. update_phys) &
          write (*,'(a,4x,a)') 'NEMO-PDAF', 'distribute_state: update physics'

     ! T
     if (id%temp > 0 .and. update_phys) then
        call state2field(state_p, &
                      tsn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_tem), &
                      sfields(id%temp)%off, sfields(id%temp)%ndims)
     end if

     ! S
     if (id%salt > 0 .and. update_phys) then
        call state2field(state_p, &
                tsn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_sal), &
                sfields(id%salt)%off, sfields(id%salt)%ndims)
     end if

     if (id%temp>0 .or. id%salt>0) then
        ! Fill halo regions
        call lbc_lnk_multi('distribute_state_pdaf', tsn(:, :, :, jp_tem), 'T', &
             1., tsn(:, :, :, jp_sal), 'T', 1.)

        ! Update before fields
        tsb(:,:,:,jp_tem) = tsn(:,:,:,jp_tem)
        tsb(:,:,:,jp_sal) = tsn(:,:,:,jp_sal)
     end if

     ! U
     if (id%uvel > 0 .and. update_phys) then
        call state2field(state_p, &
                      un(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                      sfields(id%uvel)%off, sfields(id%uvel)%ndims)
     end if

     ! V
     if (id%vvel > 0 .and. update_phys) then
        call state2field(state_p, &
                      vn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                      sfields(id%vvel)%off, sfields(id%vvel)%ndims)
     end if

     if (id%uvel>0 .or. id%vvel>0) then
        ! Fill halo regions
        call lbc_lnk_multi('distribute_state_pdaf', un, 'U', -1., vn, 'V', -1.)

        ! Update before fields
        ub = un
        vb = vn
     end if

     ! BGC
#if defined key_top
     do i = 1, jptra
        updatebgc1: if (sv_bgc1(i) .and. sfields(id%bgc1(i))%update) then

           if (mype==0) write (*,'(a,1x,a,1x,a)') 'NEMO-PDAF', &
                'distribute_state, update ERGOM variable ', sfields(id%bgc1(i))%variable

           ! Update 3 phytoplankton variables
           id_var=id%bgc1(i)
           call state2field(state_p, &
                trn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, sfields(id_var)%jptrc), &
                sfields(id_var)%off, sfields(id_var)%ndims)

           ! Fill halo regions
           call lbc_lnk('distribute_state_pdaf', trn(:, :, :, sfields(id_var)%jptrc), 'T', &
                1._pwp)
           trb(:, :, :, sfields(id_var)%jptrc) = trn(:, :, :, sfields(id_var)%jptrc)

        end if updatebgc1
     end do

    do i = 1, jptra2
       if (sv_bgc2(i)) then
          id_var=id%bgc2(i)
          select case (i)
          case (1)
             call state2field(state_p, &
                  xpco2(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                  sfields(id_var)%off, sfields(id_var)%ndims)
             ! Fill halo regions
             call lbc_lnk('distribute_state_pdaf', xpco2, 'T', 1._pwp)
          case (2)
             call state2field(state_p, &
                  xph(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                  sfields(id_var)%off, sfields(id_var)%ndims)
             ! Fill halo regions
             call lbc_lnk('distribute_state_pdaf', xph, 'T', 1._pwp)
          case (3)
             call state2field(state_p, &
                  xchl(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                  sfields(id_var)%off, sfields(id_var)%ndims)
             ! Fill halo regions
             call lbc_lnk('distribute_state_pdaf', xchl, 'T', 1._pwp)
          case (4)
             call state2field(state_p, &
                  xnetpp(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                  sfields(id_var)%off, sfields(id_var)%ndims)
             ! Fill halo regions
             call lbc_lnk('distribute_state_pdaf', xnetpp, 'T', 1._pwp)
          end select
       end if
    end do
#endif

!  else direct
!     ! ***********************************************
!     ! Compute increment for state vector 2d variables
!     ! ***********************************************
!
!     if (mype==0) write (*,'(a,4x,a)') 'NEMO-PDAF', 'distribute state increment'
!
!     ! SSH
!     if (id%ssh > 0) then
!        cnt = sfields(id%ssh)%off + 1
!        do j = 1 + j0, nj_p + j0
!           do i = 1 + i0, ni_p + i0
!              ssh_iau_pdaf(i, j) = state_p(cnt) - sshb(i, j)
!              cnt = cnt + 1
!           end do
!        end do
!
!        ! Fill halo regions
!        call lbc_lnk('distribute_state_pdaf', ssh_iau_pdaf, 'T', 1.)
!     end if
!
!
!     ! ***********************************************
!     ! Compute increment for state vector 3d variables
!     ! ***********************************************
!     ! T
!     if (id%temp > 0) then
!        cnt = sfields(id%temp)%off + 1
!        do k = 1, nk_p
!           do j = 1 + j0, nj_p + j0
!              do i = 1 + i0, ni_p + i0
!                 t_iau_pdaf(i, j, k) = state_p(cnt) - tsb(i, j, k, jp_tem)
!                 cnt = cnt + 1
!              end do
!           end do
!        end do
!     end if
!
!     ! S
!     if (id%salt > 0) then
!        cnt = sfields(id%salt)%off + 1
!        do k = 1, nk_p
!           do j = 1 + j0, nj_p + j0
!              do i = 1 + i0, ni_p + i0
!                 s_iau_pdaf(i, j, k) = state_p(cnt) - tsb(i, j, k, jp_sal)
!                 cnt = cnt + 1
!              end do
!           end do
!        end do
!     end if
!
!     if (id%temp>0 .or. id%salt>0) then
!        ! Fill halo regions
!        call lbc_lnk_multi('distribute_state_pdaf', t_iau_pdaf(:, :, :), 'T', &
!             1., s_iau_pdaf(:, :, :), 'T', 1.)
!     end if
!
!     ! U
!     if (id%uvel > 0) then
!        cnt = sfields(id%uvel)%off + 1
!        do k = 1, nk_p
!           do j = 1 + j0, nj_p + j0
!              do i = 1 + i0, ni_p + i0
!                 u_iau_pdaf(i, j, k) = state_p(cnt) - ub(i, j, k)
!                 cnt = cnt + 1
!              end do
!           end do
!        end do
!     end if
!
!     ! V
!     if (id%vvel > 0) then
!        cnt = sfields(id%vvel)%off + 1
!        do k = 1, nk_p
!           do j = 1 + j0, nj_p + j0
!              do i = 1 + i0, ni_p + i0
!                 v_iau_pdaf(i, j, k) = state_p(cnt) - vb(i, j, k)
!                 cnt = cnt + 1
!              end do
!           end do
!        end do
!     end if
!
!     if (id%uvel>0 .or. id%vvel>0) then
!        ! Fill halo regions
!        call lbc_lnk_multi('distribute_state_pdaf', u_iau_pdaf, 'U', -1., &
!             v_iau_pdaf, 'V', -1.)
!     end if
!     call div_damping_filter
!
!     do id_var = 1, n_trc
!        cnt = sfields(id%trcs(id_var))%off + 1
!        do k = 1, nk_p
!           do j = 1 + j0, nj_p + j0
!              do i = 1 + i0, ni_p + i0
!                 bgc_iau_pdaf(i, j, k, id_var) = state_p(cnt) - trb(i, j, k, sfields(id_var)%jptrc)
!                 cnt = cnt + 1
!              end do
!           end do
!        end do
!        ! Fill halo regions
!        call lbc_lnk('distribute_state_pdaf', bgc_iau_pdaf(:, :, :, id_var), 'T', 1.)
!     end do
!  end if direct

end subroutine distribute_state_pdaf
