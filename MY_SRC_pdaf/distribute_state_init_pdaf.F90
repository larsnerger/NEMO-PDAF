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
#if defined key_top
  use mod_statevector_pdaf, &
       only: sfields, id, n_trc, n_bgc1, n_bgc2, &
       jptra, jptra2, sv_bgc1, sv_bgc2, &
       id_dia, id_fla, id_cya, update_phys, update_phyto, update_nophyto
  use mod_nemo_pdaf, &
       only: ni_p, nj_p, nk_p, i0, j0, jp_tem, jp_sal, trb, &
             sshb, tsb, ub, vb, lbc_lnk, lbc_lnk_multi, &
             sshn, tsn, un, vn, trn, neuler, &
             xph, xpco2, xchl, xnetpp
#else
  use mod_statevector_pdaf, &
       only: sfields, id
  use mod_nemo_pdaf, &
       only: ni_p, nj_p, nk_p, i0, j0, jp_tem, jp_sal, &
             sshb, tsb, ub, vb, lbc_lnk, lbc_lnk_multi, &
             sshn, tsn, un, vn, neuler
#endif
  use mod_assimilation_pdaf, &
       only: ens_restart, assim_flag
  use mod_aux_pdaf, &
       only: state2field, transform_field_mv
  implicit none

! *** Arguments ***
  integer, intent(in) :: dim_p               !< PE-local state dimension
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector

! *** Local variables ***
  integer :: i, j, k, cnt ! Counters
#if defined key_top
  integer :: id_var       ! Counters
#endif
  logical :: firststep = .true.   ! Flag for first timestep
  integer :: verbose      ! Control verbosity


! ************************************
! Distribute state vector 2d variables
! ************************************

  ! Note: The loop limits account for the halo offsets i0 and j0
    ! if (mype==0) write (*,'(a,4x,a)') 'NEMO-PDAF', 'distribute state - DEACTIVATED'

  ! Aply field transformations
  if (mype==0) then
     verbose = 1
  else
     verbose = 0
  end if

  call transform_field_mv(2, state_p, 0, verbose) !21

  coldstart: if (.not. ens_restart) then

     if (mype==0) write (*,'(a,4x,a)') 'NEMO-PDAF', 'distribute state at initial time'

     ! SSH

     if (id%ssh > 0) then
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

     ! T
     if (id%temp > 0) then
        call state2field(state_p, &
                      tsn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_tem), &
                      sfields(id%temp)%off, sfields(id%temp)%ndims)
     end if

     ! S
     if (id%salt > 0) then
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
     if (id%uvel > 0) then
        call state2field(state_p, &
                      un(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                      sfields(id%uvel)%off, sfields(id%uvel)%ndims)
     end if

     ! V
     if (id%vvel > 0) then
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
       if (sv_bgc1(i)) then
         id_var=id%bgc1(i)
         call state2field(state_p, &
                       trn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, sfields(id_var)%jptrc), &
                       sfields(id_var)%off, sfields(id_var)%ndims)
         ! Fill halo regions
         call lbc_lnk('distribute_state_pdaf', trn(:, :, :, sfields(id_var)%jptrc), 'T', 1._pwp)
         trb(:, :, :, sfields(id_var)%jptrc) = trn(:, :, :, sfields(id_var)%jptrc)
       end if
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

     ! Set Euler step
     neuler = 0
     assim_flag = 1

  else coldstart

     ! Ensemble restart using restart fields from NEMO

     if (mype==0) write (*,'(a,4x,a)') 'NEMO-PDAF', 'Ensemble restart - distribute_state inactive'

     ! Set Euler step
     neuler = 0
     assim_flag = 1

  end if coldstart

end subroutine distribute_state_init_pdaf
