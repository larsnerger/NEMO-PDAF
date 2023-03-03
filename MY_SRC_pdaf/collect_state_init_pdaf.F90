!> Collecting the statevector variables for ensemble initialization
!!
!! The routine has to initialize the statevector of PDAF
!! from the fields of the model. This variant is used
!! to initialize the ensemble central state. In contrast 
!! to collect_state_pdaf it does no contain a call to 
!! transform_field_mv.
!!
!! The routine is executed by each process that is
!! participating in the model integrations.
!!
!! **Calling Sequence**
!!
!!  - Called from:* `init_ens_pdaf`
!!
subroutine collect_state_init_pdaf(dim_p, state_p)

  use mod_kind_pdaf
#if defined key_top
  use mod_statevector_pdaf, &
       only: sfields, id, n_trc, n_bgc1, n_bgc2, &
       jptra, jptra2, sv_bgc1, sv_bgc2
  use mod_nemo_pdaf, &
       only: ni_p, nj_p, nk_p, i0, j0, &
       jp_tem, jp_sal, ndastp, &
       trb, sshb, tsb, ub, vb, &
       xph, xpco2, xchl, xnetpp
#else
  use mod_statevector_pdaf, &
       only: sfields, id
  use mod_nemo_pdaf, &
       only: ni_p, nj_p, nk_p, i0, j0, &
       jp_tem, jp_sal, ndastp, &
       sshb, tsb, ub, vb
#endif
  use mod_aux_pdaf, &
       only: field2state


  implicit none

! *** Arguments ***
  integer, intent(in) :: dim_p               !< PE-local state dimension
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector

! *** Local variables ***
  real(pwp) :: missing_value = 1.0e20
  integer :: i       ! Counters


  ! *********************************
  ! Collect state vector 2d variables
  ! *********************************


  ! Note: The loop limits account for the halo offsets i0 and j0

  ! SSH
  if (id%ssh > 0) then
     call field2state(sshb(1+i0:ni_p+i0, 1+j0:nj_p+j0), &
                      state_p, &
                      sfields(id%ssh)%off, sfields(id%ssh)%ndims, missing_value)
  end if


  ! *********************************
  ! Collect state vector 3d variables
  ! *********************************

  ! T
  if (id%temp > 0) then
     call field2state(tsb(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_tem), &
                      state_p, &
                      sfields(id%temp)%off, sfields(id%temp)%ndims, missing_value)
  end if

  ! S
  if (id%salt > 0) then
     call field2state(tsb(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_sal), &
                      state_p, &
                      sfields(id%salt)%off, sfields(id%salt)%ndims, missing_value)
  end if

  ! U
  if (id%uvel > 0) then
       call field2state(ub(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                      state_p, &
                      sfields(id%uvel)%off, sfields(id%uvel)%ndims, missing_value)
  end if

  ! V
  if (id%vvel > 0) then
      call field2state(vb(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                   state_p, &
                   sfields(id%vvel)%off, sfields(id%vvel)%ndims, missing_value)
  end if

#if defined key_top
  ! BGC
  do i = 1, jptra
    if (sv_bgc1(i)) then
      call field2state(trb(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, sfields(id%bgc1(i))%jptrc), &
             state_p, &
             sfields(id%bgc1(i))%off, sfields(id%bgc1(i))%ndims, missing_value)
    end if
  end do

  do i = 1, jptra2
    if (sv_bgc2(i)) then
      select case (i)
      case (1)
         call field2state(xpco2(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
              state_p, &
              sfields(id%bgc2(i))%off, sfields(id%bgc2(i))%ndims, missing_value)
      case (2)
         call field2state(xph(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
              state_p, &
              sfields(id%bgc2(i))%off, sfields(id%bgc2(i))%ndims, missing_value)
      case (3)
         call field2state(xchl(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
              state_p, &
              sfields(id%bgc2(i))%off, sfields(id%bgc2(i))%ndims, missing_value)
      case (4)
         call field2state(xnetpp(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
              state_p, &
              sfields(id%bgc2(i))%off, sfields(id%bgc2(i))%ndims, missing_value)
      end select
    end if
  end do
#endif

end subroutine collect_state_init_pdaf
