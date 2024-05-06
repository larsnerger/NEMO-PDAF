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
  use mod_statevector_pdaf, &
       only: sfields, id
  use mod_nemo_pdaf, &
       only: ni_p, nj_p, nk_p, i0, j0, &
       jp_tem, jp_sal, sshn, tsn, un, vn
#if defined key_top
  use mod_statevector_pdaf, &
       only: n_trc, n_bgc_diag, n_bgc_diag, &
       jpbgc_prog, jpbgc_diag, sv_bgc_prog, sv_bgc_diag
  use mod_nemo_pdaf, &
       only: trb, trn, sshn, tsn, un, vn, &
       xph, xpco2, xchl, xnetpp
#endif
  use mod_aux_pdaf, &
       only: field2state

  implicit none

! *** Arguments ***
  integer, intent(in) :: dim_p               !< PE-local state dimension
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector

! *** Local variables ***
  integer :: i                ! Counter


  ! *********************************
  ! Collect state vector 2d variables
  ! *********************************

  ! Note: The calls account for the halo offsets i0 and j0

  ! SSH
  if (id%ssh > 0) then
     call field2state(sshn(1+i0:ni_p+i0, 1+j0:nj_p+j0), &
          state_p, &
          sfields(id%ssh)%off, sfields(id%ssh)%ndims)
  end if


  ! *********************************
  ! Collect state vector 3d variables
  ! *********************************

  ! T
  if (id%temp > 0) then
     call field2state(tsn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_tem), &
          state_p, &
          sfields(id%temp)%off, sfields(id%temp)%ndims)
  end if

  ! S
  if (id%salt > 0) then
     call field2state(tsn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_sal), &
          state_p, &
          sfields(id%salt)%off, sfields(id%salt)%ndims)
  end if

  ! U
  if (id%uvel > 0) then
       call field2state(un(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
            state_p, &
            sfields(id%uvel)%off, sfields(id%uvel)%ndims)
  end if

  ! V
  if (id%vvel > 0) then
      call field2state(vn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
           state_p, &
           sfields(id%vvel)%off, sfields(id%vvel)%ndims)
  end if

#if defined key_top
  ! BGC
  do i = 1, jpbgc_prog
     if (sv_bgc_prog(i)) then
        call field2state(trn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, sfields(id%bgc_prog(i))%jptrc), &
             state_p, &
             sfields(id%bgc_prog(i))%off, sfields(id%bgc_prog(i))%ndims)
     end if
  end do

! Diagnostic variables are not initialized at the initial time
!   do i = 1, jpbgc_diag
!      if (sv_bgc_diag(i)) then
!         select case (i)
!         case (1)
!            call field2state(xpco2(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
!                 state_p, &
!                 sfields(id%bgc_diag(i))%off, sfields(id%bgc_diag(i))%ndims)
!         case (2)
!            call field2state(xph(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
!                 state_p, &
!                 sfields(id%bgc_diag(i))%off, sfields(id%bgc_diag(i))%ndims)
!         case (3)
!            call field2state(xchl(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
!                 state_p, &
!                 sfields(id%bgc_diag(i))%off, sfields(id%bgc_diag(i))%ndims)
!         case (4)
!            call field2state(xnetpp(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
!                 state_p, &
!                 sfields(id%bgc_diag(i))%off, sfields(id%bgc_diag(i))%ndims)
!         end select
!      end if
!   end do
#endif

end subroutine collect_state_init_pdaf
