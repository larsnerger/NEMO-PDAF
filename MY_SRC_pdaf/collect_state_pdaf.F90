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
       only: mype_model, task_id
  use mod_statevector_pdaf, &
       only: sfields, id
  use mod_nemo_pdaf, &
       only: ni_p, nj_p, nk_p, i0, j0, &
       jp_tem, jp_sal, ndastp, &
       sshn, tsn, un, vn
#if defined key_top
  use mod_statevector_pdaf, &
       only: n_trc, n_bgc1, n_bgc2, &
       jptra, jptra2, sv_bgc1, sv_bgc2
  use mod_nemo_pdaf, &
       only: trn, xph, xpco2, xchl
  use mod_assimilation_pdaf, &
       only: netppsum
#endif
  use mod_aux_pdaf, &
       only: field2state, transform_field_mv

  implicit none

! *** Arguments ***
  integer, intent(in) :: dim_p               !< PE-local state dimension
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector

! *** Local variables ***
  integer :: i                ! Counter
  integer :: verbose          ! Control verbosity


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
  do i = 1, jptra
     if (sv_bgc1(i)) then
        call field2state(trn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, sfields(id%bgc1(i))%jptrc), &
             state_p, &
             sfields(id%bgc1(i))%off, sfields(id%bgc1(i))%ndims)
     end if
  end do

  do i = 1, jptra2
     if (sv_bgc2(i)) then
        select case (i)
        case (1)
           call field2state(xpco2(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                state_p, &
                sfields(id%bgc2(i))%off, sfields(id%bgc2(i))%ndims)
        case (2)
           call field2state(xph(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                state_p, &
                sfields(id%bgc2(i))%off, sfields(id%bgc2(i))%ndims)
        case (3)
           call field2state(xchl(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                state_p, &
                sfields(id%bgc2(i))%off, sfields(id%bgc2(i))%ndims)
        case (4)
           ! For netpp we use the daily sum computed in assimilate_pdaf
           call field2state(netppsum(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                state_p, &
                sfields(id%bgc2(i))%off, sfields(id%bgc2(i))%ndims)
        end select
     end if
  end do
#endif

  ! Aply field transformations
  if (mype_model==0 .and. task_id==1) then
     verbose = 1
  else
     verbose = 0
  end if

  call transform_field_mv(1, state_p, 0, verbose)

end subroutine collect_state_pdaf
