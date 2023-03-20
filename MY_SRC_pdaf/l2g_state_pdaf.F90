!> Initialize full state from local analysis
!!
!! The routine is called during the loop over all
!! local analysis domains in `PDAF_X_update`
!! after the analysis and ensemble transformation
!! on a single local analysis domain. It has to
!! initialize elements of the PE-local full state
!! vector from the provided analysis state vector
!! on the local analysis domain.
!!
subroutine l2g_state_pdaf(step, domain_p, dim_l, state_l, dim_p, state_p)

  use mod_kind_pdaf
  use mod_assimilation_pdaf, &
       only: id_lstate_in_pstate
  use mod_statevector_pdaf, &
       only: n_fields, sfields, sfields_l

  implicit none

! *** Arguments ***
  integer, intent(in) :: step                !< Current time step
  integer, intent(in) :: domain_p            !< Current local analysis domain
  integer, intent(in) :: dim_l               !< Local state dimension
  integer, intent(in) :: dim_p               !< PE-local full state dimension
  real(pwp), intent(in)    :: state_l(dim_l) !< State vector on local analysis domain
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local full state vector 

! *** local variables *** 
  integer :: i, ifield                       ! Counters


! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************

  ! This generic initialization uses id_lstate_in_pstate
  ! and ifields_l which have been initialized in init_dim_l_pdaf

  do ifield = 1, n_fields

     ! Update DA-active variables
     if (sfields(ifield)%update) then
        do i = sfields_l(ifield)%off+1, sfields_l(ifield)%off + sfields_l(ifield)%dim
           state_p(id_lstate_in_pstate(i)) = state_l(i)
        end do
     end if

  end do

end subroutine l2g_state_pdaf
