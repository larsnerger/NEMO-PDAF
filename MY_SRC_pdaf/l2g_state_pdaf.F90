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
subroutine l2g_state_pdaf(step, domain_p_all, dim_l, state_l, dim_p, state_p)

  use mod_kind_pdaf
  use mod_assimilation_pdaf, &
       only: id_lstate_in_pstate, type_sweep, isweep, cda_phy, cda_bio
  use mod_statevector_pdaf, &
       only: n_fields, sfields, sfields_l
  use mod_nemo_pdaf, &
       only: nwet
use mod_parallel_pdaf, only: mype_filter
  implicit none

! *** Arguments ***
  integer, intent(in) :: step                !< Current time step
  integer, intent(in) :: domain_p_all        !< Current local analysis domain
  integer, intent(in) :: dim_l               !< Local state dimension
  integer, intent(in) :: dim_p               !< PE-local full state dimension
  real(pwp), intent(in)    :: state_l(dim_l) !< State vector on local analysis domain
  real(pwp), intent(inout) :: state_p(dim_p) !< PE-local full state vector 

! *** local variables *** 
  integer :: i, ifield           ! Counters
  logical :: update_cda          ! Whether to perform strongly-coupled DA update


! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************

  do ifield = 1, n_fields

     ! Determine whether to apply update according to CDA settings
     if (type_sweep(isweep)==sfields(ifield)%type) then
        ! Type of sweep equal to variable type
        update_cda = .true.
     else
        ! Check for strgonly coupled DA configuration
        if (type_sweep(isweep)=='phy' .and. trim(cda_phy)=='strong') then
           update_cda = .true.
        elseif (type_sweep(isweep)=='bio' .and. trim(cda_bio)=='strong') then
           update_cda = .true.
        else
           update_cda = .false.
        end if
     end if

     ! Update fields which are active and fulfill the CDA criteria
     if (sfields(ifield)%update .and. update_cda) then
        do i = sfields_l(ifield)%off+1, sfields_l(ifield)%off + sfields_l(ifield)%dim
           state_p(id_lstate_in_pstate(i)) = state_l(i)
        end do
     end if

  end do

end subroutine l2g_state_pdaf
