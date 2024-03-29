!$Id$
!>  Restrict a model state to a local analysis domain
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during the loop over all
!! local analysis domains in PDAF_X_update
!! before the analysis on a single local analysis 
!! domain.  It has to initialize elements of the 
!! state vector for the local analysis domains from
!! the PE-local full state vector.
!!
!! Generic implementation using index vector 
!! ID_LSTATE_IN_PSTATE.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
subroutine g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

   use mod_assimilation_pdaf, &
        only:  id_lstate_in_pstate

  implicit none

! *** Arguments ***
  integer, intent(in) :: step           !< Current time step
  integer, intent(in) :: domain_p       !< Current local analysis domain
  integer, intent(in) :: dim_p          !< PE-local full state dimension
  integer, intent(in) :: dim_l          !< Local state dimension
  real, intent(in)    :: state_p(dim_p) !< PE-local full state vector 
  real, intent(out)   :: state_l(dim_l) !< State vector on local analysis domain

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update    (as U_g2l_state)
! Called by: PDAF_letkf_update    (as U_g2l_state)
! Called by: PDAF_lestkf_update   (as U_g2l_state)
! Called by: PDAF_lnetf_update    (as U_g2l_state)
!EOP

! *** local variables *** 
  integer(4) :: i 

  
! *************************************
! *** Initialize local state vector ***
! *************************************
  
  ! Generic initialization using ID_LSTATE_IN_PSTATE set in INIT_DIM_L_PDAF
  do i = 1, dim_l
     state_l(i) = state_p(id_lstate_in_pstate(i))
  enddo

end subroutine g2l_state_pdaf
