!$Id$
!BOP
!
! !ROUTINE: finalize_pdaf --- Finalize PDAF
!
! !INTERFACE:
SUBROUTINE finalize_pdaf()

! !DESCRIPTION:
! This routine call MPI_finalize
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, &
       ONLY: finalize_parallel, mype_world, npes_world

  IMPLICIT NONE    
  
! !CALLING SEQUENCE:
! Called by: main program
!EOP

! *** Show allocated memory for PDAF ***
  IF (mype_world==0) CALL PDAF_print_info(2)
  IF (npes_world>1) CALL PDAF_print_info(11)  ! Globally allocated memory, needs PDAF >2.0

! *** Print PDAF timings onto screen ***
  IF (mype_world==0) CALL PDAF_print_info(3)

! *** Deallocate PDAF arrays ***
  CALL PDAF_deallocate()

! *** Finalize parallel MPI region - if not done by model ***
!  CALL finalize_parallel()

END SUBROUTINE finalize_pdaf
