!/*@@
!  @file      timer.f90
!  @date      November 2000
!  @author    Lars Nerger
!  @desc 
!             F90 module for timing purposes
!  @enddesc 
!@@*/

module timer
! ******************************************************
! ***          Module to perform timing.             ***
! ***  It uses the intrinsic subroutine SYSTEM_CLOCK ***
! ***                                                ***
! *** The module provides a subroutine and two       ***
! *** functions to time regions and to access this   ***
! *** timing information.                            ***
! ***                                                ***
! *** SUBROUTINE TIMEIT(integer ID,char(len=3) flag) ***
! ***    Usage:                                      ***
! ***    CALL TIMEIT(N,'ini')                        ***
! ***         Initializes N counters                 ***
! ***                                                ***
! ***    CALL TIMEIT(M,'new')                        ***
! ***         Start timing region for counter M      ***
! ***                                                ***
! ***    CALL TIMEIT(M,'old')                        ***
! ***         End timing region for counter M        ***
! ***                                                ***
! *** REAL FUNCTION TIME_TEMP(integer ID)            ***
! ***    Usage:                                      ***
! ***    time = time_temp(M)                         ***
! ***         Duration of the last region of         ***
! ***         counter M in seconds                   ***
! ***                                                ***
! *** REAL FUNCTION TIME_TOT(integer ID)             ***
! ***    Usage:                                      ***
! ***    time = time_tot(M)                          ***
! ***         Duration of the all regions of         ***
! ***         counter M in seconds                   ***
! ******************************************************
  
  implicit none
  save

  private
  integer :: t_rate
  integer,allocatable,dimension(:) :: t_start,t_end
  real,allocatable,dimension(:) :: t_total,t_temp
  
  public :: timeit
  public :: time_tot,time_temp

  contains
    subroutine timeit(timerID,flag)

      implicit none

      integer :: timerID
      character(len=3) :: flag

      if (flag == 'ini') then
        if (.not.(allocated(t_start))) then
          allocate(t_start(timerID),t_end(timerID))
          allocate(t_total(timerID),t_temp(timerID))
        end if
        
        t_total = 0.0
      end if

      if (flag == 'new') then
        call system_clock(t_start(timerID))
      end if

      if (flag == 'old') then
        call system_clock(t_end(timerID),t_rate)
        t_temp(timerID) = float(t_end(timerID)-t_start(timerID)) &
             /float(t_rate)
        t_total(timerID) = t_total(timerID) + float(t_end(timerID) - &
             t_start(timerID)) / float(t_rate)
      end if

    end subroutine timeit

    real function time_temp(timerID)

      integer :: timerID

      time_temp = t_temp(timerID)
    end function time_temp

    real function time_tot(timerID)

      integer :: timerID

      time_tot = t_total(timerID)
    end function time_tot


end module timer
