subroutine param_perturb(stddev)
    
  use mod_kind_pdaf
  use mod_parallel_pdaf, &
       only: mype_model, task_id

#if defined key_top
  use bioparam, &
       only: alphap, alphaf, alphab, min_odial, min_oflagl, min_ocyanl, &
       wdz, wpz, wbz, mezgraz, mizgraz, rp0, rf0, rb0, q10_recs
#endif

  implicit none

! *** Arguments ***
  real(pwp), intent(in) :: stddev    ! Relative standard deviation for lognormal perturbation

! *** Local variables ***
  integer, save :: iseed(4)


! *************************************************
! *** Perturb model parameters for ensemble run ***
! *************************************************

#if defined key_top

  ! Set seed
  iseed(1)=1
  iseed(2)=25
  iseed(3)=17+5*task_id
  iseed(4)=9

  if (mype_model==0 .and. task_id==1) then
     write (*,*) 'NEMO-PDAF ', 'Perturb ERGOM parameters: '
     write (*,*) 'NEMO-PDAF ', '     alphap, alphaf, alphab, min_odial, min_oflagl, min_ocyanl, '
     write (*,*) 'NEMO-PDAF ', '     wdz, wpz, wbz, mezgraz, mizgraz, rp0, rf0, rb0, q10_recs'
     write (*,*) 'NEMO-PDAF ', 'relative stddev for lognormal perturbations', stddev
     write (*,*) 'NEMO-PDAF ', 'Original parameter values', alphap, alphaf, alphab, min_odial, min_oflagl, min_ocyanl, &
       wdz, wpz, wbz, mezgraz, mizgraz, rp0, rf0, rb0, q10_recs
  end if

  ! Apply perturbations
  call perturb_lognormal(alphap, stddev, iseed)     ! Half saturation constant dia
  call perturb_lognormal(alphaf, stddev, iseed)     ! Half saturation constant flag
  call perturb_lognormal(alphab, stddev, iseed)     ! Half saturation constant cya
  call perturb_lognormal(rp0, stddev, iseed)        ! Max uptake rate at T0 dia
  call perturb_lognormal(rf0, stddev, iseed)        ! Max uptake rate at T0 flag
  call perturb_lognormal(rb0, stddev, iseed)        ! Max uptake rate at T0 cya
  call perturb_lognormal(min_odial, stddev, iseed)  ! Min opt. dia light W/m^2
  call perturb_lognormal(min_oflagl, stddev, iseed) ! Min opt. flag light W/m^2
  call perturb_lognormal(min_ocyanl, stddev, iseed) ! Min opt. cyan light W/m^2
  if (wdz<0.0) then
     wdz = -wdz
     call perturb_lognormal(wdz, stddev, iseed)     ! Sinking vel. detritus m/d
     wdz = - wdz
  end if
  if (wpz<0.0) then
     wpz = -wpz
     call perturb_lognormal(wpz, stddev, iseed)     ! Sinking vel. dia [m/d]
     wpz = - wpz
  end if
  call perturb_lognormal(wbz, stddev, iseed)        ! Rising vel. cyano [m/d]
  call perturb_lognormal(mezgraz, stddev, iseed)    ! Mesozoopl. grazing constant
  call perturb_lognormal(mizgraz, stddev, iseed)    ! Microzoopl. grazing constant
  call perturb_lognormal(q10_recs, stddev, iseed)   ! Recycling temp dep sed.[/C]

  if (mype_model==0) &
       write (*,*) 'NEMO-PDAF ', 'Perturbed parameters:', alphap, alphaf, alphab, min_odial, min_oflagl, min_ocyanl, &
       wdz, wpz, wbz, mezgraz, mizgraz, rp0, rf0, rb0, q10_recs
#endif

end subroutine param_perturb



subroutine perturb_lognormal(value, stddev, iseed)

  use mod_kind_pdaf

  implicit none

! *** Arguments ***
  real(pwp), intent(inout) :: value       ! value to be perturbed
  real(pwp), intent(in)    :: stddev      ! Standard deviation of lognormal distribution
  integer, intent(in)      :: iseed(4)    ! Seed for dlarnv

! *** Local variables ***
  real(pwp) :: sigma2     ! Variance
  real(pwp) :: logval     ! Logrithmic mean value of input value
  real(pwp) :: rndval     ! Normal random value

  ! Generate random number
  CALL dlarnv(3, iseed, 1, rndval)

  sigma2 = log(1.0d0 + stddev*stddev)
  logval = log(value) -0.5d0 * sigma2

  value = exp(logval + sqrt(sigma2) * rndval) ! Eq. (A10) from Ciavatta et al. (2016)

end subroutine perturb_lognormal
