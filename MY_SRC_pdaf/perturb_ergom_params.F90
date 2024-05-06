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


subroutine param_scale(stddev, scaletype)
    
  use mod_kind_pdaf
  use mod_parallel_pdaf, &
       only: mype_model, task_id

#if defined key_top
  use bioparam, &
       only: alphap, alphaf, alphab, min_odial, min_oflagl, min_ocyanl, &
       wdz, wpz, wbz, mezgraz, mizgraz, rp0, rf0, rb0, q10_recs, &
       sigma_b, deltao, nb, lnmez, lnmiz
#endif

  implicit none

! *** Arguments ***
  real(pwp), intent(in) :: stddev    ! Relative standard deviation for lognormal perturbation
  integer, intent(in) :: scaletype   ! Type of parameter scaling

! *** Local variables ***
  real(pwp) :: istddev, rstddev


! *************************************************
! *** Perturb model parameters for ensemble run ***
! *************************************************

#if defined key_top

  istddev = 1.0 / (1.0 + stddev)
  rstddev = 1.0 + stddev

  if (mype_model==0 .and. task_id==1) then
     write (*,*) 'NEMO-PDAF ', 'Apply scaling factors to ERGOM parameters: '
     write (*,*) 'NEMO-PDAF ', '     alphap, alphaf, alphab, '
     write (*,*) 'NEMO-PDAF ', '     wdz, wpz, wbz, mezgraz, mizgraz, rp0, rf0, rb0'
     write (*,*) 'NEMO-PDAF ', 'scaling factor to multiply or divide', rstddev
     write (*,'(a,a,11f8.3)') 'NEMO-PDAF ', 'Original parameter values', alphap, alphaf, alphab, &
          wdz, wpz, wbz, mezgraz, mizgraz, rp0, rf0, rb0
     if (scaletype==2) then
        write (*,*) 'NEMO-PDAF ', '     sigma_b, deltao, nb, '
        write (*,*) 'NEMO-PDAF ', '     lnmez, lnmiz'
        write (*,*) 'NEMO-PDAF ', 'scaling factor to multiply or divide', rstddev
        write (*,'(a,a,11f8.3)') 'NEMO-PDAF ', 'Original parameter values', sigma_b, &
             deltao, nb, lnmez, lnmiz
     end if
  end if

  ! Apply perturbations
  if (scaletype==1) then
     if (task_id==1) then
     ! original parameters
     end if
     if (task_id==2) alphap = istddev * alphap
     if (task_id==3) alphap = rstddev * alphap
     if (task_id==4) alphaf = istddev * alphaf
     if (task_id==5) alphaf = rstddev * alphaf
     if (task_id==6) alphab = istddev * alphab
     if (task_id==7) alphab = rstddev * alphab
     if (task_id==8) wdz = istddev * wdz
     if (task_id==9) wdz = rstddev * wdz
     if (task_id==10) wpz = istddev * wpz
     if (task_id==11) wpz = rstddev * wpz
     if (task_id==12) wbz = istddev * wbz
     if (task_id==13) wbz = rstddev * wbz
     if (task_id==14) mezgraz = istddev * mezgraz
     if (task_id==15) mezgraz = rstddev * mezgraz
     if (task_id==16) mizgraz = istddev * mizgraz
     if (task_id==17) mizgraz = rstddev * mizgraz
     if (task_id==18) rp0 = istddev * rp0
     if (task_id==19) rp0 = rstddev * rp0
     if (task_id==20) rf0 = istddev * rf0
     if (task_id==21) rf0 = rstddev * rf0
     if (task_id==22) rb0 = istddev * rb0
     if (task_id==23) rb0 = rstddev * rb0
     if (task_id==24) then
        alphap = istddev * alphap
        alphaf = istddev * alphab
        alphab = istddev * alphaf
     end if
     if (task_id==25) then
        alphap = rstddev * alphap
        alphaf = rstddev * alphab
        alphab = rstddev * alphaf
     end if
     if (task_id==26) then
        rp0 = istddev * rp0
        rf0 = istddev * rf0
        rb0 = istddev * rb0
     end if
     if (task_id==27) then
        rp0 = rstddev * rp0
        rf0 = rstddev * rf0
        rb0 = rstddev * rb0
     end if
     if (task_id==28) then
        wdz = istddev * wdz
        wpz = istddev * wpz
        wbz = istddev * wbz
     end if
     if (task_id==29) then
        wdz = rstddev * wdz
        wpz = rstddev * wpz
        wbz = rstddev * wbz
     end if
     if (task_id==30) then
        mezgraz = rstddev * mezgraz
        mizgraz = rstddev * mizgraz
     end if
     if (mype_model==0) &
          write (*,'(a,i4,1x,a,11f8.3)') 'NEMO-PDAF task', task_id, 'scaled parameters:', alphap, alphaf, alphab, &
          wdz, wpz, wbz, mezgraz, mizgraz, rp0, rf0, rb0
  elseif (scaletype==2) then
     if (task_id==1) then
     ! Use parameter values from 1D model
        rp0 = 1.275
        rb0 = 0.682
        rf0 = 0.359
        alphap = 0.268
        alphaf = 0.102
        deltao = 0.023
        sigma_b = 0.03
        nb = 0.011
        mizgraz = 0.567
        if (mype_model==0) &
             write (*,'(a,i4,1x,a,11f8.3)') 'NEMO-PDAF task', task_id, 'scaled parameters:', alphap, alphaf, alphab, &
             wdz, wpz, wbz, mezgraz, mizgraz, rp0, rf0, rb0
     end if
     if (task_id==2) sigma_b = istddev * sigma_b
     if (task_id==3) sigma_b = rstddev * sigma_b
     if (task_id==4) deltao = istddev * deltao
     if (task_id==5) deltao = rstddev * deltao
     if (task_id==6) nb = istddev * nb
     if (task_id==7) nb = rstddev * nb
     if (task_id==8) lnmez = istddev * lnmez
     if (task_id==9) lnmez = rstddev * lnmez
     if (task_id==10) lnmiz = istddev * lnmiz
     if (task_id==11) lnmiz = rstddev * lnmiz
     if (mype_model==0) &
          write (*,'(a,a,11f8.3)') 'NEMO-PDAF ', 'Scaled parameters sigma', sigma_b, &
          deltao, nb, lnmez, lnmiz
  elseif (scaletype==3) then
     if (task_id==2) then
     ! Use parameter values from 1D model
        rp0 = 1.385 !1.554
        rb0 = 0.679 !0.522
        rf0 = 0.379 !0.467
        alphap = 0.248 !0.303
        alphab = 0.413 !0.559
        alphaf = 0.105 !0.128
        deltao = 0.023
        sigma_b = 0.029 !0.025
        nb = 0.011 !0.012
        mizgraz = 0.539 !0.624
     end if
     if (task_id==3) then
        mezgraz = istddev * mezgraz
        mizgraz = istddev * mizgraz
     end if
     if (mype_model==0) &
          write (*,'(a,i4,1x,a,11f8.3)') 'NEMO-PDAF task', task_id, 'scaled parameters:', alphap, alphaf, alphab, &
          wdz, wpz, wbz, mezgraz, mizgraz, rp0, rf0, rb0
  end if

#endif

end subroutine param_scale
