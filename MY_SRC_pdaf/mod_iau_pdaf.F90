!> Module for using NEMO's ASMINC module with PDAF
!!
!! This module provides the interfacing with NEMO's
!! ASMINC module. Settings of ASMINC are set on the
!! basis of DA settings for PDAF. both the direct
!! initialization and IAU of ASMINC can be used.
!!
module mod_iau_pdaf

  use mod_kind_pdaf
  use par_oce, &
       only: jpi, jpj, jpk, jpkm1, jp_tem, jp_sal
  use mod_parallel_pdaf, &
       only: mype_ens
  use mod_nemo_pdaf, &
       only: nitend, nit000
  use mod_statevector_pdaf, &
       only: update_ssh, update_temp, update_salt, update_vel
  use asminc, &
       only: ln_bkgwri, ln_trainc, ln_dyninc, ln_sshinc, &
       ln_asmdin, ln_asmiau, nitbkg, nitdin, nitiaustr, nitiaufin, &
       niaufn, ln_salfix, salfixmin, nn_divdmp, &
       tra_asm_inc, dyn_asm_inc, ssh_asm_inc
  use asmpar, &
       only: nitdin_r, nitiaustr_r, nitiaufin_r
#if defined key_top
  use mod_statevector_pdaf, only: n_trc, jptra, sv_bgc1, sfields, id
  use asminc, &
       only: n_update_bgc, ids_update_bgc, bgc_bkginc, &
       nitdinbgc, nitibgcstr, nitibgcfin, niaufnbgc, &
       ln_trcinc, ln_bgcdin, ln_bgciau, &
       ln_oxyinc, ln_no3inc, ln_nh4inc, ln_po4inc, &
       ln_flainc, ln_diainc, ln_cyainc
  use asmpar, &
       only: nitdinbgc_r, nitibgcstr_r, nitibgcfin_r
#endif

   implicit none
   save
   public

   ! Variables that control the behavior of increments
   ! Next to these, the DA parameters set for PDAF influence the increments
   logical :: do_asmiau = .false.      ! Whether to apply IAU for ocean physics; if FALSE, DIN is used
   logical :: do_bgciau = .false.      ! Whether to apply IAU for BGC; if FALSE, DIN is used
   integer :: steps_asmiau = 1         ! Number of steps for physics IAU 
   integer :: steps_bgciau = 1         ! Number of steps for BGC IAU
   integer :: shape_asmiau = 0         ! Shape of physics IAU function: (0) constant; (1) hat (NEMO: niaufn)
   integer :: shape_bgciau = 0         ! Shape of BGC IAU function: (0) constant; (1) hat (NEMO: niaufnbgc)
   integer :: iter_divdmp = 0          ! Number of iterations of divergence damping operator

! *** Local variables ***
   integer, private :: next_inc = 0

contains

!> Routine to initialises variables for NEMO's ASMINC module
!!
!! The routine initializes variables for the ASMINC module
!! of NEMO. Without PDAF they would be read from a namelist.
!! The routine also initializes the increment array for the
!! BGC data assimilation uypdate.
!!
   subroutine asm_inc_init_pdaf(delt_obs)

     implicit none

! *** Arguments ***
     integer, intent(in) :: delt_obs

! *** Local variables
     integer :: id_bgc1
     integer :: id_var


! *******************************************
! *** Initialize Settings for NEMO ASMINC ***
! *******************************************

     ! Set initial increment step
     next_inc = delt_obs

     ! Set switches and parameters for ASMINC - physics

     ln_bkgwri  = .false.   !  Logical switch for writing out background state
     if (.not. (update_ssh .or. update_temp .or. update_salt .or. update_vel)) then
        ln_trainc  = .false.   !  Logical switch for applying tracer increments
        ln_dyninc  = .false.   !  Logical switch for applying velocity increments
        ln_sshinc  = .false.   !  Logical switch for applying SSH increments
        ln_asmdin  = .false.      !  Logical switch for Direct Initialization (DI)
        ln_asmiau  = .false.   !  Logical switch for Incremental Analysis Updating (IAU)
     else
        if (update_temp .or. update_salt) then
           ln_trainc  = .true.    !  Logical switch for applying tracer increments
        else
           ln_trainc  = .false.   !  Logical switch for applying tracer increments
        end if
        if (update_vel) then
           ln_dyninc  = .true.    !  Logical switch for applying velocity increments
        else
           ln_dyninc  = .false.   !  Logical switch for applying velocity increments
        end if
        if (update_ssh) then
           ln_sshinc  = .true.    !  Logical switch for applying SSH increments
        else
           ln_sshinc  = .false.   !  Logical switch for applying SSH increments
        end if
        if (do_asmiau) then
           ln_asmiau  = .true.    !  Logical switch for Incremental Analysis Updating (IAU)
           ln_asmdin  = .false.      !  Logical switch for Direct Initialization (DI)
        else
           ln_asmiau  = .false.   !  Logical switch for Incremental Analysis Updating (IAU)
           ln_asmdin  = .true.      !  Logical switch for Direct Initialization (DI)
        end if
     end if
     nitbkg     = 0         !  Timestep of background in [0,nitend-nit000-1]
     nitdin     = delt_obs  !  Timestep of background for DI in [0,nitend-nit000-1]
     nitiaustr = delt_obs               !  Timestep of start of BGC IAU interval
     nitiaufin = delt_obs+steps_asmiau  !  Timestep of end of BGC IAU interval
     niaufn     = shape_asmiau          !  Type of IAU weighting function
     ln_salfix  = .false.               !  Logical switch for ensuring that the sa > salfixmin
     salfixmin  = -9999                 !  Minimum salinity after applying the increments
     nn_divdmp  = iter_divdmp           !  Number of iterations of divergence damping operator


#if defined key_top

     ! Switches and parameters for ASM BGC increments

     ln_trcinc  = .true.        !  Logical switch for applying BGC increments

     if (do_bgciau) then        ! Switch between IAU and direct insertion
        ln_bgciau  = .true.     !  Logical switch for Incremental Analysis Updating (IAU)
        ln_bgcdin = .false.     !  Logical switch for applying DI for BGC
     else
        ln_bgciau  = .false.    !  Logical switch for Incremental Analysis Updating (IAU)
        ln_bgcdin = .true.      !  Logical switch for applying DI for BGC
     end if
     niaufnbgc  = shape_bgciau  !  Type of BGC IAU weighting function
     
     ! The steps are initialized here and later changed in cycled DA
     nitdinbgc = delt_obs                 !  Timestep for DI for BGC
     nitibgcstr = delt_obs                !  Timestep of start of BGC IAU interval
     nitibgcfin = delt_obs+steps_bgciau-1 !  Timestep of end of BGC IAU interval

     ! The folloing switches should be set for consistency
     ! they are not used with PDAF
     ln_oxyinc = .false.     !: No oxygen concentration increment
     ln_no3inc = .false.     !: No nitrate concentration increment
     ln_nh4inc = .false.     !: No ammonium concentration increment
     ln_po4inc = .false.     !: No phosphate concentration increment
     ln_flainc = .false.     !: No flagellate concentration increment
     ln_diainc = .false.     !: No diatom concentration increment
     ln_cyainc = .false.     !: No cyano concentration increment
#endif


! *************************************************
! *** Prepare increment array for BCG variables ***
! *************************************************

#if defined key_top

    ! Count number of prognostic BGC fields that are updated by the DA
    allocate(ids_update_bgc(jptra))
    do id_bgc1 = 1, jptra

       if (sv_bgc1(id_bgc1)) then

          id_var=id%bgc1(id_bgc1)

          if (sfields(id_var)%update) then
             n_update_bgc = n_update_bgc + 1
             ids_update_bgc(n_update_bgc) = id_bgc1
          end if
       end if
    end do

    ! Allocate increment array for BGC variables
    allocate(bgc_bkginc(jpi, jpj, jpk, n_update_bgc))
    bgc_bkginc = 0.0

    ! Update ln_trcinc if no BGC updates are done by PDAF
    if (n_update_bgc==0) ln_trcinc  = .false.    !  Logical switch for applying tracer increments

#endif


! *********************
! *** Screen output ***
! *********************

    if (mype_ens==0) then
       write (*,'(/a,4x,a)') 'NEMO-PDAF', '******** Setup for NEMO ASM ********'
       if (update_temp .or. update_salt .or. update_ssh .or. update_vel) then 
          write (*,'(a,4x,a)') 'NEMO-PDAF', '--- Apply increment for NEMO physics'
          if (update_temp) write (*,'(a,9x,a)') 'NEMO-PDAF', '--- update temperature'
          if (update_salt) write (*,'(a,9x,a)') 'NEMO-PDAF', '--- update salinity'
          if (update_ssh) write (*,'(a,9x,a)') 'NEMO-PDAF', '--- update SSH'
          if (update_vel) write (*,'(a,9x,a)') 'NEMO-PDAF', '--- update velocities'
          if (ln_asmdin) write (*,'(a,9x,a)') 'NEMO-PDAF', '--- Use DIN for NEMO fields'
          if (ln_asmiau) then
             write (*,'(a,8x,a)') 'NEMO-PDAF', '--- Use IAU for NEMO fields'
             write (*,'(a,8x,a,i9)') 'NEMO-PDAF', '--- Number of IAU steps', steps_asmiau
          end if
          if (nn_divdmp>0) &
               write (*,'(a,8x,a,i9)') 'NEMO-PDAF', '--- Apply divergence damping, iterations', nn_divdmp
       else
          write (*,'(a,4x,a)') 'NEMO-PDAF', '--- No increment for NEMO physics'
       end if

#if defined key_top
       if (ln_trcinc) then
          write (*,'(a,4x,a)') 'NEMO-PDAF', '--- Apply increment for BGC'
          write (*,'(a,4x,a,i5)') 'NEMO-PDAF', '--- Number of updated BCG fields:', n_update_bgc
          if (ln_bgciau) then
             write (*,'(a,8x,a)') 'NEMO-PDAF', '--- Use IAU for BGC fields'
             write (*,'(a,8x,a,i9)') 'NEMO-PDAF', '--- Number of BGC IAU steps', steps_bgciau
          else
             write (*,'(a,8x,a)') 'NEMO-PDAF', '--- Use DIN for BGC fields'
          end if

       else
          write (*,'(a,4x,a,i5)') 'NEMO-PDAF', '--- No increment for BCG fields'
       end if
#endif
       write (*,'(a,4x,a)') 'NEMO-PDAF', '************************************'
    end if

     ! Update increment step
     call update_asm_step_pdaf()

   end subroutine asm_inc_init_pdaf



!-------------------------------------------------------------------------------
!> Routine to store analysis step for NEMO-ASM
!!
!! For direct initialization, the increment is applied on the
!! time step after the analysis update is computed by PDAF
   subroutine store_asm_step_pdaf(nextinc)

     implicit none

     integer, intent(in) :: nextinc      ! Time of (starting) next assimilation increment

     ! Store step of (start) of next increments
     next_inc = nextinc + 1

     call update_asm_step_pdaf()

   end subroutine store_asm_step_pdaf



!-------------------------------------------------------------------------------
!> Routine to update analysis step values for NEMO-ASM
!!
!! This routine is called in assimilate_pdaf to
!! update the time step information for applying
!! assimilation increments.
!!
   subroutine update_asm_step_pdaf()

     implicit none

! *** Settings for ASMINC physics ***

     nitdin       = next_inc                           ! Time step of the background state for direct initialization
     nitdin_r     = nitdin    + nit000 - 1             ! Background time for DI referenced to nit000
     nitiaustr   = next_inc                            ! Timestep of start of IAU interval
     nitiaufin   = next_inc+steps_asmiau - 1           ! Timestep of end of IAU interval
     nitiaustr_r = nitiaustr + nit000 - 1              ! Start of IAU interval referenced to nit000
     nitiaufin_r = nitiaufin + nit000 - 1              ! End of IAU interval referenced to nit000

     if (mype_ens==0) then
        if (ln_trainc .or. ln_dyninc .or. ln_sshinc) then
           if (ln_asmiau) then
              write (*,'(a,5x,a,2i)') 'NEMO-PDAF', '--- set IAU steps for ASMINC: ', nitiaustr_r, nitiaufin_r
           else
              write (*,'(a,5x,a,i)') 'NEMO-PDAF', '--- set DIN step for ASMINC: ', nitdin_r
           end if
        end if
     end if


#if defined key_top
! *** Settings for ASMINC BGC ***

     nitdinbgc    = next_inc                           ! Time step of direct init for BGC
     nitdinbgc_r  = nitdinbgc    + nit000 - 1          ! Background time for DI referenced to nit000
     nitibgcstr   = next_inc                           ! Timestep of start of BGC IAU interval
     nitibgcfin   = next_inc+steps_bgciau - 1          ! Timestep of end of BGC IAU interval
     nitibgcstr_r = nitibgcstr + nit000 - 1            ! Start of BGC IAU interval referenced to nit000
     nitibgcfin_r = nitibgcfin + nit000 - 1            ! End of BGC IAU interval referenced to nit000

     if (mype_ens==0) then
        if (ln_trcinc) then
           if (ln_bgciau) then
              write (*,'(a,5x,a,2i)') 'NEMO-PDAF', '--- set BGC IAU steps for ASMINC: ', nitibgcstr_r, nitibgcfin_r
           else
              write (*,'(a,5x,a,i)') 'NEMO-PDAF', '--- set BGC DIN step for ASMINC: ', nitdinbgc_r
           end if
        end if
     end if
#endif

   end subroutine update_asm_step_pdaf



!-------------------------------------------------------------------------------
!> Routine to update the bkginc arrays of NEMO-ASM
!!
!! This routine updates the increment arrays for the
!! ASMINC module of NEMO.
!!
   subroutine update_bkginc_pdaf(dim_p, state_p, verbose)

     use mod_statevector_pdaf, &
          only: sfields, id
     use mod_nemo_pdaf, &
          only: ni_p, nj_p, nk_p, i0, j0, &
          jp_tem, jp_sal, lbc_lnk, lbc_lnk_multi, &
          sshn, tsn, un, vn, tmask
     use asminc, &
          only: ssh_bkginc, t_bkginc, s_bkginc, u_bkginc, v_bkginc, &
          ssh_bkg, t_bkg, s_bkg, u_bkg, v_bkg
#if defined key_top
     use mod_statevector_pdaf, &
          only: jptra, sv_bgc1
     use mod_nemo_pdaf, &
          only: trb, trn
#endif
     use mod_aux_pdaf, &
          only: state2field_inc

     implicit none

! *** Arguments ***
     integer, intent(in) :: dim_p               !< PE-local state dimension
     real(pwp), intent(inout) :: state_p(dim_p) !< PE-local state vector
     integer, intent(in) :: verbose             !< Verbosity flag

! *** Local variables ***
     integer :: i_bgcinc, i_trn     ! Indices
     integer :: id_var              ! Index


! *************************************************
! *** Prepare increment arrays for NEMO physics ***
! *************************************************

     physics: if (update_temp .or. update_salt .or. update_ssh .or. update_vel) then

        ! Ensure that the increment arrays are allocated and set to zero
        if (.not. allocated(ssh_bkginc)) allocate(ssh_bkginc(jpi,jpj))
        if (.not. allocated(t_bkginc)) allocate(t_bkginc(jpi,jpj,jpk))
        if (.not. allocated(s_bkginc)) allocate(s_bkginc(jpi,jpj,jpk))
        if (.not. allocated(u_bkginc)) allocate(u_bkginc(jpi,jpj,jpk))
        if (.not. allocated(v_bkginc)) allocate(v_bkginc(jpi,jpj,jpk))
        ssh_bkginc = 0.0_pwp
        t_bkginc = 0.0_pwp
        s_bkginc = 0.0_pwp
        u_bkginc = 0.0_pwp
        v_bkginc = 0.0_pwp

        ! SSH
        if (id%ssh > 0 .and. update_ssh) then
           call state2field_inc(state_p, sshn(1+i0:ni_p+i0, 1+j0:nj_p+j0), &
                ssh_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0), sfields(id%ssh)%off, sfields(id%ssh)%ndims)

           ! Fill halo regions
           call lbc_lnk('distribute_state_pdaf', ssh_bkginc, 'T', 1.)
        end if

        ! T
        if (id%temp > 0 .and. update_temp) then
           call state2field_inc(state_p, &
                tsn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_tem), &
                t_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                sfields(id%temp)%off, sfields(id%temp)%ndims)
        end if

        ! S
        if (id%salt > 0 .and. update_salt) then
           call state2field_inc(state_p, &
             tsn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_sal), &
             s_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
             sfields(id%salt)%off, sfields(id%salt)%ndims)
        end if

        if (id%temp>0 .or. id%salt>0) then
           ! Fill halo regions
           call lbc_lnk_multi('distribute_state_pdaf', t_bkginc, 'T', &
                1., s_bkginc, 'T', 1.)
        end if

        ! U
        if (id%uvel > 0 .and. update_vel) then
           call state2field_inc(state_p, &
                un(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                u_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                sfields(id%uvel)%off, sfields(id%uvel)%ndims)
        end if

        ! V
        if (id%vvel > 0 .and. update_vel) then
           call state2field_inc(state_p, &
                vn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                v_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                sfields(id%vvel)%off, sfields(id%vvel)%ndims)
        end if

        if (id%uvel>0 .or. id%vvel>0) then
           ! Fill halo regions
           call lbc_lnk_multi('distribute_state_pdaf', u_bkginc, 'U', -1., v_bkginc, 'V', -1.)

           ! Apply divergence damping to the velocity increment
           call div_damping_filter()
        end if


        ! For direct initialization we intialize the background arrays
        ! and apply the increments here. They are set the Euler flag
        if (ln_asmdin) then

           
           ! Ensure that the background arrays are allocated
           if (.not. allocated(ssh_bkg)) allocate(ssh_bkg(jpi,jpj))
           if (.not. allocated(t_bkg)) allocate(t_bkg(jpi,jpj,jpk))
           if (.not. allocated(s_bkg)) allocate(s_bkg(jpi,jpj,jpk))
           if (.not. allocated(u_bkg)) allocate(u_bkg(jpi,jpj,jpk))
           if (.not. allocated(v_bkg)) allocate(v_bkg(jpi,jpj,jpk))

           ! Initialize background fields for use in asm_inc routines
           ssh_bkg = sshn
           t_bkg = tsn(:,:,:,jp_tem)
           s_bkg = tsn(:,:,:,jp_sal)
           u_bkg = un
           v_bkg = vn
           
           if (mype_ens==0) &
                write (*,'(a,4x,a)') 'NEMO-PDAF', '--- Apply full increment in ASMINC'

           ! Apply assimilation increment
           IF( ln_trainc )   CALL tra_asm_inc(next_inc)      ! Tracers
           IF( ln_dyninc )   CALL dyn_asm_inc(next_inc)      ! Dynamics
           IF( ln_sshinc )   CALL ssh_asm_inc(next_inc)      ! SSH

        else
           if (mype_ens==0) &
                write (*,'(a,4x,a)') 'NEMO-PDAF', '--- Store increments for IAU'
        end if

     end if physics


! *************************************************
! *** Prepare increment array for BCG variables ***
! *************************************************

#if defined key_top
     do i_bgcinc = 1, n_update_bgc

        ! Get index in BGC tracer array
        i_trn = ids_update_bgc(i_bgcinc)

        if (sv_bgc1(i_trn)) then

           ! Get field index in state vector
           id_var=id%bgc1(i_trn)

           if (sfields(id_var)%update) then

              ! Update prognostic BGC variables
              call state2field_inc(state_p, &
                   trn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, i_trn), &
                   bgc_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, i_bgcinc), &
                   sfields(id_var)%off, sfields(id_var)%ndims)

              ! Fill halo regions
              call lbc_lnk('distribute_state_pdaf', bgc_bkginc(:, :, :, i_bgcinc), 'T', &
                   1._pwp)

              bgc_bkginc(:,:,:, i_bgcinc) = bgc_bkginc(:,:,:,i_bgcinc) * tmask(:,:,:)
           end if
        end if
     end do
#endif

   end subroutine update_bkginc_pdaf


!-------------------------------------------------------------------------------
!< Routine to apply divergence damping
!!
!! The functionality inside the routine is copied from NEMO's
!! routine asm_inc_init
!!
   subroutine div_damping_filter()

     use asminc, &
          only: u_bkginc, v_bkginc
     use par_oce, only: jpkm1, jpjm1, jpim1
     use dom_oce, only: e1v, e3v_n, e2u, e3u_n, e3t_n, &
          r1_e1u, r1_e2v, umask, vmask
     use lbclnk, only: lbc_lnk

     implicit none

! *** Local variables ***
     integer :: ji, jj, jt, jk                       ! Counters
     real(pwp), allocatable ::   zhdiv(:,:)   ! 2D workspace

   !! * Substitutions
#include "vectopt_loop_substitute.h90"

! ***************************************
! *** Apply divergence damping filter ***
! ***************************************

     divdmp: if ( ln_dyninc .and. nn_divdmp > 0 ) then    ! Apply divergence damping filter

        allocate( zhdiv(jpi,jpj) ) 

        do jt = 1, nn_divdmp

           do jk = 1, jpkm1           ! zhdiv = e1e1 * div
              zhdiv(:,:) = 0.0
              do jj = 2, jpjm1
                 do ji = fs_2, fs_jpim1   ! vector opt.
                    zhdiv(ji,jj) = (  e2u(ji  ,jj) * e3u_n(ji  ,jj,jk) * u_bkginc(ji  ,jj,jk)    &
                         &            - e2u(ji-1,jj) * e3u_n(ji-1,jj,jk) * u_bkginc(ji-1,jj,jk)    &
                         &            + e1v(ji,jj  ) * e3v_n(ji,jj  ,jk) * v_bkginc(ji,jj  ,jk)    &
                         &            - e1v(ji,jj-1) * e3v_n(ji,jj-1,jk) * v_bkginc(ji,jj-1,jk)  ) / e3t_n(ji,jj,jk)
                 end do
              end do
              call lbc_lnk( 'asminc', zhdiv, 'T', 1. )   ! lateral boundary cond. (no sign change)

              do jj = 2, jpjm1
                 do ji = fs_2, fs_jpim1   ! vector opt.
                    u_bkginc(ji,jj,jk) = u_bkginc(ji,jj,jk)                         &
                         + 0.2_pwp * ( zhdiv(ji+1,jj) - zhdiv(ji  ,jj) ) * r1_e1u(ji,jj) * umask(ji,jj,jk)
                    v_bkginc(ji,jj,jk) = v_bkginc(ji,jj,jk)                         &
                         + 0.2_pwp * ( zhdiv(ji,jj+1) - zhdiv(ji,jj  ) ) * r1_e2v(ji,jj) * vmask(ji,jj,jk) 
                 end do
              end do
           end do

        end do

        deallocate( zhdiv ) 

     endif divdmp

   end subroutine div_damping_filter

!-------------------------------------------------------------------------------
!> Routine to deallocate the BGC increment array
!!
   subroutine asm_inc_deallocate_pdaf

     implicit none
     
#if defined key_top
     if (allocated(bgc_bkginc)) deallocate(bgc_bkginc)
#endif

   end subroutine asm_inc_deallocate_pdaf


end module mod_iau_pdaf
