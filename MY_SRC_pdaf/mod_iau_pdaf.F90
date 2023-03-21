!>##Using the Incremental Analysis Method
!>The material in this module is **heavily** based on a similar
!>implementation for the NEMOVAR data assimilation system. Please
!>refer to the `ASM` subdirectory in the `OCE` source code for
!>precise details.
!>
module mod_iau_pdaf

  use mod_kind_pdaf
  use par_oce, &
       only: jpi, jpj, jpk, jpkm1, jp_tem, jp_sal
  use mod_parallel_pdaf, &
       only: mype_ens
  use mod_nemo_pdaf, &
       only: nitend, nit000
  use mod_assimilation_pdaf, &
       only: delt_obs
  use mod_statevector_pdaf, &
       only: update_phys
  use asminc, &
       ONLY: ln_bkgwri, ln_trainc, ln_dyninc, ln_sshinc, &
       ln_asmdin, ln_asmiau, nitbkg, nitdin, nitiaustr, nitiaufin, &
       niaufn, ln_salfix, salfixmin, nn_divdmp
#if defined key_top
  use mod_statevector_pdaf, only: n_trc, jptra, sv_bgc1, sfields, id
  use asminc, &
       only: n_update_bgc, ids_update_bgc, bgc_bkginc, nitdinbgc, &
       ln_trcinc, ln_bgcdi, nitdinbgc, ln_oxyinc, ln_no3inc, &
       ln_nh4inc, ln_po4inc, ln_flainc, ln_diainc, ln_cyainc, &
       nitibgcstr, nitibgcfin
  use asmpar, &
       only: nitdinbgc_r
#endif

   implicit none
   save

! *** Local variables ***
   integer :: next_inc = 0

contains

!> Routine to initialises variables for NEMO's ASMINC module
!!
!! The routine initializes variables for the ASMINC module
!! of NEMO. Without PDAF they would be read from a namelist.
!! The routine also initializes the increment array for the
!! BGC data assimilation uypdate.
!!
   subroutine asm_inc_init_pdaf

     implicit none

     integer :: id_bgc1
     integer :: id_var

     ! Set switches and parameters for ASMINC
     ln_bkgwri  = .false.   !  Logical switch for writing out background state
     if (.not. update_phys) then
        ln_trainc  = .false.   !  Logical switch for applying tracer increments
        ln_dyninc  = .false.   !  Logical switch for applying velocity increments
        ln_sshinc  = .false.   !  Logical switch for applying SSH increments
     else
        ln_trainc  = .true.    !  Logical switch for applying tracer increments
        ln_dyninc  = .true.    !  Logical switch for applying velocity increments
        ln_sshinc  = .true.    !  Logical switch for applying SSH increments
     end if
     ln_asmdin  = .false.   !  Logical switch for Direct Initialization (DI)
     ln_asmiau  = .false.   !  Logical switch for Incremental Analysis Updating (IAU)

     nitbkg     = 0         !  Timestep of background in [0,nitend-nit000-1]
     nitdin     = delt_obs  !  Timestep of background for DI in [0,nitend-nit000-1]
     nitiaustr  = 1         !  Timestep of start of IAU interval in [0,nitend-nit000-1]
     nitiaufin  = 40        !  Timestep of end of IAU interval in [0,nitend-nit000-1]
     niaufn     = 0         !  Type of IAU weighting function
     ln_salfix  = .false.   !  Logical switch for ensuring that the sa > salfixmin
     salfixmin  = -9999     !  Minimum salinity after applying the increments
     nn_divdmp  = 0         !  Number of iterations of divergence damping operator

#if defined key_top

     ! Switches and settings for BGC increments
     ln_trcinc  = .true.    !  Logical switch for applying BGC increments
     ln_bgcdi = .true.      !  Logical switch for applying DI for BGC
     nitdinbgc = delt_obs   !  Timestep for DI for BGC
!     nitibgcstr = 1         !  Timestep of start of BGC IAU interval
!     nitibgcfin = 3         !  Timestep of end of BGC IAU interval

     ! The folloing switches should be set for consistency
     ! they are not used with PDAF
     ln_oxyinc = .false.     !: No oxygen concentration increment
     ln_no3inc = .false.     !: No nitrate concentration increment
     ln_nh4inc = .false.     !: No ammonium concentration increment
     ln_po4inc = .false.     !: No phosphate concentration increment
     ln_flainc = .false.     !: No flagellate concentration increment
     ln_diainc = .false.     !: No diatom concentration increment
     ln_cyainc = .false.     !: No cyano concentration increment


! *************************************************
! *** Prepare increment array for BCG variables ***
! *************************************************

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

#endif

    if (mype_ens==0) then
       write (*,'(/a,2x,a)') 'NEMO-PDAF', '*** Setup for NEMO ASM ***'
       if (update_phys) then 
          write (*,'(/a,5x,a)') 'NEMO-PDAF', '--- Apply increment for NEMO physics'
       else
          write (*,'(/a,5x,a)') 'NEMO-PDAF', '--- No increment for NEMO physics'
       end if
#if defined key_top
       write (*,'(a,5x,a,i5)') 'NEMO-PDAF', '--- Number of updated BCG fields:', n_update_bgc
       if (ln_bgcdi) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Use DI for BGC fields'
#endif
       if (ln_asmdin) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Use DI for NEMO fields'
       if (ln_asmiau) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Use IAU for NEMO fields'
       write (*,'(/a)') ' '

 write (*,*) 'NEMO-PDAF', 'TRC-IDs of updated BGC fields', ids_update_bgc(1:n_update_bgc)
    end if

   end subroutine asm_inc_init_pdaf



!> Routine to store analysis step for NEMO-ASM
!!
!! For direct initialization, the increment is applied on the
!! time step after the analysis update is computed by PDAF
   subroutine store_asm_step_pdaf(nextinc)

     implicit none

     integer, intent(in) :: nextinc

     next_inc = nextinc + 1

   end subroutine store_asm_step_pdaf


!> Routine to update analysis step for NEMO-ASM
!!
!! This routine is called in assimilate_pdaf to
!! update the time step information for applying
!! assimilation increments
   subroutine update_asm_step_pdaf()

     implicit none

     nitdin = next_inc
     nitdinbgc = next_inc
     nitdinbgc_r    = nitdinbgc    + nit000 - 1            ! Background time for DI referenced to nit000
!     nitibgcstr_r = nitibgcstr + nit000 - 1            ! Start of IAU interval referenced to nit000
!     nitibgcfin_r = nitibgcfin + nit000 - 1            ! End of IAU interval referenced to nit000
!     iiauperbgc   = nitibgcfin_r - nitibgcstr_r + 1     ! IAU interval length

   end subroutine update_asm_step_pdaf



!> Routine to update the bkginc arrays of NEMO-ASM
   subroutine update_bkginc_pdaf(dim_p, state_p, verbose)

     use mod_statevector_pdaf, &
          only: sfields, id
     use mod_nemo_pdaf, &
          only: ni_p, nj_p, nk_p, i0, j0, &
          jp_tem, jp_sal, lbc_lnk, lbc_lnk_multi, &
          sshn, tsn, un, vn, tmask
     use asminc, &
          only: ssh_bkginc, t_bkginc, s_bkginc, u_bkginc, v_bkginc
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
     integer :: id_var              ! Indices


     physics: if (update_phys > 0) then

        ! SSH
        if (id%ssh > 0) then
           call state2field_inc(state_p, sshn(1+i0:ni_p+i0, 1+j0:nj_p+j0), &
                ssh_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0), sfields(id%ssh)%off, sfields(id%ssh)%ndims)

           ! Fill halo regions
           call lbc_lnk('distribute_state_pdaf', ssh_bkginc, 'T', 1.)
        end if

        ! T
        if (id%temp > 0) then
           call state2field_inc(state_p, &
                tsn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p, jp_tem), &
                t_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                sfields(id%temp)%off, sfields(id%temp)%ndims)
        end if

        ! S
        if (id%salt > 0) then
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
        if (id%uvel > 0) then
           call state2field_inc(state_p, &
                un(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                u_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                sfields(id%uvel)%off, sfields(id%uvel)%ndims)
        end if

        ! V
        if (id%vvel > 0) then
           call state2field_inc(state_p, &
                vn(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                v_bkginc(1+i0:ni_p+i0, 1+j0:nj_p+j0, 1:nk_p), &
                sfields(id%vvel)%off, sfields(id%vvel)%ndims)
        end if

        if (id%uvel>0 .or. id%vvel>0) then
           ! Fill halo regions
           call lbc_lnk_multi('distribute_state_pdaf', u_bkginc, 'U', -1., v_bkginc, 'V', -1.)

        end if

     end if physics

     ! BGC
#if defined key_top
     do i_bgcinc = 1, n_update_bgc

        ! Get index in BGC tracer array
        i_trn = ids_update_bgc(i_bgcinc)

        if (sv_bgc1(i_trn)) then

           ! Get field index in state vector
           id_var=id%bgc1(i_trn)

           if (sfields(id_var)%update) then
     if (verbose==1) write (*,*) 'NEMO-PDAF', 'increment IDs', i_bgcinc, i_trn, id_var, sfields(id_var)%variable

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


   subroutine bgc3d_asm_inc( kt )
!       use dom_oce, only: rdt
!       use lib_mpp, only: ctl_stop
!       USE trc, ONLY:           & ! passive tracer variables
!       & trn,                &
!       & trb
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dyn_asm_inc  ***
      !!
      !! ** Purpose : Apply generic 3D biogeochemistry assimilation increments.
      !!
      !! ** Action  :
      !!----------------------------------------------------------------------
      integer,  intent(IN) :: kt      ! Current time step
      !
      integer   :: it              ! Index
      integer   :: jptrc
      integer   :: j, id_var
      real(pwp) :: zincwgt         ! IAU weight for current time step
      !!----------------------------------------------------------------------

      !--------------------------------------------------------------------
      ! Incremental Analysis Updating
      !--------------------------------------------------------------------

!       IF (kt <= nitiaufin  .AND. kt >= nitiaustr) THEN
! 
!          it = kt - nitiaustr + 1
!          zincwgt = wgtiau(it)   ! IAU weight for the current time step
!          ! note this is not a tendency so should not be divided by rdt
! 
!          ! Update the 3D BGC variables
!          ! Add directly to trn and trb, rather than to tra, because tra gets
!          ! reset to zero at the start of trc_stp, called after this routine
!          ! Don't apply increments if they'll take concentrations negative
! 
! #if defined key_top
!          do j = 1, jptra
!            if (sv_bgc1(j)) then
!             id_var = id%bgc1(j)
!             jptrc = sfields(id_var)%jptrc
!             WHERE( bgc_iau_pdaf(:,:,:,j) > 0.0_pwp .OR. &
!                  & trn(:,:,:,jptrc) + bgc_iau_pdaf(:,:,:,j) * zincwgt > 0.0_pwp )
!                trn(:,:,:,jptrc) = trn(:,:,:,jptrc) + bgc_iau_pdaf(:,:,:,j) * zincwgt
!                trb(:,:,:,jptrc) = trb(:,:,:,jptrc) + bgc_iau_pdaf(:,:,:,j) * zincwgt
!             END WHERE
!           end if
!          end do
! #else
!          CALL ctl_stop ( ' bgc3d_asm_inc: no compatible BGC model defined' )
! #endif
! 
!       ENDIF
      !
   end subroutine bgc3d_asm_inc

!> Routine to deallocate for BGC increment array
!!
   subroutine asm_inc_deallocate_pdaf

     implicit none
     
#if defined key_top
     if (allocated(bgc_bkginc)) deallocate(bgc_bkginc)
#endif

   end subroutine asm_inc_deallocate_pdaf


end module mod_iau_pdaf
