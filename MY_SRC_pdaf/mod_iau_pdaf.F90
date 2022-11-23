!>##Using the Incremental Analysis Method
!>The material in this module is **heavily** based on a similar
!>implementation for the NEMOVAR data assimilation system. Please
!>refer to the `ASM` subdirectory in the `OCE` source code for
!>precise details.
!>
MODULE mod_iau_pdaf

   USE mod_kind_pdaf
#if defined key_top
   USE mod_statevector_pdaf, only: n_trc, sfields, id
#endif
   USE par_oce, &
      ONLY: jpi, jpj, jpk, jpkm1, jp_tem, jp_sal

   IMPLICIT NONE
   SAVE

   !> Array to store the ssh IAU
   REAL(pwp), DIMENSION(:, :), ALLOCATABLE :: ssh_iau_pdaf
   !> Array to store the T, S IAU
   REAL(pwp), DIMENSION(:, :, :), ALLOCATABLE :: t_iau_pdaf, s_iau_pdaf
   !> Array to store the U, V IAU
   REAL(pwp), DIMENSION(:, :, :), ALLOCATABLE :: u_iau_pdaf, v_iau_pdaf
#if defined key_top
   REAL(pwp), DIMENSION(:, :, :, :), ALLOCATABLE :: bgc_iau_pdaf
#endif
   !>
   integer :: niaufn = 0     !: Type of IAU weighing function: = 0   Constant weighting
   !                         !: = 1   Linear hat-like, centred in middle of IAU interval
   integer :: iiauper = 1    ! Length of IAU
   integer :: nitiaustr      ! Timestep of start of IAU interval in [0,nitend-nit000-1]
   integer :: nitiaufin      ! Timestep of end of IAU interval in [0,nitend-nit000-1]
   INTEGER :: nn_divdmp = 0  !: Apply divergence damping filter nn_divdmp times
   real(pwp), ALLOCATABLE :: wgtiau(:) ! IAU weights for each time step

   namelist /iau_nml/ iiauper, niaufn, nn_divdmp

#  include "vectopt_loop_substitute.h90"
contains

   !>This routine initialises the arrays for the IAU. The arrays
   !>**must* be initialised to zero, as PDAF does not compute an
   !>analysis update for all gridpoints. (In this case, no increment
   !>should be added to the model at this gridpoint.)
   !>
   !>The arrays are deallocated in `finalise_pdaf`.
   SUBROUTINE asm_inc_init_pdaf
      integer :: imid
      real(pwp) :: znorm
      integer :: jt
      ! 2D variables
      ALLOCATE (ssh_iau_pdaf(jpi, jpj))
      ssh_iau_pdaf = 0._pwp

      ! 3D variables
      ALLOCATE (t_iau_pdaf(jpi, jpj, jpk), s_iau_pdaf(jpi, jpj, jpk))
      t_iau_pdaf = 0._pwp
      s_iau_pdaf = 0._pwp
      ALLOCATE (u_iau_pdaf(jpi, jpj, jpk), v_iau_pdaf(jpi, jpj, jpk))
      u_iau_pdaf = 0._pwp
      v_iau_pdaf = 0._pwp

#if defined key_top
      ! tracers
      ALLOCATE (bgc_iau_pdaf(jpi, jpj,jpk, n_trc))
#endif

      ! giving the incremental length
      ALLOCATE( wgtiau( iiauper ) )
      wgtiau(:) = 0._pwp
      IF( niaufn == 0 ) THEN           ! Constant IAU forcing
         !                             !---------------------------------------------------------
         DO jt = 1, iiauper
            wgtiau(jt) = 1.0 / REAL( iiauper, pwp )
         END DO
         !                             !---------------------------------------------------------
      ELSEIF ( niaufn == 1 ) THEN      ! Linear hat-like, centred in middle of IAU interval
         !                             !---------------------------------------------------------
         ! Compute the normalization factor
         znorm = 0._pwp
         IF( MOD( iiauper, 2 ) == 0 ) THEN   ! Even number of time steps in IAU interval
            imid = iiauper / 2
            DO jt = 1, imid
               znorm = znorm + REAL( jt )
            END DO
            znorm = 2.0 * znorm
         ELSE                                ! Odd number of time steps in IAU interval
            imid = ( iiauper + 1 ) / 2
            DO jt = 1, imid - 1
               znorm = znorm + REAL( jt )
            END DO
            znorm = 2.0 * znorm + REAL( imid )
         ENDIF
         znorm = 1.0 / znorm
         !
         DO jt = 1, imid - 1
            wgtiau(jt) = REAL( jt ) * znorm
         END DO
         DO jt = imid, iiauper
            wgtiau(jt) = REAL( iiauper - jt + 1 ) * znorm
         END DO
         !
      ENDIF
   END SUBROUTINE asm_inc_init_pdaf

   !>This routine is almost identical to a similar routine
   !>from NEMOVAR. It applies the IAU to the dynamical fields.
   !>
   !>*Called from:* `step.F90`
   SUBROUTINE dyn_asm_inc_pdaf(kt)
      USE oce, &
         ONLY: ua, va
      use dom_oce, only: rdt
      !> Current time step
      INTEGER, INTENT(IN) :: kt

      !> Counter
      INTEGER :: jk
      integer :: it
      real(pwp) :: zincwgt
      ! Check whether to update the dynamic tendencies
      IF (kt <= nitiaufin  .AND. kt >= nitiaustr) THEN
         it = kt - nitiaustr + 1
         zincwgt = wgtiau(it) / rdt
         DO jk = 1, jpkm1
            ua(:, :, jk) = ua(:, :, jk) + u_iau_pdaf(:, :, jk)*zincwgt
            va(:, :, jk) = va(:, :, jk) + v_iau_pdaf(:, :, jk)*zincwgt
         END DO
      END IF
   END SUBROUTINE dyn_asm_inc_pdaf

   !>This routine is almost identical to a similar routine
   !>from NEMOVAR. It applies the IAU to the tracer fields.
   !>
   !>*Called from:* `step.F90`
   SUBROUTINE tra_asm_inc_pdaf(kt)
      USE eosbn2, &
         ONLY: eos_fzp
      USE dom_oce, &
         ONLY: gdept_n
      USE oce, &
         ONLY: tsn, tsa
      use dom_oce, only: rdt

      !> Current time step
      INTEGER, INTENT(IN) :: kt

      !> Counter
      INTEGER  :: jk
      integer :: it
      real(pwp) :: zincwgt
      !> ! 3d freezing point values
      !> Nick: Taken from NEMOVAR. Will this lead to stack overflow?
      REAL(pwp), DIMENSION(jpi, jpj, jpk) :: fzptnz

      real(pwp) :: salfixmin=0.0

      ! Check whether to update the tracer tendencies
      IF (kt <= nitiaufin .AND. kt >= nitiaustr) THEN
         ! Freezing point calculation taken from oc_fz_pt (but calculated for
         ! all depths). Used to prevent the applied increments taking the
         ! temperature below the local freezing point.
         DO jk = 1, jpkm1
            CALL eos_fzp(tsn(:, :, jk, jp_sal), fzptnz(:, :, jk), gdept_n(:, :, jk))
         END DO
         it = kt - nitiaustr + 1
         zincwgt = wgtiau(it) / rdt
         ! Do not apply nonnegative increments.
         ! Do not apply increments if the temperature will fall below freezing
         ! or if the salinity will fall below a specified minimum value.
         DO jk = 1, jpkm1
            WHERE (t_iau_pdaf(:, :, jk) > 0.0_pwp .OR. &
                   tsn(:, :, jk, jp_tem) + tsa(:, :, jk, jp_tem) + t_iau_pdaf(:, :, jk) * wgtiau(it) &
                   > fzptnz(:, :, jk))
               tsa(:, :, jk, jp_tem) = tsa(:, :, jk, jp_tem) + t_iau_pdaf(:, :, jk)*zincwgt
            END WHERE

            WHERE (s_iau_pdaf(:, :, jk) > 0.0_pwp .OR. &
                   tsn(:, :, jk, jp_sal) + tsa(:, :, jk, jp_sal) + s_iau_pdaf(:, :, jk) * wgtiau(it) &
                   > salfixmin)
               tsa(:, :, jk, jp_sal) = tsa(:, :, jk, jp_sal) + s_iau_pdaf(:, :, jk)
            END WHERE

         END DO
      END IF

   END SUBROUTINE tra_asm_inc_pdaf


   !>This routine is almost identical to a similar routine
   !>from NEMOVAR. It applies the IAU to the ssh fields.
   !>
   !>*Called from:* `step.F90`, `ssh_asm_div_pdaf`
   SUBROUTINE ssh_asm_inc_pdaf(kt)
      use dom_oce, only: rdt

      !> Current time step
      INTEGER, INTENT(IN) :: kt

      !> Counter
      integer :: it
      real(pwp) :: zincwgt

      ! Check whether to update the dynamic tendencies
      IF (kt <= nitiaufin  .AND. kt >= nitiaustr) THEN
         it = kt - nitiaustr + 1
         zincwgt = wgtiau(it) / rdt
         ssh_iau_pdaf(:,:) = ssh_iau_pdaf(:,:) * zincwgt
      ELSE IF( kt == nitiaufin+1 ) THEN
         ssh_iau_pdaf = 0._pwp
      END IF

   END SUBROUTINE ssh_asm_inc_pdaf

   !>This routine is almost identical to a similar routine
   !>from NEMOVAR. It applies the IAU to the ssh divergence term.
   !>
   !>*Called from:* `divhor.F90`
   SUBROUTINE ssh_asm_div_pdaf(kt, phdivn)
      USE mod_parallel_pdaf, &
         ONLY: abort_parallel
      USE dom_oce, &
         ONLY: e3t_n, ht_n, ln_linssh, tmask, ssmask

      !> Current time step
      INTEGER, INTENT(IN) :: kt
      !> Horizontal divergence
      REAL(pwp), DIMENSION(:, :, :), INTENT(inout) :: phdivn

      integer :: jk
      real(pwp), dimension(:, :), ALLOCATABLE :: ztim ! local array

      ! Check whether to update the tracer tendencies
      CALL ssh_asm_inc_pdaf( kt )
      ! NEMO-PDAF currently only implemented for linear free surface.
      IF (ln_linssh) THEN
         phdivn(:, :, 1) = phdivn(:, :, 1) - &
                           ssh_iau_pdaf(:, :)/e3t_n(:, :, 1)*tmask(:, :, 1)
      ELSE
         ALLOCATE( ztim(jpi,jpj) )
         ztim(:,:) = ssh_iau_pdaf(:,:) / ( ht_n(:,:) + 1.0 - ssmask(:,:) )
         DO jk = 1, jpkm1
            phdivn(:,:,jk) = phdivn(:,:,jk) - ztim(:,:) * tmask(:,:,jk)
         END DO
         !
         DEALLOCATE(ztim)
      END IF
   END SUBROUTINE ssh_asm_div_pdaf

   SUBROUTINE init_IAU_pdaf(kt)
      integer, intent(in) :: kt
      nitiaustr = kt
      nitiaufin = kt + iiauper - 1
   END SUBROUTINE init_IAU_pdaf

   !>This routine is almost identical to a similar routine
   !>from NEMOVAR. It applies the divergence damping to the velocity.
   !>
   subroutine div_damping_filter
      use par_oce, only: jpkm1, jpjm1, jpim1
      use dom_oce, only: e1v, e3v_n, e2u, e3u_n, e3t_n, &
                         r1_e1u, r1_e2v, umask, vmask
      use lbclnk, only: lbc_lnk

      integer :: jt, jk, jj, ji  ! counter
      REAL(pwp), ALLOCATABLE, DIMENSION(:,:) ::   zhdiv   ! 2D workspace

      ALLOCATE( zhdiv(jpi,jpj) )
      !
      DO jt = 1, nn_divdmp
         !
         DO jk = 1, jpkm1           ! zhdiv = e1e1 * div
            zhdiv(:,:) = 0._pwp
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zhdiv(ji,jj) = (  e2u(ji  ,jj) * e3u_n(ji  ,jj,jk) * u_iau_pdaf(ji  ,jj,jk)    &
                     &            - e2u(ji-1,jj) * e3u_n(ji-1,jj,jk) * u_iau_pdaf(ji-1,jj,jk)    &
                     &            + e1v(ji,jj  ) * e3v_n(ji,jj  ,jk) * v_iau_pdaf(ji,jj  ,jk)    &
                     &            - e1v(ji,jj-1) * e3v_n(ji,jj-1,jk) * v_iau_pdaf(ji,jj-1,jk)  ) / e3t_n(ji,jj,jk)
               END DO
            END DO
            CALL lbc_lnk( 'mod_iau_pdaf', zhdiv, 'T', 1._pwp)   ! lateral boundary cond. (no sign change)
            !
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  u_iau_pdaf(ji,jj,jk) = u_iau_pdaf(ji,jj,jk) &
                     &               + 0.2_pwp * ( zhdiv(ji+1,jj) - zhdiv(ji  ,jj) ) * r1_e1u(ji,jj) * umask(ji,jj,jk)
                  v_iau_pdaf(ji,jj,jk) = v_iau_pdaf(ji,jj,jk) &
                     &               + 0.2_pwp * ( zhdiv(ji,jj+1) - zhdiv(ji,jj  ) ) * r1_e2v(ji,jj) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !
      END DO
      !
      DEALLOCATE( zhdiv )
   end subroutine div_damping_filter


   SUBROUTINE bgc3d_asm_inc( kt )
      use dom_oce, only: rdt
      use lib_mpp, only: ctl_stop
      USE trc, ONLY:           & ! passive tracer variables
      & trn,                &
      & trb      
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dyn_asm_inc  ***
      !!
      !! ** Purpose : Apply generic 3D biogeochemistry assimilation increments.
      !!
      !! ** Action  :
      !!----------------------------------------------------------------------
      INTEGER,  INTENT(IN) :: kt      ! Current time step
      !
      INTEGER   :: it              ! Index
      integer   :: jptrc
      integer   :: j, id_var
      REAL(pwp) :: zincwgt         ! IAU weight for current time step
      !!----------------------------------------------------------------------

      !--------------------------------------------------------------------
      ! Incremental Analysis Updating
      !--------------------------------------------------------------------

      IF (kt <= nitiaufin  .AND. kt >= nitiaustr) THEN

         it = kt - nitiaustr + 1
         zincwgt = wgtiau(it)   ! IAU weight for the current time step
         ! note this is not a tendency so should not be divided by rdt

         ! Update the 3D BGC variables
         ! Add directly to trn and trb, rather than to tra, because tra gets
         ! reset to zero at the start of trc_stp, called after this routine
         ! Don't apply increments if they'll take concentrations negative

#if defined key_top
         do j = 1, n_trc
            id_var = id%trcs(j)
            jptrc = sfields(id_var)%jptrc
            WHERE( bgc_iau_pdaf(:,:,:,j) > 0.0_pwp .OR. &
                 & trn(:,:,:,jptrc) + bgc_iau_pdaf(:,:,:,j) * zincwgt > 0.0_pwp )
               trn(:,:,:,jptrc) = trn(:,:,:,jptrc) + bgc_iau_pdaf(:,:,:,j) * zincwgt
               trb(:,:,:,jptrc) = trb(:,:,:,jptrc) + bgc_iau_pdaf(:,:,:,j) * zincwgt
            END WHERE
         end do
#else
         CALL ctl_stop ( ' bgc3d_asm_inc: no compatible BGC model defined' )
#endif

      ENDIF
      !
   END SUBROUTINE bgc3d_asm_inc
END MODULE mod_iau_pdaf
