MODULE trcupdate_ergom
#if defined key_top
   !!======================================================================
   !!                         ***  MODULE trcsms_ergom  ***
   !! TOP :   Main module of the ERGOM tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !! trc_sms_ergom       : ERGOM model main routine
   !! trc_sms_ergom_alloc : allocate arrays specific to ERGOM sms
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc
   USE trcbc!, only : trc_bc
   
   USE bioparam
   USE biolight
   USE trdtrc
   USE biomod 
   USE sms_ergom
   USE iom             ! incase   CALL iom_put( "VARA" , vara )
   USE zdf_oce         ! diffusivity
   USE lbclnk          ! 
!  USE sbcmod, only: ln_blk_core
   USE sbcrnf, ONLY: h_rnf   ! Height the runoff is distributed over.
   USE sbcrnf, ONLY: nk_rnf  ! Number of levels the runoff is distributed over.
   USE sbcblk
   USE sbcdcy         ! surface boundary condition: diurnal cycle
   !USE sbc_oce         ! adopt atmospheric forcing from SBC routines
   !USE sbc_blk    !only : atm_bio_u10, atm_bio_v10, atm_bio_swr, atm_bio_hum
   USE ice,       only : at_i!,icethi
   !USE ldftra_oce
   USE biosun,    only : shortwaverad 
   USE constants, only : grarad !,zeitstu
   USE dom_oce,   only : nsec_year, nmonth, nyear
   USE phycst,    only : r1_rau0

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ergom_update       ! called by trctrp.F90 module

   ! Defined HERE the arrays specific to ERGOM sms and ALLOCATE them 
   ! in trc_sms_ergom_alloc

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_ergom.F90 5075 2015-02-11 10:50:34Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ergom_update( kt )
      !!----------------------------------------------------------------------
      !!                     ***  trc_ergom_update  ***
      !!
      !! ** Purpose :   Managment of the call to Biological  
      !!               routines of ERGOM bio-model ( update step of ERGOM) 
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER ::   i, jn, ji, jj, jk, jpkoce,jpko, ikt, bkt, kh,jjm1,jim1,ikbot        ! dummy loop index
      !INTEGER :: zeitstu
      INTEGER(4) ::   kmx, iw2l, ndtp , rivmap     ! dummy indexes for hbm
      INTEGER(4),SAVE, DIMENSION(:), ALLOCATABLE :: kh2 ,mcol
      INTEGER(4),SAVE, DIMENSION(:,:), ALLOCATABLE :: ind , idx,msrf 
      REAL(wp):: xkoor, ykoor, dtbio     !ERGOM VARS
      REAL(wp):: ubt, vbt         !bottom velocities
      REAL(wp):: npr,cloud,humid,wu,wv        !atmospheric forcing
      REAL(wp) ::   za1, zb, zt, zsc, zv, zo2s, zdep ! river oxygen
      LOGICAL ::  casus
      REAL(wp),SAVE, DIMENSION(:), ALLOCATABLE :: h ,t,s,light,dispv,netpp,kpar
      REAL(wp),SAVE, DIMENSION(:), ALLOCATABLE :: ph, pco2, rho
      !REAL(wp),SAVE, DIMENSION(:), ALLOCATABLE :: nitrific, denitrific, uptake !3D diagn output
      !REAL(wp), SAVE :: ben_nitr !2D diagnostics output of benthic
      REAL(wp),SAVE, DIMENSION(:,:), ALLOCATABLE :: ts,benthos, components
      REAL(wp),SAVE, DIMENSION(:,:), ALLOCATABLE :: srfflx
      REAL(wp), DIMENSION(jpi,jpj)     ::   z2d   ! 2D workspace
      REAL(wp)                         :: r_chl_to_n, frac  ! for CHL calculation
      INTEGER                          :: jtk, len_kpar 
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('trc_ergom_update')
      ! 
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' trc_ergom_update:  ERGOM model update'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
         
         ntbio=jp_bgc
         ALLOCATE( t(0:jpk),s(0:jpk) , ts(0:jpk,1:ntbio), benthos(0:jpk,1:ntben))
         ALLOCATE( light(0:jpk), ind(2,1), idx(1,jpk), components(0:jpk,ntbio))
         ALLOCATE( kh2(jpk), mcol(2), srfflx(0:jpk,1:ntsrfflx), dispv(0:jpk))
         ALLOCATE( netpp(0:jpk), kpar(0:jpk), msrf(1,1), h(0:jpk))
         ALLOCATE( rho(0:jpk), ph(0:jpk), pco2(0:jpk))    
         !ALLOCATE( nitrific(0:jpk), denitrific(0:jpk), uptake(0:jpk)) 
         mcol=0
      ENDIF
      
      IF(lwp) WRITE(numout,*) ' trc_ergom_update:  ERGOM model kt',kt 
      
      ! Initialize as zeros
      xlight(1:jpi,1:jpj,1:jpk) = 0._wp
      xnetpp(1:jpi,1:jpj,1:jpk) = 0._wp
      xsecchi(1:jpi,1:jpj,1:jpk)= 0._wp
      xph(1:jpi,1:jpj,1:jpk)    = 0._wp
      xpco2(1:jpi,1:jpj,1:jpk)  = 0._wp
      xchl(1:jpi,1:jpj,1:jpk)   = 0._wp
      xphy(1:jpi,1:jpj,1:jpk)   = 0._wp
      xzoo(1:jpi,1:jpj,1:jpk)   = 0._wp
      xkpar(1:jpi,1:jpj,1:jpk)  = 0._wp
   
      IF( kt == nit000 ) THEN  ! initialize (if statement can be deleted for monthly changing of atm_pco2)            
           CALL init_atm_pco2_from_month(nmonth, nyear)
      ENDIF

      
! TODO:  wrap biomod.f90
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( tmask(ji,jj,1) == 1 ) THEN 
               ikt = mikt(ji,jj)
               bkt = mbkt(ji,jj)
               ! SURFACE VAR
               xkoor =  glamt(ji,jj) * grarad
               ykoor =  gphit(ji,jj) * grarad
               kh    =  mbkt(ji,jj)  
               !COLUMN VAR
               h(0) = 0.
               h(1:) = e3t_n(ji,jj,:)
               ! bottom velocity as calculated as in zdfgls.F90:  
               !  maybe we could use    ustarb2 !: Squared bottom  velocity scale at T-points
               jim1=max(1,ji-1) 
               jjm1=max(1,jj-1)
               ubt= ( un(ji,jj,mbku(ji,jj)) + un(jim1,jj,mbku(jim1,jj))  ) * ( 1._wp - 0.5_wp * umask(ji,jj,1) * umask(jim1,jj,1)  )
               vbt= ( vn(ji,jj,mbkv(ji,jj)) + vn(ji,jjm1,mbkv(ji,jjm1))  ) * ( 1._wp - 0.5_wp * vmask(ji,jj,1) * vmask(ji,jjm1,1)  )    
                !! TODO cases if drying occurs
                !! TODO rewrite CALL bio_dynamics( ...  , such that  cloud and humid > atm_forcing_array
                !!  NB:     In ergom biomod  if cloud == -1  then swr =humid  else  swr is calculated from cloud and humid
                 cloud =  -1 
		 humid =  max(atm_bio_swr(ji,jj),0._wp) 
      
                wu = atm_bio_u10(ji,jj) 
                wv = atm_bio_v10(ji,jj) 
                iw2l=1
                mcol=1+1
                idx=1
                kh2(:)=bkt
                IF (at_i  (ji,jj) .gt. 0.5) THEN   ! ice total fractional area (ice concentration) 
                   casus=.false.  ! asume ice 
                ELSE
                   casus=.true.   ! asume no ice 
                ENDIF

                ind(1,1)=ji
                ind(2,1)=jj
                kmx=jpk  
                
		!light
		light(0) = 0._wp
		light(1:)= 0._wp
		kpar(0)  = 0._wp
		kpar(1:) = 0._wp
                
		! temperature, salinity, denstiy
                t(0)  = 0._wp
                t(1:) = tsn(ji,jj,:,jp_tem)
                s(0)  = 0._wp
                s(1:) = tsn(ji,jj,:,jp_sal)
                rho(0) = 0._wp
                rho(1:) = rhop(ji,jj,:) !potential volumic mass (kg m-3)
                
		! Pellagic biovariable,
                ts(0,1:)  = 0._wp
                ts(1:,1:) = tra(ji,jj,:,:) 
                components= ts
                dispv(0)  = 0._wp
                dispv(1:) = avt(ji,jj,:)
                dispv =max(dispv,1.e-12)
               !  Should do same as call bio_limit
                DO i=1 , ntbio          
                   IF (i .ne. idx_oxy)  THEN
                      components(ikt:bkt,i) =max(ts(ikt:bkt,i),0.) /Nnorm
                   ELSE
                      components(ikt:bkt,i) =ts(ikt:bkt,i) /Onorm
                   ENDIF
                ENDDO
                ! Sediments
                benthos(1,:) = seda(ji,jj,:) /Nnorm 
                IF( kt == nit000 ) THEN  ! initialize             
               !    CALL init_atm_pco2(nmonth,zeitstu, nyear)
		   CALL init_biolight(iw2l,mcol,kh2,ind,kmx,ntbio,idx,xkoor, &
                                      ykoor,h,casus,cloud,humid,components,  &
                                      light,kpar)                 
                ENDIF

                msrf=1
                npr=np_ratio(ji,jj)
		!COUPLING TIMESTEP
                IF (ln_top_euler) THEN !for explicit Euler BGC transport 1Dt
                        dtbio=rn_rdt
                ELSE
                        dtbio=2._wp*rn_rdt !for leap-frog BGC transport 2Dt
                ENDIF
                ndtp=1
                rivmap=1
                srfflx=0.
 
                CALL bio_dynamics(iw2l, kmx, msrf, mcol, ind, kh2, h,ubt, vbt, wu, wv,&
                                 t, s, rho, humid, cloud, components, benthos,       &
                                 xkoor, ykoor, dispv, casus, npr,   &
                                 dtbio, idx, ndtp, rivmap, light,   &
                                 srfflx, netpp, kpar, ph, pco2 )
                                 !srfflx, netpp, kpar, ph, pco2, nitrific, denitrific, &
                                 !uptake, ben_nitr)
               
                DO i=1 , ntbio          
                   IF (i .ne. idx_oxy)  THEN
                      ts(ikt:bkt,i) =max(components(ikt:bkt,i),0.) *Nnorm
                   ELSE
                      ts(ikt:bkt,i)=components(ikt:bkt,i) *Onorm
                   ENDIF
                ENDDO

                tra(ji,jj,:,:)           = ts(1:,1:)   
                seda(ji,jj,:)            = benthos(1,:) * Nnorm 
                xlight(ji,jj,:)          = light(1:)
                xnetpp(ji,jj,:)          = netpp(1:) * Nnorm * rfc2n * mweig_C

                len_kpar = SIZE(kpar(1:))

                DO jtk = 1, len_kpar
                   
                   IF (kpar(jtk) /= 0.0) THEN

                      xsecchi(ji,jj,jtk)         = secchi_coeff * kpar(jtk)**secchi_exp    

                   ELSE

                      xsecchi(ji,jj,jtk) = 0.010_wp
                      
                   END IF

                END DO



         xph(ji,jj,:)             = ph(1:)
         xpco2(ji,jj,:)           = pco2(1:)
         xkpar(ji,jj,:)		 = kpar(1:)
		
		
                !IF (jp_dia3d > 0) trc3d(ji,jj,:,1)         = nitrific(1:)
                !IF (jp_dia3d > 1) trc3d(ji,jj,:,2)         = uptake(1:)
                !IF (jp_dia3d > 2) trc3d(ji,jj,:,3)         = denitrific(1:)
                !IF (jp_dia2d > 0) trc2d(ji,jj,1) = ben_nitr
	 	
		
            ENDIF   ! tmask
         ENDDO  ! ji
      ENDDO  ! jj
      
      ! 
      IF ( iom_use ( "PP" ))   CALL iom_put( "PP"    , xnetpp(:,:,:)  )
      IF ( iom_use ("SECCHI")) CALL iom_put( "SECCHI", xsecchi(:,:,:) )
      IF ( iom_use ("PH"))     CALL iom_put( "PH"    , xph(:,:,:)     )
      IF ( iom_use ("PCO2"))   CALL iom_put( "PCO2"  , xpco2(:,:,:)   )
      IF ( iom_use ("LIGHT"))  CALL iom_put( "LIGHT" , xlight(:,:,:)  ) 
      IF ( iom_use ("KD"))     CALL iom_put( "KD"    , xkpar(:,:,:)  ) 

      IF ( iom_use ( "sbo" )) then
	 DO jj = 1, jpj
            DO ji = 1, jpi
               ikbot = mbkt(ji,jj)
               z2d(ji,jj) = trn(ji,jj,ikbot,idx_oxy)
            END DO
         END DO
         CALL iom_put( "sbo", z2d )                ! bottom oxygen
      ENDIF

      !! IRI moved from trcwri_ergom on 19/10-2021
      ! write the diagnostic concentrations in the file
      ! moved together with rest of ERGOM update step to
      ! trcupdate_ergom on 15/03-2021
      ! ---------------------------------------
      flush(numout)

      ! For PDAF: Always compute xchl for the data assimilation
      frac = max_chl_to_n / (min_chl_to_n*light_max)
      DO jk = 1,jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tmask(ji,jj,1) == 1) THEN
                  r_chl_to_n = max( min_chl_to_n, max_chl_to_n*(one-frac*xlight(ji,jj,jk)) )
                  xchl(ji,jj,jk) = (trn(ji,jj,jk,idx_dia) + trn(ji,jj,jk,idx_flag) + &
                       & trn(ji,jj,jk,idx_cyano)) * r_chl_to_n
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      IF ( iom_use ( "CHL" )) THEN
         CALL iom_put( "CHL",    xchl(:,:,:) )
      ENDIF
      !!
      
      IF ( iom_use ( "PHYC" )) THEN
         DO jk = 1,jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( tmask(ji,jj,1) == 1) THEN
                     xphy(ji,jj,jk) = (trn(ji,jj,jk,idx_dia) + trn(ji,jj,jk,idx_flag) + &
                     & trn(ji,jj,jk,idx_cyano)) * rfc2n
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         CALL iom_put( "PHYC",    xphy(:,:,:) )
      ENDIF
      
      IF ( iom_use ( "ZOOC" )) THEN
         DO jk = 1,jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( tmask(ji,jj,1) == 1) THEN
                     xzoo(ji,jj,jk) = (trn(ji,jj,jk,idx_miz) + trn(ji,jj,jk,idx_mez)) * rfc2n
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         CALL iom_put( "ZOOC",    xzoo(:,:,:) )
      ENDIF

!Todo: Implement the indexing jn
      ! Save the trends in the mixed layer
!      IF( l_trdtrc ) THEN
!          DO jn = jp_erg0, jp_erg1
!            ztrmyt(:,:,:) = tra(:,:,:,jn)
!            CALL trd_trc( ztrmyt, jn, jptra_sms, kt )   ! save trends
!          END DO
!          DEALLOCATE( ztrmyt )
!      END IF
                
      !
      IF( ln_timing )  CALL timing_stop('trc_ergom_update')
      !
   END SUBROUTINE trc_ergom_update
#endif
END MODULE trcupdate_ergom
