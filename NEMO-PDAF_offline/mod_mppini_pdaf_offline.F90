!>Model domain decomposition for model parallelization
!!
!! The routine is called in the initialize subroutine
!! for offline pdaf such that the offline version 
!! can still have a model decomposition.
!! The code is heavily based on the code in NEMO in the
!! src/OCE/LBC/mppini.F90
!!
module mod_mppini_pdaf_offline
use mod_kind_pdaf, only: pwp
use mod_parallel_pdaf, only: abort_parallel
implicit none
contains
   subroutine get_decomposition(mype, jpiglo, jpjglo, knbij, nlei, nlej, nimpp, njmpp)
      ! arguments
      INTEGER, INTENT(in)  ::  mype            ! rank of current model (mype_model)       
      INTEGER, INTENT(in)  ::  knbij           ! total number of subdomains (npes_model)
      INTEGER, INTENT(in)  ::  jpiglo, jpjglo  ! global number gridpoints of the domain
      INTEGER, intent(out) ::  nlei, nlej      ! local number of grid points in i, j direction
      INTEGER, intent(out) ::  nimpp, njmpp    ! local starting index in global domain

      ! local variables
      INTEGER              ::  knbi, knbj                 ! number of subdomains along i and j
      INTEGER              ::  iin(knbij), ijn(knbij)     ! conversion between npes and i, j index
      INTEGER, allocatable ::  klci(:, :), klcj(:, :)     ! equivalent to nlci, nlcj in nemo
      INTEGER, allocatable ::  kimppt(:, :), kjmppt(:, :) ! nimpp and njmpp

      ! domain partition
      call mpp_init_bestpartition( jpiglo, jpjglo, knbij, knbi, knbj )

      ! get local domain information
      allocate(klci(knbi, knbj), klcj(knbi, knbj))
      allocate(kimppt(knbi, knbj), kjmppt(knbi, knbj))
      call mpp_basic_decomposition( jpiglo, jpjglo, knbi, knbj, kimppt, kjmppt, klci, klcj)

      ! get conversion array
      call get_proc_index(knbi, knbj, iin, ijn)

      ! assign to nemo values
      nlei = klci(iin(mype+1), ijn(mype+1))
      nlej = klcj(iin(mype+1), ijn(mype+1))
      nimpp = kimppt(iin(mype+1), ijn(mype+1))
      njmpp = kjmppt(iin(mype+1), ijn(mype+1))
      deallocate(kimppt, kjmppt, klci, klcj)
   end subroutine get_decomposition

   SUBROUTINE mpp_init_bestpartition( jpiglo, jpjglo, knbij, knbi, knbj )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE mpp_init_bestpartition  ***
      !!
      !! ** Purpose :
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)    ::   jpiglo
      INTEGER, INTENT(in)    ::   jpjglo
      INTEGER, INTENT(in   ) ::   knbij         ! total number if subdomains (knbi*knbj)
      INTEGER, INTENT(  out) ::   knbi, knbj    ! number if subdomains along i and j (knbi and knbj)
      !
      INTEGER :: ji, jj, ii
      INTEGER :: iszitst, iszjtst
      INTEGER :: isziref, iszjref
      INTEGER :: inbij
      INTEGER :: inbimax, inbjmax, inbijmax
      INTEGER :: isz1
      INTEGER, DIMENSION(  :), ALLOCATABLE :: inbi0, inbj0   ! number of subdomains along i,j
      INTEGER, DIMENSION(  :), ALLOCATABLE :: iszi0, iszj0   ! max size of the subdomains along i,j
      INTEGER, DIMENSION(  :), ALLOCATABLE :: inbi1, inbj1, inbij1   ! number of subdomains along i,j
      INTEGER, DIMENSION(  :), ALLOCATABLE :: iszi1, iszj1, iszij1   ! max size of the subdomains along i,j
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: llmsk2d                 ! max size of the subdomains along i,j

      !!----------------------------------------------------------------------
      !


      inbij = knbij

      inbijmax = inbij

      !
      ALLOCATE(inbi0(inbijmax),inbj0(inbijmax),iszi0(inbijmax),iszj0(inbijmax))
      !
      inbimax = 0
      inbjmax = 0
      isziref = jpiglo*jpjglo+1
      iszjref = jpiglo*jpjglo+1
      !
      ! get the list of knbi that gives a smaller jpimax than knbi-1
      ! get the list of knbj that gives a smaller jpjmax than knbj-1
      DO ji = 1, inbijmax
         iszitst = ( jpiglo + (ji-1) ) / ji
         IF( iszitst < isziref ) THEN
            isziref = iszitst
            inbimax = inbimax + 1
            inbi0(inbimax) = ji
            iszi0(inbimax) = isziref
         ENDIF
         iszjtst = ( jpjglo + (ji-1) ) / ji
         IF( iszjtst < iszjref ) THEN
            iszjref = iszjtst
            inbjmax = inbjmax + 1
            inbj0(inbjmax) = ji
            iszj0(inbjmax) = iszjref
         ENDIF
      END DO

      ! combine these 2 lists to get all possible knbi*knbj <  inbijmax
      ALLOCATE( llmsk2d(inbimax,inbjmax) )
      DO jj = 1, inbjmax
         DO ji = 1, inbimax
            IF ( inbi0(ji) * inbj0(jj) <= inbijmax ) THEN   ;   llmsk2d(ji,jj) = .TRUE.
            ELSE                                            ;   llmsk2d(ji,jj) = .FALSE.
            ENDIF
         END DO
      END DO
      isz1 = COUNT(llmsk2d)
      ALLOCATE( inbi1(isz1), inbj1(isz1), iszi1(isz1), iszj1(isz1) )
      ii = 0
      DO jj = 1, inbjmax
         DO ji = 1, inbimax
            IF( llmsk2d(ji,jj) .EQV. .TRUE. ) THEN
               ii = ii + 1
               inbi1(ii) = inbi0(ji)
               inbj1(ii) = inbj0(jj)
               iszi1(ii) = iszi0(ji)
               iszj1(ii) = iszj0(jj)
            END IF
         END DO
      END DO
      DEALLOCATE( inbi0, inbj0, iszi0, iszj0 )
      DEALLOCATE( llmsk2d )

      ALLOCATE( inbij1(isz1), iszij1(isz1) )
      inbij1(:) = inbi1(:) * inbj1(:)
      iszij1(:) = iszi1(:) * iszj1(:)

      ! if therr is no land and no print
      ! get the smaller partition which gives the smallest subdomain size
      ii = MINLOC(inbij1, mask = iszij1 == MINVAL(iszij1), dim = 1)
      knbi = inbi1(ii)
      knbj = inbj1(ii)
      DEALLOCATE( inbi1, inbj1, inbij1, iszi1, iszj1, iszij1 )
   END SUBROUTINE mpp_init_bestpartition

   SUBROUTINE mpp_basic_decomposition( jpiglo, jpjglo, knbi, knbj, kimppt, kjmppt, klci, klcj)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE mpp_basic_decomposition  ***
      !!                    
      !! ** Purpose :   Lay out the global domain over processors.
      !!
      !! ** Method  :   Global domain is distributed in smaller local domains.
      !!
      !! ** Action : - set for all knbi*knbj domains:
      !!                    kimppt     : longitudinal index
      !!                    kjmppt     : latitudinal  index
      !!                    klci       : first dimension
      !!                    klcj       : second dimension
      !!----------------------------------------------------------------------
      INTEGER,                       INTENT(in   ) ::   jpiglo, jpjglo
      INTEGER,                       INTENT(in   ) ::   knbi, knbj
      INTEGER, DIMENSION(knbi,knbj), INTENT(  out) ::   kimppt, kjmppt
      INTEGER, DIMENSION(knbi,knbj), INTENT(  out) ::   klci, klcj
      !
      integer ::   kimax, kjmax
      INTEGER ::   iresti, irestj, ijpjmin
      integer :: ji, jj

      !!----------------------------------------------------------------------
      !
      kimax = ( jpiglo + (knbi-1) ) / knbi    ! first  dim.
      kjmax = ( jpjglo + (knbj-1) ) / knbj    ! second dim.
      !
      !  1. Dimension arrays for subdomains
      ! -----------------------------------
      !  Computation of local domain sizes klci() klcj()
      !  These dimensions depend on global sizes knbi,knbj and jpiglo,jpjglo
      !  The subdomains are squares lesser than or equal to the global
      !  dimensions divided by the number of processors minus the overlap array.
      !
      iresti = 1 + MOD( jpiglo -1 , knbi )
      irestj = 1 + MOD( jpjglo -1 , knbj )
      !
      !  Need to use kimax and kjmax here since jpi and jpj not yet defined
      klci(1:iresti      ,:) = kimax
      klci(iresti+1:knbi ,:) = kimax-1
      IF( MINVAL(klci) < 3 ) THEN
         WRITE(*,*) '   mpp_basic_decomposition: minimum value of jpi must be >= 3'
         WRITE(*,*) '   We have ', MINVAL(klci)
      ENDIF

      ijpjmin = 3
      klcj(:,      1:irestj) = kjmax
      klcj(:, irestj+1:knbj) = kjmax-1
      IF( MINVAL(klcj) < ijpjmin ) THEN
         WRITE(*,*) '   mpp_basic_decomposition: minimum value of jpj must be >= ', ijpjmin
         WRITE(*,*) '   We have ', MINVAL(klcj)
         CALL abort_parallel()
      ENDIF

      !  2. Index arrays for subdomains
      ! -------------------------------
      kimppt(:,:) = 1
      kjmppt(:,:) = 1
      !
      IF( knbi > 1 ) THEN
         DO jj = 1, knbj
            DO ji = 2, knbi
               kimppt(ji,jj) = kimppt(ji-1,jj) + klci(ji-1,jj)
            END DO
         END DO
      ENDIF
      !
      IF( knbj > 1 )THEN
         DO jj = 2, knbj
            DO ji = 1, knbi
               kjmppt(ji,jj) = kjmppt(ji,jj-1) + klcj(ji,jj-1)
            END DO
         END DO
      ENDIF
   END SUBROUTINE mpp_basic_decomposition

   subroutine get_proc_index(jpni, jpnj, iin, ijn)
      integer, intent(in) :: jpni, jpnj
      integer, intent(out) :: iin(jpni*jpnj)
      integer, intent(out) :: ijn(jpni*jpnj)
      integer :: jarea
      integer :: iarea0
      integer :: ii, ij
      integer :: icont

      icont = -1
      DO jarea = 1, jpni*jpnj
         iarea0 = jarea - 1
         ii = 1 + MOD(iarea0,jpni)
         ij = 1 +     iarea0/jpni
         icont = icont + 1
         iin(icont+1) = ii
         ijn(icont+1) = ij
      END DO
   end subroutine get_proc_index
end module mod_mppini_pdaf_offline