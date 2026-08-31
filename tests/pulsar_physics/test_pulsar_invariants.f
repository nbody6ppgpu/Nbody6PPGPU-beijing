C     Standalone invariant checks for the PULSARINIT/PULSAREVO physics
C     in src/Main/pulsar.F, linked directly against the built
C     pulsar.o/ran2.o (no lifecycle hook calls these yet, so this is
C     the only real verification available before Stage 4 wires them
C     up). Exercises the plan's section 10.4 invariants:
C       - birth distributions: finite, positive P/B/Pdot for every
C         supported PSR_PMODE/PSR_BMODE combination (normal and AIC)
C       - isolated spin-down: B decays monotonically toward PSR_BMIN
C         and never below it; P and Pdot stay finite and positive
C       - stable RLOF accretion: no negative mass/period/field, no NaN
C     Fixed-form Fortran so common6.h/pulsar_ctrl.h INCLUDE safely.
      PROGRAM TEST_PULSAR_INVARIANTS
C     No IMPLICIT NONE here: common6.h supplies IMPLICIT REAL*8
C     (A-H,O-Z) which the whole codebase (and this file's own I-N
C     range names below) relies on; every local variable used here is
C     given an explicit type regardless.
      INCLUDE 'common6.h'
      INCLUDE 'pulsar_ctrl.h'

      INTEGER NFAIL, PMODE, BI, BMODE, BMODES(4)
      REAL*8 PULSARINITB, PULSARINITP, PULSARINITPD, PSRMASS
      REAL*8 B_I, B_F, PERIOD_I, PERIOD_F, PDOT_I, PDOT_F
      REAL*8 MASS, MASS_COMP, RAD_COMP, PSRMDOT, DTM, NSINCL
      REAL*8 B_PREV, DTM_STEP, BMIN_TESLA
      REAL*8 B_ISO, P_ISO, PD_ISO
      INTEGER ISTEP
      LOGICAL CE

      DATA BMODES /1, 2, 4, 5/

      NFAIL = 0
      IDUM1 = -12345

      PSR_SPINA = 0.01D0
      PSR_SPINB = 0.5D0
      PSR_BMAGA = 12.0D0
      PSR_BMAGB = 13.0D0
      PSRM_SPINA = 0.01D0
      PSRM_SPINB = 0.5D0
      PSRM_BMAGA = 12.0D0
      PSRM_BMAGB = 13.0D0
      PSR_BMIN = 8.0D0
      PSR_RAD = 10.0D0
      PSR_TAU = 1.0D0
      PSR_KAPPA = 1.0D-3
      PSR_EP = 1.0D0
      PSR_ACC_CE = 0
      PSRMASS = 1.4D0

C     --- Birth invariants: every supported PMODE x BMODE combo ---
      DO 10 PMODE = 1, 4
         DO 20 BI = 1, 4
            BMODE = BMODES(BI)
            PSR_PMODE = PMODE
            PSR_BMODE = BMODE
            CALL PULSARINIT(PULSARINITB, PULSARINITP, PULSARINITPD,
     &           PSRMASS, .FALSE.)
            CALL CHKPOS(1, PMODE, BMODE, PULSARINITB, NFAIL)
            CALL CHKPOS(2, PMODE, BMODE, PULSARINITP, NFAIL)
            CALL CHKPOS(3, PMODE, BMODE, PULSARINITPD, NFAIL)
   20    CONTINUE
   10 CONTINUE

C     Same, for the AIC (PSRM_*) birth path.
      DO 30 PMODE = 1, 4
         DO 40 BI = 1, 4
            BMODE = BMODES(BI)
            PSRM_PMODE = PMODE
            PSRM_BMODE = BMODE
            CALL PULSARINIT(PULSARINITB, PULSARINITP, PULSARINITPD,
     &           PSRMASS, .TRUE.)
            CALL CHKPOS(4, PMODE, BMODE, PULSARINITB, NFAIL)
            CALL CHKPOS(5, PMODE, BMODE, PULSARINITP, NFAIL)
            CALL CHKPOS(6, PMODE, BMODE, PULSARINITPD, NFAIL)
   40    CONTINUE
   30 CONTINUE

C     --- Isolated spin-down: B -> PSR_BMIN (never below), P/Pdot ok ---
      PSR_PMODE = 1
      PSR_BMODE = 1
      PSR_SPINA = 0.1D0
      PSR_BMAGA = 13.0D0
      CALL PULSARINIT(PULSARINITB, PULSARINITP, PULSARINITPD,
     &     PSRMASS, .FALSE.)

      B_I = PULSARINITB
      PERIOD_I = PULSARINITP
      PDOT_I = PULSARINITPD
      MASS = PSRMASS
      MASS_COMP = 0.0D0
      RAD_COMP = 0.0D0
      PSRMDOT = 0.0D0
      NSINCL = 1.5708D0
      CE = .FALSE.
      DTM_STEP = 1.0D6 * 3.1536D13
      B_PREV = B_I
      BMIN_TESLA = 10.0D0**PSR_BMIN * 1.0D-4

      DO 50 ISTEP = 1, 50
         CALL PULSAREVO(MASS, MASS_COMP, RAD_COMP, B_I, PERIOD_I,
     &        PDOT_I, PSRMDOT, DTM_STEP, NSINCL, CE, 1,
     &        B_F, PERIOD_F, PDOT_F)
         CALL CHKPOS(7, ISTEP, 0, B_F, NFAIL)
         CALL CHKPOS(8, ISTEP, 0, PERIOD_F, NFAIL)
         CALL CHKPOS(9, ISTEP, 0, PDOT_F, NFAIL)
         IF (B_F.GT.B_PREV+1.0D-30) THEN
            WRITE(6,*) 'FAIL: isolated B increased at step', ISTEP,
     &           B_PREV, B_F
            NFAIL = NFAIL + 1
         ENDIF
         IF (B_F.LT.BMIN_TESLA-1.0D-30) THEN
            WRITE(6,*) 'FAIL: isolated B below PSR_BMIN at step',
     &           ISTEP, B_F, BMIN_TESLA
            NFAIL = NFAIL + 1
         ENDIF
         B_PREV = B_F
         B_I = B_F
         PERIOD_I = PERIOD_F
         PDOT_I = PDOT_F
   50 CONTINUE

C     --- Stable RLOF accretion: no negative mass/period/field, no NaN
      PSR_PMODE = 1
      PSR_BMODE = 1
      PSR_SPINA = 0.1D0
      PSR_BMAGA = 12.0D0
      CALL PULSARINIT(PULSARINITB, PULSARINITP, PULSARINITPD,
     &     PSRMASS, .FALSE.)

      B_I = PULSARINITB
      PERIOD_I = PULSARINITP
      PDOT_I = PULSARINITPD
      MASS = PSRMASS
      MASS_COMP = 1.0D0
      RAD_COMP = 1.0D0
      PSRMDOT = 1.0D-9
      NSINCL = 1.5708D0
      CE = .FALSE.
      DTM_STEP = 1.0D4 * 3.1536D13

      DO 60 ISTEP = 1, 20
         CALL PULSAREVO(MASS, MASS_COMP, RAD_COMP, B_I, PERIOD_I,
     &        PDOT_I, PSRMDOT, DTM_STEP, NSINCL, CE, 2,
     &        B_F, PERIOD_F, PDOT_F)
         CALL CHKPOS(10, ISTEP, 0, B_F, NFAIL)
         CALL CHKPOS(11, ISTEP, 0, PERIOD_F, NFAIL)
         CALL CHKPOS(12, ISTEP, 0, PDOT_F, NFAIL)
         B_I = B_F
         PERIOD_I = PERIOD_F
         PDOT_I = PDOT_F
   60 CONTINUE

C     --- PSR_ACC_CE=0 means "CE behaves exactly like isolated
C     evolution" (plan section 6.3): with the same starting state,
C     ce=.TRUE.+psrmdot>0 must give bit-identical output to
C     ce=.FALSE.+psrmdot=0, since both should land in the same first
C     branch of PULSAREVO.
      PSR_PMODE = 1
      PSR_BMODE = 1
      PSR_SPINA = 0.1D0
      PSR_BMAGA = 13.0D0
      CALL PULSARINIT(PULSARINITB, PULSARINITP, PULSARINITPD,
     &     PSRMASS, .FALSE.)

      B_I = PULSARINITB
      PERIOD_I = PULSARINITP
      PDOT_I = PULSARINITPD
      MASS = PSRMASS
      MASS_COMP = 0.0D0
      RAD_COMP = 0.0D0
      NSINCL = 1.5708D0
      DTM_STEP = 1.0D6 * 3.1536D13

      PSRMDOT = 0.0D0
      CE = .FALSE.
      CALL PULSAREVO(MASS, MASS_COMP, RAD_COMP, B_I, PERIOD_I,
     &     PDOT_I, PSRMDOT, DTM_STEP, NSINCL, CE, 3,
     &     B_F, PERIOD_F, PDOT_F)
      B_ISO = B_F
      P_ISO = PERIOD_F
      PD_ISO = PDOT_F

C     PULSAREVO mutates its B_i/period_i/pdot_i dummies in place while
C     splitting the isolated branch into substeps, so the actuals
C     (B_I/PERIOD_I/PDOT_I here) must be reset to the true starting
C     state before reusing them for the second call.
      B_I = PULSARINITB
      PERIOD_I = PULSARINITP
      PDOT_I = PULSARINITPD
      MASS = PSRMASS
      PSRMDOT = 1.0D-9
      CE = .TRUE.
      CALL PULSAREVO(MASS, MASS_COMP, RAD_COMP, B_I, PERIOD_I,
     &     PDOT_I, PSRMDOT, DTM_STEP, NSINCL, CE, 4,
     &     B_F, PERIOD_F, PDOT_F)

      IF (B_F.NE.B_ISO .OR. PERIOD_F.NE.P_ISO
     &     .OR. PDOT_F.NE.PD_ISO) THEN
         WRITE(6,*) 'FAIL: PSR_ACC_CE=0 CE result differs from',
     &        ' isolated evolution:', B_ISO, B_F, P_ISO, PERIOD_F,
     &        PD_ISO, PDOT_F
         NFAIL = NFAIL + 1
      ENDIF

      IF (NFAIL.EQ.0) THEN
         WRITE(6,*) 'All pulsar physics invariant checks passed.'
         CALL EXIT(0)
      ELSE
         WRITE(6,*) NFAIL, 'invariant check(s) FAILED.'
         CALL EXIT(1)
      ENDIF

      END

      SUBROUTINE CHKPOS(TAG, A, B, X, NFAIL)
C     Fails if X is not finite and strictly positive. TAG identifies
C     which check this is; A/B are context (mode numbers or step).
      IMPLICIT NONE
      INTEGER TAG, A, B, NFAIL
      REAL*8 X
      LOGICAL BAD

      BAD = .FALSE.
      IF (X.NE.X) BAD = .TRUE.
      IF (X.LE.0.0D0) BAD = .TRUE.
      IF (ABS(X).GT.HUGE(1.0D0)/1.0D3) BAD = .TRUE.
      IF (BAD) THEN
         WRITE(6,*) 'FAIL: check', TAG, 'A=', A, 'B=', B,
     &        'not finite/positive:', X
         NFAIL = NFAIL + 1
      ENDIF
      RETURN
      END
