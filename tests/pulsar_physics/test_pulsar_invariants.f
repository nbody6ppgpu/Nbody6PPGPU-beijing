C     Standalone PULSARINIT/PULSAREVO invariant and registration checks,
C     linked directly against the built pulsar.o/ran2.o. Exercises:
C       - non-NS -> NS registration, AIC classification, and dedup
C       - KZ(29)=0/1/2 transition-registration behavior
C       - finalized merger registration, rebinding, and slot retirement
C       - the plan's section 10.4 physics invariants:
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
      REAL*8 B_PREV, DTM_STEP, BMIN_TESLA, B_I_PRECALL
      REAL*8 B_ISO, P_ISO, PD_ISO
      REAL*8 EVENTTIME1, EVENTTIME2
      REAL*8 EVENTXM1, EVENTBM1, EVENTPER1, EVENTPD1
      REAL*8 EVENTAX1, EVENTA01, EVENTBM01, EVENTXM01, EVENTTAU1
      REAL*8 EVENTXM2, EVENTBM2, EVENTPER2, EVENTPD2
      REAL*8 EVENTAX2, EVENTA02, EVENTBM02, EVENTXM02, EVENTTAU2
      INTEGER ISTEP, EVENTCODE1, EVENTCODE2, EVENTNAME1, EVENTNAME2
      INTEGER EVENTIOS, PTEST
      LOGICAL CE, REGISTD
      CHARACTER*256 EVENT1, EVENT2

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

C     --- RLOF sign/units bracket for psrmdot, at a realistic step ---
C     The loop above uses DTM_STEP=1e4 Myr as a numerical stress test
C     (deliberately large enough to collapse B to PSR_BMIN in one call
C     and confirm that stays finite/positive, not NaN) -- too extreme
C     for a burial-direction/magnitude check, since it saturates B to
C     BMIN regardless of psrmdot's sign or units. Use a step closer to
C     roche.f's actual per-call DTM (a stellar-evolution substep, order
C     Myr or less) instead.
C     roche.f's two live call sites both compute
C     PSRMDOT = DM22_or_DM2/(DTM*myr_in_sec); a wrong sign or a missed
C     *TSTAR/myr_in_sec unit factor there would violate this bracket
C     even though the loop above's positivity/no-NaN checks would not
C     catch it. NSACCRETE's B update is
C     B_f = BMIN + EXP(-psrmdot*dtm/KAPPA)*(B_i-BMIN) -- accretion
C     buries the field toward BMIN, it does not grow it. With
C     psrmdot, dtm > 0 this must hold:
C       B_f < B_i          (burial, not growth)
C       B_f > BMIN_TESLA   (does not overshoot the floor in one step --
C                            empirically confirmed to catch psrmdot too
C                            large by ~myr_in_sec)
C       B_f not within round-off of B_i (a weak backstop against a
C                            near-zero/near-cancelled psrmdot; NOT
C                            calibrated to catch a moderate, e.g.
C                            10x-1e4x, magnitude error -- see below)
C     NOTE on what this bracket does NOT catch: pulsar.F has its own
C     branch-selection threshold `acc=1e-20` (Msun/s) gating which of
C     ISOSPINDOWN/NSACCRETE runs. A sign-flipped psrmdot, or one
C     undershot by more than ~1e17 (crossing below acc), reroutes
C     silently to isolated evolution -- which also decreases B via a
C     different mechanism (ISOSPINDOWN), so it can still pass the
C     "B_f < B_i" check above for the wrong reason. Confirmed this
C     empirically: neither a sign-flipped nor a ~3e13x-undershot
C     psrmdot triggers any of the three checks below. Only a magnitude
C     error that stays within the accretion branch (i.e. psrmdot stays
C     positive and above 1e-20) is actually exercised here; this is
C     the realistic error shape for roche.f's call sites anyway (both
C     compute PSRMDOT = DM22_or_DM2/(DTM*myr_in_sec) from
C     already-positive BSE quantities, so a missing *TSTAR or
C     *myr_in_sec factor makes PSRMDOT too LARGE by that factor, not
C     negative or near-zero) -- see source-map.md for this finding.
C     PULSAREVO mutates its B_i/period_i/pdot_i dummies in place during
C     its internal sub-stepping (source-map.md item 7's aliasing
C     footgun -- confirmed here to also hit the accretion branch, not
C     just isolated/CE), so B_I itself is not reliable as "the value
C     before this call" once the call returns -- snapshot it first.
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
C     PULSAR_KAPPA=PSR_KAPPA*MSOL_TO_KG=1e-3*1.989e30 kg, dtm=1 Myr;
C     PSRMDOT=1.0D-9 (used in the stress-test loop above, chosen to
C     combine with its 1e4 Myr step) makes psrmdot*dtm/KAPPA ~3e16 --
C     saturates EXP() to 0 even at dtm=1 Myr. Need a realistic
C     accretion rate here instead: 1.5D-17 Msun/s (~5e-10 Msun/yr,
C     sub-Eddington for a NS) gives psrmdot*dtm/KAPPA ~0.5, a partial
C     (non-saturated, non-negligible) decay -- the actual regime this
C     bracket needs to discriminate in.
      PSRMDOT = 1.5D-17
      NSINCL = 1.5708D0
      CE = .FALSE.
      DTM_STEP = 1.0D0 * 3.1536D13
      B_I_PRECALL = B_I
      CALL PULSAREVO(MASS, MASS_COMP, RAD_COMP, B_I, PERIOD_I,
     &     PDOT_I, PSRMDOT, DTM_STEP, NSINCL, CE, 16,
     &     B_F, PERIOD_F, PDOT_F)
      IF (B_F.GE.B_I_PRECALL) THEN
         WRITE(6,*) 'FAIL: RLOF bracket: B did not decrease',
     &        ' (burial) B_i,B_f=', B_I_PRECALL, B_F
         NFAIL = NFAIL + 1
      ENDIF
      IF (B_F.LE.BMIN_TESLA) THEN
         WRITE(6,*) 'FAIL: RLOF bracket: B collapsed to/below BMIN',
     &        ' in one step, B_f,BMIN=', B_F, BMIN_TESLA
         NFAIL = NFAIL + 1
      ENDIF
      IF (ABS(B_F-B_I_PRECALL).LT.1.0D-12*ABS(B_I_PRECALL)) THEN
         WRITE(6,*) 'FAIL: RLOF bracket: B did not move (psrmdot',
     &        ' units too small?) B_i,B_f=', B_I_PRECALL, B_F
         NFAIL = NFAIL + 1
      ENDIF

C     --- Non-positive dtm must not silently produce Inf/NaN (it did,
C     before the dtm.LE.0 guard was added to PULSAREVO: EXP(-step/TAU)
C     with step<0 exceeds 1, B grows, and Pfsquared can go negative,
C     so SQRT gives NaN). Lifecycle hooks compute dtm as a TPHYS
C     difference, which forward-time integration should keep
C     non-negative, but this is exactly the kind of invariant that
C     should not be allowed to fail silently if it ever is violated.
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
      DTM_STEP = -1.0D6 * 3.1536D13
      CALL PULSAREVO(MASS, MASS_COMP, RAD_COMP, B_I, PERIOD_I,
     &     PDOT_I, PSRMDOT, DTM_STEP, NSINCL, CE, 5,
     &     B_F, PERIOD_F, PDOT_F)
      CALL CHKPOS(13, 0, 0, B_F, NFAIL)
      CALL CHKPOS(14, 0, 0, PERIOD_F, NFAIL)
      CALL CHKPOS(15, 0, 0, PDOT_F, NFAIL)
      IF (B_F.NE.B_I .OR. PERIOD_F.NE.PERIOD_I
     &     .OR. PDOT_F.NE.PDOT_I) THEN
         WRITE(6,*) 'FAIL: negative dtm did not return inputs',
     &        ' unchanged:', B_I, B_F, PERIOD_I, PERIOD_F,
     &        PDOT_I, PDOT_F
         NFAIL = NFAIL + 1
      ENDIF

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

C     --- NS transition registration and AIC distribution selection ---
C     Use distinct fixed periods so the selected distribution is
C     observable without depending on random draws.
      PSR_PMODE = 1
      PSR_BMODE = 1
      PSR_SPINA = 0.5D0
      PSR_BMAGA = 12.0D0
      PSRM_PMODE = 1
      PSRM_BMODE = 1
      PSRM_SPINA = 0.01D0
      PSRM_BMAGA = 9.0D0
      NSCOUNT = 0
      KZ(50) = 0
      rank = 0

C     Disabled pulsar evolution must not register anything.
      KZ(29) = 0
      CALL PSRREG_TRANSITION(101,1.30D0,12,13,10.0D0,REGISTD)
      IF (REGISTD.OR.NSCOUNT.NE.0) THEN
         WRITE(6,*) 'FAIL: KZ(29)=0 registered an NS transition'
         NFAIL = NFAIL + 1
      ENDIF

C     A type that has not become an NS must not be registered.
      KZ(29) = 2
      CALL PSRREG_TRANSITION(101,1.30D0,12,12,10.0D0,REGISTD)
      IF (REGISTD.OR.NSCOUNT.NE.0) THEN
         WRITE(6,*) 'FAIL: non-NS transition registered a pulsar'
         NFAIL = NFAIL + 1
      ENDIF

C     A non-AIC NS transition uses the ordinary PSR_* distribution.
      CALL PSRREG_TRANSITION(102,1.40D0,8,13,11.0D0,REGISTD)
      IF (.NOT.REGISTD.OR.NSCOUNT.NE.1) THEN
         WRITE(6,*) 'FAIL: ordinary NS transition not registered'
         NFAIL = NFAIL + 1
      ELSEIF (ABS(PERIODNS(1)-PSR_SPINA).GT.1.0D-12) THEN
         WRITE(6,*) 'FAIL: ordinary NS used wrong birth distribution'
         NFAIL = NFAIL + 1
      ENDIF

C     With KZ(29)=2, ONeWD -> NS uses the PSRM_* distribution.
      CALL PSRREG_TRANSITION(103,1.26D0,12,13,12.0D0,REGISTD)
      IF (.NOT.REGISTD.OR.NSCOUNT.NE.2) THEN
         WRITE(6,*) 'FAIL: AIC NS transition not registered'
         NFAIL = NFAIL + 1
      ELSEIF (ABS(PERIODNS(2)-PSRM_SPINA).GT.1.0D-12) THEN
         WRITE(6,*) 'FAIL: AIC NS used wrong birth distribution'
         NFAIL = NFAIL + 1
      ENDIF

C     Repeated lifecycle hooks must not allocate a second pulsar slot.
      CALL PSRREG_TRANSITION(103,1.26D0,12,13,12.0D0,REGISTD)
      IF (REGISTD.OR.NSCOUNT.NE.2) THEN
         WRITE(6,*) 'FAIL: duplicate NS transition registered twice'
         NFAIL = NFAIL + 1
      ENDIF

C     With KZ(29)=1, even ONeWD -> NS uses the regular distribution.
      KZ(29) = 1
      CALL PSRREG_TRANSITION(104,1.26D0,12,13,13.0D0,REGISTD)
      IF (.NOT.REGISTD.OR.NSCOUNT.NE.3) THEN
         WRITE(6,*) 'FAIL: KZ(29)=1 AIC transition not registered'
         NFAIL = NFAIL + 1
      ELSEIF (ABS(PERIODNS(3)-PSR_SPINA).GT.1.0D-12) THEN
         WRITE(6,*) 'FAIL: KZ(29)=1 used AIC birth distribution'
         NFAIL = NFAIL + 1
      ENDIF

C     A merger-produced NS without an ONeWD uses the regular distribution.
      KZ(29) = 2
      CALL PSRREG_MERGER(201,202,202,1.35D0,11,1,13,14.0D0,
     &     REGISTD)
      IF (.NOT.REGISTD.OR.NSCOUNT.NE.4) THEN
         WRITE(6,*) 'FAIL: ordinary merger NS not registered'
         NFAIL = NFAIL + 1
      ELSEIF (ABS(PERIODNS(4)-PSR_SPINA).GT.1.0D-12) THEN
         WRITE(6,*) 'FAIL: ordinary merger used wrong distribution'
         NFAIL = NFAIL + 1
      ENDIF

C     An ONeWD progenitor selects the AIC distribution after a merger.
      CALL PSRREG_MERGER(203,204,204,1.26D0,12,11,13,15.0D0,
     &     REGISTD)
      IF (.NOT.REGISTD.OR.NSCOUNT.NE.5) THEN
         WRITE(6,*) 'FAIL: AIC merger NS not registered'
         NFAIL = NFAIL + 1
      ELSEIF (ABS(PERIODNS(5)-PSRM_SPINA).GT.1.0D-12) THEN
         WRITE(6,*) 'FAIL: AIC merger used wrong distribution'
         NFAIL = NFAIL + 1
      ENDIF

C     A merger whose final type is not an NS must not register a pulsar.
      CALL PSRREG_MERGER(205,206,205,2.0D0,12,11,14,16.0D0,
     &     REGISTD)
      IF (REGISTD.OR.NSCOUNT.NE.5) THEN
         WRITE(6,*) 'FAIL: non-NS merger registered a pulsar'
         NFAIL = NFAIL + 1
      ENDIF

C     If survivor-name selection discards an existing NS name, carry its
C     state into the final NAME instead of initializing a second pulsar.
      CALL PSRREG_TRANSITION(210,1.40D0,8,13,16.0D0,REGISTD)
      PERIODNS(6) = 0.123D0
      BMAGNS(6) = 4.56D8
      PDOTNS(6) = 7.89D-18
      CALL PSRREG_MERGER(211,210,211,1.55D0,12,13,13,17.0D0,
     &     REGISTD)
      IF (REGISTD.OR.NSCOUNT.NE.6.OR.NAMENS(6).NE.211) THEN
         WRITE(6,*) 'FAIL: merger did not rebind existing pulsar'
         NFAIL = NFAIL + 1
      ELSEIF (ABS(PERIODNS(6)-0.123D0).GT.1.0D-12.OR.
     &        ABS(BMAGNS(6)-4.56D8).GT.1.0D-6.OR.
     &        ABS(PDOTNS(6)-7.89D-18).GT.1.0D-28.OR.
     &        ABS(AGENSX(6)-16.0D0).GT.1.0D-12) THEN
         WRITE(6,*) 'FAIL: merger reset existing pulsar state'
         NFAIL = NFAIL + 1
      ELSEIF (ABS(XMNS(6)-1.55D0).GT.1.0D-12.OR.
     &        NSTYPE(6).NE.13) THEN
         WRITE(6,*) 'FAIL: merger did not update pulsar remnant'
         NFAIL = NFAIL + 1
      ENDIF

C     If an existing NS is destroyed in a non-NS merger, its registry
C     slot must be removed rather than left under a vanished NAME.
      CALL PSRREG_MERGER(211,212,211,2.0D0,13,1,14,18.0D0,
     &     REGISTD)
      IF (REGISTD.OR.NSCOUNT.NE.5) THEN
         WRITE(6,*) 'FAIL: non-NS remnant did not retire NS slot'
         NFAIL = NFAIL + 1
      ENDIF
      DO ISTEP = 1, NSCOUNT
         IF (NAMENS(ISTEP).EQ.211) THEN
            WRITE(6,*) 'FAIL: destroyed NS NAME remains registered'
            NFAIL = NFAIL + 1
         ENDIF
      END DO

C     A double-NS merger that remains an NS must retain one progenitor
C     state, update the remnant, and compact away the other slot.
      CALL PSRREG_TRANSITION(220,1.30D0,8,13,19.0D0,REGISTD)
      CALL PSRREG_TRANSITION(221,1.40D0,8,13,20.0D0,REGISTD)
      PERIODNS(6) = 0.220D0
      PERIODNS(7) = 0.221D0
      BMAGNS(7) = 2.21D8
      PDOTNS(7) = 2.21D-18
      CALL PSRREG_MERGER(220,221,221,2.60D0,13,13,13,21.0D0,
     &     REGISTD)
      IF (REGISTD.OR.NSCOUNT.NE.6.OR.NAMENS(6).NE.221) THEN
         WRITE(6,*) 'FAIL: double-NS merger left wrong registry'
         NFAIL = NFAIL + 1
      ELSEIF (ABS(PERIODNS(6)-0.221D0).GT.1.0D-12.OR.
     &        ABS(BMAGNS(6)-2.21D8).GT.1.0D-6.OR.
     &        ABS(PDOTNS(6)-2.21D-18).GT.1.0D-28.OR.
     &        ABS(XMNS(6)-2.60D0).GT.1.0D-12) THEN
         WRITE(6,*) 'FAIL: double-NS survivor state was not preserved'
         NFAIL = NFAIL + 1
      ENDIF
      DO ISTEP = 1, NSCOUNT
         IF (NAMENS(ISTEP).EQ.220) THEN
            WRITE(6,*) 'FAIL: double-NS merger left orphan slot'
            NFAIL = NFAIL + 1
         ENDIF
      END DO

C     A double-NS merger with a non-NS remnant retires both slots.
      CALL PSRREG_TRANSITION(230,1.30D0,8,13,22.0D0,REGISTD)
      CALL PSRREG_TRANSITION(231,1.40D0,8,13,23.0D0,REGISTD)
      CALL PSRREG_MERGER(230,231,231,2.70D0,13,13,14,24.0D0,
     &     REGISTD)
      IF (REGISTD.OR.NSCOUNT.NE.6) THEN
         WRITE(6,*) 'FAIL: non-NS double merger kept pulsar slots'
         NFAIL = NFAIL + 1
      ENDIF
      DO ISTEP = 1, NSCOUNT
         IF (NAMENS(ISTEP).EQ.230.OR.NAMENS(ISTEP).EQ.231) THEN
            WRITE(6,*) 'FAIL: non-NS merger left orphan NS NAME'
            NFAIL = NFAIL + 1
         ENDIF
      END DO

C     A NAME migration emits adjacent code-5 snapshots with old and new
C     identities.  The old record precedes the mass update; both records
C     retain all other state and use the event-bank physical time.
      CALL PSRREG_TRANSITION(240,1.40D0,8,13,25.0D0,REGISTD)
      PTEST = NSCOUNT
      BMAGNS(PTEST) = 2.40D8
      PERIODNS(PTEST) = 0.240D0
      PDOTNS(PTEST) = 2.40D-18
      AGENSX(PTEST) = 25.0D0
      AGENS0(PTEST) = 2.50D0
      BMAGNS0(PTEST) = 2.40D7
      XMNS0(PTEST) = 1.24D0
      PSR_TAU = 3.50D0
      OPEN(UNIT=204,STATUS='SCRATCH',FORM='FORMATTED')
      KZ(50) = 1
      TTOT = 2.0D0
      TSTAR = 3.0D0
      CALL PSRREG_MERGER(241,240,241,1.60D0,1,13,13,26.0D0,
     &     REGISTD)
      REWIND 204
      READ(204,'(A)',IOSTAT=EVENTIOS) EVENT1
      IF (EVENTIOS.NE.0) THEN
         WRITE(6,*) 'FAIL: missing old-NAME merger event'
         NFAIL = NFAIL + 1
      ELSE
         READ(204,'(A)',IOSTAT=EVENTIOS) EVENT2
         IF (EVENTIOS.NE.0) THEN
            WRITE(6,*) 'FAIL: missing new-NAME merger event'
            NFAIL = NFAIL + 1
         ELSE
            READ(EVENT1(8:),*) EVENTCODE1, EVENTTIME1,
     &           EVENTNAME1, EVENTXM1, EVENTBM1, EVENTPER1,
     &           EVENTPD1, EVENTAX1, EVENTA01, EVENTBM01,
     &           EVENTXM01, EVENTTAU1
            READ(EVENT2(8:),*) EVENTCODE2, EVENTTIME2,
     &           EVENTNAME2, EVENTXM2, EVENTBM2, EVENTPER2,
     &           EVENTPD2, EVENTAX2, EVENTA02, EVENTBM02,
     &           EVENTXM02, EVENTTAU2
            IF (EVENTCODE1.NE.5.OR.EVENTCODE2.NE.5.OR.
     &          EVENTNAME1.NE.240.OR.EVENTNAME2.NE.241) THEN
               WRITE(6,*) 'FAIL: merger NAME event pair is wrong'
               NFAIL = NFAIL + 1
            ENDIF
            IF (ABS(EVENTTIME1-6.0D0).GT.1.0D-10.OR.
     &          ABS(EVENTTIME2-6.0D0).GT.1.0D-10) THEN
               WRITE(6,*) 'FAIL: merger NAME event time is wrong'
               NFAIL = NFAIL + 1
            ENDIF
            IF (ABS(EVENTXM1-1.40D0).GT.1.0D-10.OR.
     &          ABS(EVENTXM2-1.60D0).GT.1.0D-10) THEN
               WRITE(6,*) 'FAIL: merger NAME event mass is wrong'
               NFAIL = NFAIL + 1
            ENDIF
            IF (ABS(EVENTBM1-2.40D8).GT.1.0D0.OR.
     &          ABS(EVENTPER1-0.240D0).GT.1.0D-10.OR.
     &          ABS(EVENTPD1-2.40D-18).GT.1.0D-26.OR.
     &          ABS(EVENTAX1-25.0D0).GT.1.0D-10.OR.
     &          ABS(EVENTA01-2.50D0).GT.1.0D-10.OR.
     &          ABS(EVENTBM01-2.40D7).GT.1.0D0.OR.
     &          ABS(EVENTXM01-1.24D0).GT.1.0D-10.OR.
     &          ABS(EVENTTAU1-3.50D0).GT.1.0D-10.OR.
     &          ABS(EVENTBM2-2.40D8).GT.1.0D0.OR.
     &          ABS(EVENTPER2-0.240D0).GT.1.0D-10.OR.
     &          ABS(EVENTPD2-2.40D-18).GT.1.0D-26.OR.
     &          ABS(EVENTAX2-25.0D0).GT.1.0D-10.OR.
     &          ABS(EVENTA02-2.50D0).GT.1.0D-10.OR.
     &          ABS(EVENTBM02-2.40D7).GT.1.0D0.OR.
     &          ABS(EVENTXM02-1.24D0).GT.1.0D-10.OR.
     &          ABS(EVENTTAU2-3.50D0).GT.1.0D-10) THEN
               WRITE(6,*) 'FAIL: merger NAME event state is wrong'
               NFAIL = NFAIL + 1
            ENDIF
         ENDIF
      ENDIF
      CLOSE(204)
      KZ(50) = 0

C     With no NAME change, reconciliation writes one code-5 event.
      CALL PSRREG_TRANSITION(250,1.50D0,8,13,27.0D0,REGISTD)
      PTEST = NSCOUNT
      BMAGNS(PTEST) = 2.50D8
      PERIODNS(PTEST) = 0.250D0
      PDOTNS(PTEST) = 2.50D-18
      AGENSX(PTEST) = 27.0D0
      AGENS0(PTEST) = 2.70D0
      BMAGNS0(PTEST) = 2.50D7
      XMNS0(PTEST) = 1.25D0
      PSR_TAU = 4.50D0
      OPEN(UNIT=204,STATUS='SCRATCH',FORM='FORMATTED')
      KZ(50) = 1
      TTOT = 4.0D0
      TSTAR = 2.0D0
      CALL PSRREG_MERGER(250,251,250,1.75D0,13,1,13,28.0D0,
     &     REGISTD)
      REWIND 204
      READ(204,'(A)',IOSTAT=EVENTIOS) EVENT1
      IF (EVENTIOS.NE.0) THEN
         WRITE(6,*) 'FAIL: missing unchanged-NAME merger event'
         NFAIL = NFAIL + 1
      ELSE
         READ(EVENT1(8:),*) EVENTCODE1, EVENTTIME1,
     &        EVENTNAME1, EVENTXM1, EVENTBM1, EVENTPER1,
     &        EVENTPD1, EVENTAX1, EVENTA01, EVENTBM01,
     &        EVENTXM01, EVENTTAU1
         IF (EVENTCODE1.NE.5.OR.EVENTNAME1.NE.250.OR.
     &       ABS(EVENTTIME1-8.0D0).GT.1.0D-10.OR.
     &       ABS(EVENTXM1-1.75D0).GT.1.0D-10.OR.
     &       ABS(EVENTBM1-2.50D8).GT.1.0D0.OR.
     &       ABS(EVENTPER1-0.250D0).GT.1.0D-10.OR.
     &       ABS(EVENTPD1-2.50D-18).GT.1.0D-26.OR.
     &       ABS(EVENTAX1-27.0D0).GT.1.0D-10.OR.
     &       ABS(EVENTA01-2.70D0).GT.1.0D-10.OR.
     &       ABS(EVENTBM01-2.50D7).GT.1.0D0.OR.
     &       ABS(EVENTXM01-1.25D0).GT.1.0D-10.OR.
     &       ABS(EVENTTAU1-4.50D0).GT.1.0D-10) THEN
            WRITE(6,*) 'FAIL: unchanged-NAME merger event is wrong'
            NFAIL = NFAIL + 1
         ENDIF
         READ(204,'(A)',IOSTAT=EVENTIOS) EVENT2
         IF (EVENTIOS.GE.0) THEN
            WRITE(6,*) 'FAIL: unchanged NAME wrote extra code-5 event'
            NFAIL = NFAIL + 1
         ENDIF
      ENDIF
      CLOSE(204)
      KZ(50) = 0

C     Removing an interior slot writes its old state before compaction,
C     moves every parallel array, and clears every old tail-slot field.
      NSCOUNT = 3
      NAMENS(2) = 602
      XMNS(2) = 1.62D0
      NSTYPE(2) = 62
      NSSTAT(2) = 6020
      BMAGNS(2) = 6.02D8
      PERIODNS(2) = 0.602D0
      PDOTNS(2) = 6.02D-18
      AGENSX(2) = 62.0D0
      AGENS0(2) = 6.20D0
      BMAGNS0(2) = 6.02D7
      XMNS0(2) = 1.02D0
      NAMENS(3) = 703
      XMNS(3) = 1.73D0
      NSTYPE(3) = 73
      NSSTAT(3) = 7030
      BMAGNS(3) = 7.03D8
      PERIODNS(3) = 0.703D0
      PDOTNS(3) = 7.03D-18
      AGENSX(3) = 73.0D0
      AGENS0(3) = 7.30D0
      BMAGNS0(3) = 7.03D7
      XMNS0(3) = 1.03D0
      PSR_TAU = 5.50D0
      OPEN(UNIT=204,STATUS='SCRATCH',FORM='FORMATTED')
      KZ(50) = 1
      TTOT = 5.0D0
      TSTAR = 2.0D0
      CALL PSRREG_REMOVE(2)
      REWIND 204
      READ(204,'(A)',IOSTAT=EVENTIOS) EVENT1
      IF (EVENTIOS.NE.0) THEN
         WRITE(6,*) 'FAIL: slot retirement wrote no code-6 event'
         NFAIL = NFAIL + 1
      ELSE
         READ(EVENT1(8:),*) EVENTCODE1, EVENTTIME1,
     &        EVENTNAME1, EVENTXM1, EVENTBM1, EVENTPER1,
     &        EVENTPD1, EVENTAX1, EVENTA01, EVENTBM01,
     &        EVENTXM01, EVENTTAU1
         IF (EVENTCODE1.NE.6.OR.EVENTNAME1.NE.602.OR.
     &       ABS(EVENTTIME1-10.0D0).GT.1.0D-10.OR.
     &       ABS(EVENTXM1-1.62D0).GT.1.0D-10.OR.
     &       ABS(EVENTBM1-6.02D8).GT.1.0D0.OR.
     &       ABS(EVENTPER1-0.602D0).GT.1.0D-10.OR.
     &       ABS(EVENTPD1-6.02D-18).GT.1.0D-26.OR.
     &       ABS(EVENTAX1-62.0D0).GT.1.0D-10.OR.
     &       ABS(EVENTA01-6.20D0).GT.1.0D-10.OR.
     &       ABS(EVENTBM01-6.02D7).GT.1.0D0.OR.
     &       ABS(EVENTXM01-1.02D0).GT.1.0D-10.OR.
     &       ABS(EVENTTAU1-5.50D0).GT.1.0D-10) THEN
            WRITE(6,*) 'FAIL: code-6 event used post-move state'
            NFAIL = NFAIL + 1
         ENDIF
      ENDIF
      CLOSE(204)
      IF (NSCOUNT.NE.2.OR.NAMENS(2).NE.703.OR.
     &    ABS(XMNS(2)-1.73D0).GT.1.0D-12.OR.
     &    NSTYPE(2).NE.73.OR.NSSTAT(2).NE.7030.OR.
     &    ABS(BMAGNS(2)-7.03D8).GT.1.0D-6.OR.
     &    ABS(PERIODNS(2)-0.703D0).GT.1.0D-12.OR.
     &    ABS(PDOTNS(2)-7.03D-18).GT.1.0D-28.OR.
     &    ABS(AGENSX(2)-73.0D0).GT.1.0D-12.OR.
     &    ABS(AGENS0(2)-7.30D0).GT.1.0D-12.OR.
     &    ABS(BMAGNS0(2)-7.03D7).GT.1.0D-6.OR.
     &    ABS(XMNS0(2)-1.03D0).GT.1.0D-12) THEN
         WRITE(6,*) 'FAIL: slot compaction missed a registry array'
         NFAIL = NFAIL + 1
      ENDIF
      IF (NAMENS(3).NE.0.OR.XMNS(3).NE.0.0D0.OR.
     &    NSTYPE(3).NE.0.OR.NSSTAT(3).NE.0.OR.
     &    BMAGNS(3).NE.0.0D0.OR.PERIODNS(3).NE.0.0D0.OR.
     &    PDOTNS(3).NE.0.0D0.OR.AGENSX(3).NE.0.0D0.OR.
     &    AGENS0(3).NE.0.0D0.OR.BMAGNS0(3).NE.0.0D0.OR.
     &    XMNS0(3).NE.0.0D0) THEN
         WRITE(6,*) 'FAIL: retired registry tail slot was not cleared'
         NFAIL = NFAIL + 1
      ENDIF
      KZ(50) = 0

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
