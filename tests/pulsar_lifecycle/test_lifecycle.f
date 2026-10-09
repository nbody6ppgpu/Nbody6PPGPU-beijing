      PROGRAM TEST_PSR_LIFECYCLE
      INCLUDE 'common6.h'
      INCLUDE 'pulsar_ctrl.h'
      INTEGER RNGWORDS(99), SAVEDRNG(99), SAVEDSEED
      COMMON/RAND2/ RNGWORDS
      INTEGER NBAD, NUNRESOLVED, CODE, EVENTNAME, IOS, J
      REAL*8 EVENTTIME
      LOGICAL REGISTD
      CHARACTER*256 RECORD
      rank = 0
      KZ(29) = 2
      KZ(50) = 0
      IDUM1 = -12345
      PSR_PMODE = 3
      PSR_BMODE = 5
      PSR_SPINA = 0.075D0
      PSR_SPINB = 0.025D0
      PSR_BMAGA = 11.5D0
      PSR_BMAGB = 0.5D0
      PSR_BMIN = 7.D0
      PSR_RAD = 10.D0
      PSR_TAU = 1000.D0
      NSCOUNT = 0
      CALL PSRREG_TRANSITION(10,1.4D0,8,13,0.D0,REGISTD)
      CALL PSRREG_TRANSITION(20,2.5D0,8,13,0.D0,REGISTD)
      CALL PSRREG_TRANSITION(30,1.3D0,8,13,0.D0,REGISTD)
      SAVEDSEED = IDUM1
      SAVEDRNG = RNGWORDS
      OPEN(204,STATUS='SCRATCH',FORM='FORMATTED')
      KZ(50) = 1
      TTOT = 2.D0
      TSTAR = 3.D0
      CALL PSRREG_TRANSITION(20,2.51D0,13,14,6.D0,REGISTD)
      IF (REGISTD.OR.NSCOUNT.NE.2) STOP 1
      IF (NAMENS(1).NE.10.OR.NAMENS(2).NE.30) STOP 2
      IF (XMNS(2).NE.1.3D0.OR.NAMENS(3).NE.0) STOP 3
C     Repetition, missing NAME and unchanged NS must not retire another.
      CALL PSRREG_TRANSITION(20,2.51D0,13,14,6.D0,REGISTD)
      CALL PSRREG_TRANSITION(99,2.51D0,13,14,6.D0,REGISTD)
      CALL PSRREG_TRANSITION(10,1.4D0,13,13,6.D0,REGISTD)
      IF (REGISTD.OR.NSCOUNT.NE.2) STOP 4
      IF (SAVEDSEED.NE.IDUM1.OR.ANY(SAVEDRNG.NE.RNGWORDS)) STOP 5
      REWIND 204
      READ(204,'(A)',IOSTAT=IOS) RECORD
      IF (IOS.NE.0) STOP 6
      READ(RECORD(8:),*) CODE,EVENTTIME,EVENTNAME
      IF (CODE.NE.6.OR.EVENTNAME.NE.20.OR.EVENTTIME.NE.6.D0)
     &   STOP 7
      READ(204,'(A)',IOSTAT=IOS) RECORD
      IF (IOS.GE.0) STOP 8
      CLOSE(204)
      KZ(50) = 0
C     Disabled mode has no registry or RNG effects.
      KZ(29) = 0
      CALL PSRREG_TRANSITION(10,2.51D0,13,14,6.D0,REGISTD)
      IF (REGISTD.OR.NSCOUNT.NE.2) STOP 9
      KZ(29) = 2
C     Audit only resolves physical positive-mass entries in 1:N.
      N = 4
      NZERO = 50
      DO J = 1, 4
         BODY(J) = 1.D0
         NAME(J) = 10*J
         KSTAR(J) = 13
      END DO
      KSTAR(2) = 14
      NAME(4) = 0
      CALL PSRREG_AUDIT(NBAD,NUNRESOLVED)
      IF (NBAD.NE.0.OR.NUNRESOLVED.NE.0) STOP 10
C     A stale BH record is a definite inconsistency.
      KSTAR(3) = 14
      CALL PSRREG_AUDIT(NBAD,NUNRESOLVED)
      IF (NBAD.NE.1) STOP 11
      KSTAR(3) = 13
C     Missing pulsar history, including a resolved KS component.
      NSCOUNT = 1
      CALL PSRREG_AUDIT(NBAD,NUNRESOLVED)
      IF (NBAD.NE.1) STOP 12
      NSCOUNT = 2
C     Absent escaper and zero-mass ghost records remain unresolved.
      NAME(3) = 31
      KSTAR(3) = 1
      CALL PSRREG_AUDIT(NBAD,NUNRESOLVED)
      IF (NBAD.NE.0.OR.NUNRESOLVED.NE.1) STOP 13
      NAME(3) = 30
      BODY(3) = 0.D0
      CALL PSRREG_AUDIT(NBAD,NUNRESOLVED)
      IF (NBAD.NE.0.OR.NUNRESOLVED.NE.1) STOP 14
      BODY(3) = 1.D0
      KSTAR(3) = 13
C     Duplicate slots and bad types must be rejected, not repaired.
      NSCOUNT = 3
      NAMENS(3) = 30
      NSTYPE(3) = 13
      CALL PSRREG_AUDIT(NBAD,NUNRESOLVED)
      IF (NBAD.NE.1.OR.NSCOUNT.NE.3) STOP 15
      NSCOUNT = 2
      NSTYPE(2) = 14
      CALL PSRREG_AUDIT(NBAD,NUNRESOLVED)
      IF (NBAD.NE.1) STOP 16
      IF (SAVEDSEED.NE.IDUM1.OR.ANY(SAVEDRNG.NE.RNGWORDS)) STOP 17
      WRITE(6,*) 'PASS: transition retirement, RNG and restart audit'
      END
