C     Test-only setup of a restart with explicitly specified pulsar state.
C     START runs with KZ(29)=0 so this is NOT support for fresh NS input.
C     The production executable subsequently reads the ordinary dump and
C     integrates an accreting binary with the real ROCHE/HRDIAG routines.
      PROGRAM SEED_PSR_RESTART
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      INCLUDE 'pulsar_ctrl.h'
      INTEGER J, KZ50SAVE
      LOGICAL REGISTD
      NAMELIST /INNBODY6/ KSTART,TCOMP,TCRTP0,
     &                   isernb,iserreg,iserks
      rank = 0
      icore = 1
      READ(5,NML=INNBODY6)
      TCRITP = TCRTP0
      CPU = TCOMP
      CALL START
      KZ(29) = 2
      CALL READPULSAR
      KZ50SAVE = KZ(50)
      KZ(50) = 0
      DO J = 1, N
         IF (KSTAR(J).EQ.13) THEN
            CALL PSRREG(NAME(J),BODY(J)*ZMBAR,13,0.D0,
     &                  .FALSE.,REGISTD)
C           Explicit current state, not inferred historical NS birth.
            PERIODNS(NSCOUNT) = 0.05D0
            BMAGNS(NSCOUNT) = 1.D8
            PDOTNS(NSCOUNT) = 1.D-15
            AGENSX(NSCOUNT) = 0.D0
         ENDIF
      END DO
      KZ(50) = KZ50SAVE
      CALL ADJUST
      CALL MYDUMP(1,1)
      WRITE(6,*) 'SEEDED RESTART: NSCOUNT=', NSCOUNT
      END
