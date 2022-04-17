      SUBROUTINE BODIES
*
*
*       Output of single bodies or binaries.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8 RCM(3),VCM(3)
*
*       Optional search for soft binaries (frequency NFIX with KZ(6) = 4).
      IF (KZ(6).EQ.0) GO TO 70
      IF (KZ(6).EQ.2) GO TO 50
      IF (KZ(6).EQ.3.AND.NPRINT.NE.1) GO TO 50
      IF (KZ(6).EQ.4.AND.NPRINT.NE.1) GO TO 70
      SIMAX = 0.01*TCR
*
      DO 40 I = IFIRST,NTOT
          IF (STEP(I).GT.SIMAX) GO TO 40
          JMIN = 0
          RJMIN2 = RSCALE**2
          NNB = LIST(1,I)
          DO 30 L = 1,NNB
              J = LIST(L+1,I)
              IF (STEP(J).GT.SIMAX) GO TO 30
              A1 = X(1,I) - X(1,J)
              A2 = X(2,I) - X(2,J)
              A3 = X(3,I) - X(3,J)
              RIJ2 = A1**2 + A2**2 + A3**2
              IF (RIJ2.LT.RJMIN2) THEN
                  RJMIN2 = RIJ2
                  JMIN = J
              END IF
   30     CONTINUE
          IF (JMIN.LE.I) GO TO 40
          RIJMIN = SQRT(RJMIN2)
          VR2 = (XDOT(1,I) - XDOT(1,JMIN))**2 +
     &          (XDOT(2,I) - XDOT(2,JMIN))**2 +
     &          (XDOT(3,I) - XDOT(3,JMIN))**2
          EREL = 0.5*VR2 - (BODY(I) + BODY(JMIN))/RIJMIN
*       Only print significant binaries.
          IF (EREL.GT.-0.1*ECLOSE) GO TO 40
          SEMI = -0.5*(BODY(I) + BODY(JMIN))/EREL
          ZN = SQRT((BODY(I) + BODY(JMIN))/SEMI**3)
          RDOT = (X(1,I) - X(1,JMIN))*(XDOT(1,I) - XDOT(1,JMIN)) +
     &           (X(2,I) - X(2,JMIN))*(XDOT(2,I) - XDOT(2,JMIN)) +
     &           (X(3,I) - X(3,JMIN))*(XDOT(3,I) - XDOT(3,JMIN))
          ECC2 = (1.0 - RIJMIN/SEMI)**2 +
     &                             RDOT**2/(SEMI*(BODY(I) + BODY(JMIN)))
          ECC = SQRT(ECC2)
          RCM(1:3) = BODY(I)*X(1:3,I) + BODY(JMIN)*X(1:3,JMIN)
          VCM(1:3) = BODY(I)*XDOT(1:3,I) + BODY(JMIN)*XDOT(1:3,JMIN)
          RCM2 = RCM(1)**2+RCM(2)**2+RCM(3)**2
          VCM2 = VCM(1)**2+VCM(2)**2+VCM(3)**2
          if(rank.eq.0)
     &    WRITE (6,36)  TIME, NAME(I), NAME(JMIN), BODY(I), BODY(JMIN),
     &  EREL, SEMI, ZN, RIJMIN, ECC, LIST(1,I), SQRT(RCM2), SQRT(VCM2)
   36     FORMAT (' significant binaries ',1P,E18.5,2I9,
     &                  7E12.4,I5, 2E12.4)
   40 CONTINUE
*
*       Output of regularized binaries (frequency NFIX with KZ(6) = 4).
   50 DO 60 JPAIR = 1,NPAIRS
          IF (H(JPAIR).GE.0.0) GO TO 60
          I = 2*JPAIR - 1
          ICM = N + JPAIR
          JMIN = I + 1
          IF (BODY(I).LE.0.0D0) GO TO 60
          EB = BODY(I)*BODY(JMIN)*H(JPAIR)/BODY(ICM)
          SEMI = -0.5*(BODY(I) + BODY(JMIN))/H(JPAIR)
          ZN = SQRT((BODY(I) + BODY(JMIN))/SEMI**3)
          RP = U(1,JPAIR)**2 + U(2,JPAIR)**2 + U(3,JPAIR)**2 +
     &                                         U(4,JPAIR)**2
          ECC2 = (1.0 - RP/SEMI)**2 +
     &                     TDOT2(JPAIR)**2/(SEMI*(BODY(I) + BODY(JMIN)))
          ECC = SQRT(ECC2)
          RI = SQRT(X(1,ICM)**2 + X(2,ICM)**2 + X(3,ICM)**2)
          VI = SQRT(XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + XDOT(3,ICM)**2)
          if(rank.eq.0)
     &    WRITE (6,35)  TIME, NAME(I), NAME(JMIN), NAME(ICM), BODY(I),
     &           BODY(JMIN), EB, SEMI, ZN, RP, ECC, LIST(1,N+JPAIR),
     &           GAMMA(JPAIR), RI, VI, KSLOW(JPAIR)
   35     FORMAT (' regularized binaries ',1P,E18.5,
     &                  3I9,7E12.4,I5,3E12.4,I9)
   60 CONTINUE
*
   70 RETURN
*
      END
