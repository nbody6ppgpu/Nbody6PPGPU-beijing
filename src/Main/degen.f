      SUBROUTINE DEGEN(I1,I2,ICASE)
*
*
*       Binary output for degenerate stars.
*       -----------------------------------
*
      INCLUDE 'common6.h'
      REAL*8 EB(KMAX),SEMI(KMAX),ECC(KMAX),RCM(KMAX),ECM(KMAX)
      CHARACTER*8  WHICH1
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
*
*
*       Skip output for start and end of merger.
      IF (IPHASE.EQ.6.OR.IPHASE.EQ.7) GO TO 50
      IF (ICASE.EQ.7) GO TO 40
*
*       See whether KS binaries contain any degenerate stars.
      NB = 0
      DO 1 IPAIR = I1,I2
          J2 = 2*IPAIR
          J1 = J2 - 1
          IF (KSTAR(J1).GT.9.OR.KSTAR(J2).GT.9) THEN
              NB = NB + 1
          END IF
    1 CONTINUE
*
      IF (NB.GT.0.AND.FIRST.AND.KSTART.EQ.1) THEN
          FIRST = .FALSE.
*
*       Print cluster scaling parameters at start of the run.
          if(rank.eq.0)then
             WRITE (4,2)  RBAR, BODYM*ZMBAR, BODY1*ZMBAR, TSCALE,
     &            NBIN0, NZERO
 2           FORMAT (/,6X,1P, 'MODEL:    RBAR =',E26.17,'  <M>[M*] =',
     &            E26.17,'  M1[M*] =',E26.17,'  TSCALE =',E26.17,0P,
     &            '  NB0 =',I12,'  N0 =',I12,//)
*     
             WRITE (4,3)
 3           FORMAT ('       ICASE                 Time[Myr]',
     &         '                  SEMI[AU]                       ECC',
     &         '                   PERI/RS                   P[days]',
     &         '                    RI[PC]                 M(I1)[M*]',
     &         '                 M(I2)[M*]      K*(I1)      K*(I2)',
     &         '     K*(ICM)    NAME(I1)    NAME(I2)',/)
          end if
      END IF
*
*       Form binding energy and central distance for degenerate stars.
      JPAIR = 0
      IPRINT = 0
      DO 20 IPAIR = I1,I2
          J2 = 2*IPAIR
          J1 = J2 - 1
          IF (KSTAR(J1).LT.10.AND.KSTAR(J2).LT.10) GO TO 20
          JPAIR = JPAIR + 1
          ICM = N + IPAIR
          CALL JPRED(ICM,TIME,TIME)
*       Avoid division by zero for merged or synchronous ghost binary.
          IF (BODY(J1).GT.0.0) THEN
              EB(JPAIR) = BODY(J1)*BODY(J2)*H(IPAIR)/
     &                                             (BODY(J1) + BODY(J2))
              SEMI(JPAIR) = -0.5d0*BODY(ICM)/H(IPAIR)
              ECC2 = (1.d0 - R(IPAIR)/SEMI(JPAIR))**2 +
     &                       TDOT2(IPAIR)**2/(BODY(ICM)*SEMI(JPAIR))
              ECC(JPAIR) = SQRT(ECC2)
*       Set zero eccentricity after common envelope stage (still large R).
              IF (R(IPAIR).GT.2.0*SEMI(JPAIR)) THEN
                  ECC(JPAIR) = 0.d0
              END IF
C              EB(JPAIR) = MAX(EB(JPAIR),-9.99999d0)
          ELSE
              EB(JPAIR) = 0.d0
              SEMI(JPAIR) = R(IPAIR)
              ECC(JPAIR) = 0.d0
          END IF
          RCM(JPAIR) = SQRT((X(1,ICM) - RDENS(1))**2 +
     &                      (X(2,ICM) - RDENS(2))**2 +
     &                      (X(3,ICM) - RDENS(3))**2)
*       Obtain binding energy (per unit mass) of c.m. motion.
          VJ2 = XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + 
     &                           XDOT(3,ICM)**2
*         POTJ = 0.0
*         DO 5 J = IFIRST,NTOT
*             IF (J.EQ.ICM) GO TO 5
*             RIJ2 = (X(1,ICM) - X(1,J))**2 + (X(2,ICM) - X(2,J))**2 +
*    &                                        (X(3,ICM) - X(3,J))**2 
*             POTJ = POTJ + BODY(J)/SQRT(RIJ2)
*   5     CONTINUE
*         IF (TIME.EQ.TADJ) THEN
*             POTJ = -PHI(ICM)
*         ELSE
              POTJ = 0.d0
*         END IF
          ECM(JPAIR) = 0.5d0*VJ2 - POTJ
*       Check for external tidal field (NB! already included on GRAPE).
*         IF (KZ(14).GT.0) THEN
*             CALL XTRNLV(ICM,ICM)
*             ECM(JPAIR) = ECM(JPAIR) + HT/(BODY(ICM) + 1.0D-20)
*         END IF
          TPHYS = TTOT*TSTAR
          TK = SEMI(JPAIR)*SQRT(ABS(SEMI(JPAIR))/(BODY(ICM) + 1.0d-20))
          TK = DAYS*TK
C          TK = MIN(TK,999999.9d0)
          RBIG = MAX(RADIUS(J1),RADIUS(J2))
          RATIO = SEMI(JPAIR)*(1.d0 - ECC(JPAIR))/RBIG
C          RATIO = MIN(RATIO,99.9d0)
          IF (SEMI(JPAIR).LT.0.0.AND.RATIO.GT.5.0) GO TO 20
          SEMI(JPAIR) = SEMI(JPAIR)*RBAR*AU
          IF (IPRINT.EQ.0.AND.ICASE.EQ.0) THEN
              if(rank.eq.0) WRITE (4,*)
          END IF
          if(rank.eq.0)then
          WRITE (4,*)  ICASE, TPHYS, SEMI(JPAIR), ECC(JPAIR), RATIO,
     &            TK, RCM(JPAIR)*RBAR, BODY(J1)*ZMBAR,
     &            BODY(J2)*ZMBAR,KSTAR(J1), KSTAR(J2),
     &            KSTAR(ICM),NAME(J1), NAME(J2)
C   10     FORMAT (I2,F8.1,F8.2,F7.3,F6.1,F9.1,F6.2,2F5.1,3I4,2I6)
          end if
          IPRINT = IPRINT + 1
   20 CONTINUE
*
*       Close file at end of main output.
      IF (IPRINT.GT.0.AND.ICASE.EQ.0) THEN
          CALL FLUSH(4)
      END IF
*
*       Search for neutron stars at main output. (seem never used)
      IF (ICASE.EQ.0) THEN
          DO 30 J = 1,N
              IF (KSTAR(J).EQ.13) THEN
                  IF (J.LT.IFIRST) THEN
                      WHICH1 = ' BINARY '
                      I = KVEC(J) + N
                      CALL JPRED(I,TIME,TIME)
                      VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + 
     &                                     XDOT(3,I)**2
                  ELSE
                      WHICH1 = ' SINGLE '
                      VI2 = XDOT(1,J)**2 + XDOT(2,J)**2 + 
     &                                     XDOT(3,J)**2
                  END IF
                  if(rank.eq.0)then
                  WRITE (33,25)  WHICH1, J, NAME(J), IFIRST, KSTAR(J),
     &                           TPHYS, SQRT(VI2)*VSTAR
   25             FORMAT (1X,A8,'NS','  I NAME IFIRST K* Time[Myr] ',
     &                 'VI[km/s] ', 3I12,I4,1P,E25.16,E16.7,0P)
                  CALL FLUSH(33)
                  end if
              ELSE IF (KSTAR(J).GT.13) THEN
                  IF (J.LT.IFIRST) THEN
                      JCM = N + KVEC(J)
                      IF (NAME(JCM).LT.0) GO TO 30
                  END IF
                  VI2 = XDOT(1,J)**2 + XDOT(2,J)**2 + 
     &                                 XDOT(3,J)**2
                  VI = SQRT(VI2)
                  WHICH1 = ' AIC    '
                  IF (KSTAR(J).EQ.14) WHICH1 = ' BH     '
                  IF (KSTAR(J).EQ.15) THEN
                     CALL JPRED(J,TIME,TIME)
*       Ensure that errant massless remnant will escape.
                      RI = SQRT(X(1,J)**2 + X(2,J)**2 + X(3,J)**2)
                      DO 26 L = 1,3
                          X0(L,J) = 1000.0*RSCALE*X(L,J)/RI
                          X(L,J) = X0(L,J)
                          X0DOT(L,J) = SQRT(0.004*ZMASS/RSCALE)*
     &                                 XDOT(L,J)/VI
                          XDOT(L,J) = X0DOT(L,J)
   26                 CONTINUE
                  ELSE
                      VI = VI*VSTAR
                      if(rank.eq.0)then
                      WRITE (34,28) WHICH1,J,NAME(J),KSTAR(J),TPHYS,VI
   28                 FORMAT (1X,A8,'I NAME K* Time[Myr] VI[km/s] ',
     &                     2I12,I4,1P,E25.16,E16.7,0P)
                      CALL FLUSH(34)
                      end if
                  END IF
              END IF
   30     CONTINUE
      END IF
*
*       Print single neutron stars at escape time.
   40 IF (ICASE.EQ.7) THEN
          J = I1
          VI2 = XDOT(1,J)**2 + XDOT(2,J)**2 + XDOT(3,J)**2
          VI = SQRT(VI2)*VSTAR
          WHICH1 = ' ESCAPE '
          if(rank.eq.0)then
          WRITE (33,25)  WHICH1, J, NAME(J), IFIRST, KSTAR(J), TPHYS, VI
          end if
      END IF
*
*       Include counter for doubly degenerate binary.
      IF (ICASE.EQ.3.OR.ICASE.EQ.4) THEN
          IPAIR = I1
          NPAIR = N + IPAIR
          J1 = 2*IPAIR - 1
          J2 = J1 + 1
          IF (KSTAR(J1).GE.10.AND.KSTAR(J2).GE.10) THEN
          RII=DSQRT((X(1,NPAIR) - RDENS(1))**2 +
     &              (X(2,NPAIR) - RDENS(2))**2 +
     &              (X(3,NPAIR) - RDENS(3))**2)
          VII=DSQRT(XDOT(1,NPAIR)**2+XDOT(2,NPAIR)**2+XDOT(3,NPAIR)**2)
              if(rank.eq.0)then
              WRITE (6,48) TTOT,NAME(J1),NAME(J2),NAME(NPAIR),KSTAR(J1),
     &                     KSTAR(J2),KSTAR(NPAIR),IPAIR,DTAU(IPAIR),
     &                     BODY(J1),BODY(J2),R(IPAIR),ECC(IPAIR),
     &                  SEMI(IPAIR),EB(IPAIR),TK,H(IPAIR),GAMMA(IPAIR),
     &                     STEP(NPAIR),LIST(1,J1),LIST(1,NPAIR),
     &                     BODY(J1)*ZMBAR,BODY(J2)*ZMBAR,
     &                     RADIUS(J1)*SU,RADIUS(J2)*SU,R(IPAIR)*SU,
     &                     RII,VII
  48      FORMAT (/,' NEW DOUBLE DEGEN TIME[NB]',1P,E15.7,' NM1,2,S=',
     &         3I10,' KW1,2,S=',3I4,' IPAIR',I9,' DTAU',E13.5,
     &         ' M1,2[NB]',2E13.5,' R12[NB]',E13.5,
     &         ' e,a,eb[NB]=',3E13.5,' P[d]=',E13.5,' H',E13.5,
     &         ' GAMMA',E13.5,' STEP(ICM)',E13.5,' NPERT',I5,
     &         ' NB(ICM)',I5,' M1,2[*]',2E13.5,' RAD1,2,S[*]',3E13.5,
     &         ' RI,VI[NB]=',2E13.5)
              end if
          END IF
      END IF
*
   50 RETURN
*
      END
