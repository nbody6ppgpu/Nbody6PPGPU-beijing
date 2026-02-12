      SUBROUTINE KSRECT(IPAIR)
*
*
*       Rectification of KS orbit.
*       --------------------------
*
      Include 'kspars.h'
      INCLUDE 'common6.h'
      DATA NKSCNT,NKSBAS/0,1/
*
*       Include some diagnostic output.
      UPR2 = 0.0
      DO 1 K = 1,4
          UPR2 = UPR2 + UDOT(K,IPAIR)**2
    1 CONTINUE
*
*       Skip rectification for small eccentricity or large perturbation.
      I = N + IPAIR
      SEMI = -0.5*BODY(I)/H(IPAIR)
      ECC2 = (1.0-R(IPAIR)/SEMI)**2+TDOT2(IPAIR)**2/(SEMI*BODY(I))
      ECC = SQRT(ECC2)
      IF(ISNAN(SEMI).OR.ISNAN(ECC2).OR.ISNAN(ECC))THEN
          i1 = 2*ipair-1
          i2 = 2*ipair
      if(rank.eq.0)print*,' Warning KSRECT begin: ',
     &    't,ipair,n12,m12,semi,ecc2,ecc,kw,r,h,gamma,tdot2=',ttot,
     &    ipair,name(i1),name(i2),body(i1)*zmbar,body(i2)*zmbar, 
     &    semi,ecc2,ecc,kstar(n+ipair),r(ipair),h(ipair),gamma(ipair),
     &    tdot2(ipair)
      STOP
      END IF
*
      NKSCNT = NKSCNT + 1
      IF(MOD(NKSCNT,NKSBAS).EQ.0)THEN
      if(rank.eq.0)then
          ICM = N + IPAIR
          I1 = 2*IPAIR - 1
          I2 = 2*IPAIR
          RI = SQRT((X(1,ICM) - RDENS(1))**2 +
     &              (X(2,ICM) - RDENS(2))**2 +
     &              (X(3,ICM) - RDENS(3))**2)
          VI = SQRT(XDOT(1,ICM)**2+XDOT(2,ICM)**2+XDOT(3,ICM)**2)
      IF(ISNAN(RI).OR.ISNAN(VI).OR.ISNAN(SEMI).OR.ISNAN(ECC2).OR.
     &     ISNAN(ECC).OR.ECC2.GT.1.D0)THEN
      IF(ISNAN(RI).OR.ISNAN(VI))print*,' Warning RI,VI NaN: ',
     &   ' X,XDOT,RDENS 1-3=',X(1:3,ICM),XDOT(1:3,ICM),RDENS(1:3)
          EB = -0.5*BODY(I1)*BODY(I2)/SEMI
          PD = TWOPI*SEMI*SQRT(DABS(SEMI)/BODY(ICM))*TSTAR*365.24D6
          WRITE (6,17)  TTOT,NAME(I1),NAME(I2),
     &         NAME(ICM),KSTAR(I1),KSTAR(I2),KSTAR(ICM),
     &         IPAIR,DTAU(IPAIR),BODY(I1),BODY(I2),
     &         R(IPAIR),ECC,SEMI,EB,PD,H(IPAIR),GAMMA(IPAIR),
     &         STEP(ICM),LIST(1,I1),LIST(1,ICM),
     &         BODY(I1)*ZMBAR,BODY(I2)*ZMBAR,
     &         RADIUS(I1)*SU,RADIUS(I2)*SU,R(IPAIR)*SU,RI,VI
  17      FORMAT (/,' BEG KSRECT TIME[NB]',1P,E17.10,' NM1,2,S=',
     &         3I10,' KW1,2,S=',3I4,' IPAIR',I9,' DTAU',E11.3,
     &         ' M1,2[NB]',2E11.3,' R12[NB]',E11.3,
     &         ' e,a,eb[NB]=',2E12.4,E11.3,' P[d]=',E11.3,' H',E11.3,
     &         ' GAMMA',1P,E11.3,' STEP(ICM)',E11.3,' NPERT',I5,
     &         ' NB(ICM)',I5,' M1,2[*]',2E11.3,' RAD1,2,S[*]',3E11.3,
     &         ' RI,VI[NB]=',2E11.3)
      END IF
      end if
      END IF
*
      IF (ECC.LE.0.01) GO TO 50
      IF (GAMMA(IPAIR).GT.0.1) GO TO 50
*
*       Include Roche rectification for large ECC (additional to roche.f).
      IF (KSTAR(I).GT.10.AND.ECC.GT.0.2) THEN
          QP = SEMI*(1.0 - ECC)
          ICIRC = -1
          CALL TCIRC(QP,ECC,2*IPAIR-1,2*IPAIR,ICIRC,TC)
          if(rank.eq.0)
     &    WRITE (6,2)  NAME(2*IPAIR-1), KSTAR(I), ECC, QP, TC,
     &                 GAMMA(IPAIR)
    2     FORMAT (' RECTIFY K*    NM K* E QP TC G',I8,I4,F8.4,1P,3E11.3)
          KSTAR(I) = 0
*     ks MPI communication
*         call ksparmpi(K_store,K_int,K_KSTAR,I,0,KSTAR(I))
      END IF
*
      ICM = N + IPAIR
      I1 = 2*IPAIR - 1
      I2 = 2*IPAIR
      RI = SQRT((X(1,ICM) - RDENS(1))**2 +
     &              (X(2,ICM) - RDENS(2))**2 +
     &              (X(3,ICM) - RDENS(3))**2)
      VI = SQRT(XDOT(1,ICM)**2+XDOT(2,ICM)**2+XDOT(3,ICM)**2)
      HI = (2.0*UPR2 - BODY(ICM))/R(IPAIR)
      ERR = (HI - H(IPAIR))/HI
      ZMU = BODY(I1)*BODY(I2)/BODY(ICM)
      DB = ZMU*(HI - H(IPAIR))
      IF (ABS(DB).GT.1.0D-08) THEN
          SEMI = -0.5*BODY(ICM)/H(IPAIR)
          RA = R(IPAIR)/SEMI
          IF (SEMI.LT.0.0) RA = R(IPAIR)
          if(rank.eq.0)WRITE (16,616)  TIME+TOFF,NAME(I1),NAME(I2),
     &         NAME(ICM),KSTAR(I1),KSTAR(I2),KSTAR(ICM),
     &         IPAIR,DTAU(IPAIR),BODY(I1),BODY(I2),
     &         R(IPAIR),ECC,SEMI,EB,PD,H(IPAIR),GAMMA(IPAIR),
     &         STEP(ICM),LIST(1,I1),LIST(1,ICM),
     &         BODY(I1)*ZMBAR,BODY(I2)*ZMBAR,
     &         RADIUS(I1)*SU,RADIUS(I2)*SU,R(IPAIR)*SU,RI,VI,DB,ERR
 616      FORMAT (' KSRECT: Time[NB]',1P,E17.10,' NM1,2,S=',
     &         3I10,' KW1,2,S=',3I4,' IPAIR',I9,' DTAU',E11.3,
     &         ' M1,2[NB]',2E11.3,' R12[NB]',E11.3,
     &         ' e,a,eb[NB]=',2E12.4,E11.3,' P[d]=',E11.3,' H',E11.3,
     &         ' GAMMA',1P,E11.3,' STEP(ICM)',E11.3,' NPERT',I5,
     &         ' NB(ICM)',I5,' M1,2[*]',2E11.3,' RAD1,2,S[*]',3E11.3,
     &         ' RI,VI[NB]=',2E11.3,' DB,ERR=',2E11.3)
          CALL FLUSH(16)
      END IF
*
*       Initialize iteration counter for difficult case (SJA 10/97).
      ITER = 0
*
*       Form square regularized velocity for the explicit binding energy.
   10 UPR2 = 0.0
      DO 15 K = 1,4
          UPR2 = UPR2 + UDOT(K,IPAIR)**2
   15 CONTINUE
*
*       Form KS scaling factors from energy and angular momentum relation.
      A1 = 0.25D0*BODY(N+IPAIR)/UPR2
*       Solve for C1 from H = (2*U'*U'*C1**2 - M)/(U*U*C2**2) with C2 = 1/C1.
      A2 = A1**2 + 0.5D0*H(IPAIR)*R(IPAIR)/UPR2
*
*       Avoid negative round-off value on second try (NB! no change in CK).
      IF (ITER.EQ.2.AND.A2.LT.0.0) A2 = 0.0D0
*
*       Check for undefined case (circular orbit or eccentric anomaly = 90).
      IF (A2.GE.0.0D0) THEN
          IF (A1.LT.1.0) THEN
*       Choose square root sign from eccentric anomaly (e*cos(E) = 1 - R/a).
              C1 = SQRT(A1 + SQRT(A2))
          ELSE
              C1 = SQRT(A1 - SQRT(A2))
          END IF
          CK = 1.0
      ELSE
*       Adopt C1*C2 = CK for difficult case (Seppo's suggestion of 1991).
          C1 = 1.0
          CK = BODY(N+IPAIR)/SQRT(-8.0D0*H(IPAIR)*R(IPAIR)*UPR2)
          if(rank.eq.0)
     &    WRITE (6,20)  IPAIR, KSTAR(N+IPAIR), ECC, R(IPAIR), H(IPAIR),
     &                  GAMMA(IPAIR), UPR2, A2, CK-1.0
   20     FORMAT (' WARNING!    KSRECT    KS K* E R H G UPR2 A2 CK-1 ',
     &                                    I8,I4,F10.4,1P,6E11.3)
          ITER = ITER + 1
      END IF
*
*       Specify KS coordinate scaling from angular momentum conservation.
      C2 = CK/C1
*
*       Transform KS variables to yield the prescribed elements.
      R(IPAIR) = 0.0D0
*     UPR2 = 0.0D0
      TD2 = 0.0D0
      DO 25 K = 1,4
          U(K,IPAIR) = C2*U(K,IPAIR)
          UDOT(K,IPAIR) = C1*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
          TD2 = TD2 + U(K,IPAIR)*UDOT(K,IPAIR)
*         UPR2 = UPR2 + UDOT(K,IPAIR)**2
   25 CONTINUE
      TDOT2(IPAIR) = 2.0*TD2
*
*       Include diagnostic output.
*     HI = (2.0*UPR2 - BODY(N+IPAIR))/R(IPAIR)
*     if(rank.eq.0)
*    &WRITE (16,30)  IPAIR, IPHASE, KSTAR(N+IPAIR), R(IPAIR),
*    &               GAMMA(IPAIR), (HI-H(IPAIR))/HI
*  30 FORMAT (' KSRECT:    KS IPH K* R G DH/H ',3I4,1P,3E11.3)
*     CALL FLUSH(16)
*
*       Improve solution by second iteration in case of CK procedure.
      ITER = ITER + 1
      IF (ITER.EQ.2) GO TO 10
*
*     IF(MOD(NKSCNT,NKSBAS).EQ.0)THEN
*     NKSCNT = 0
*     if(rank.eq.0)then
*         print*,' Output KSRECT begin: nksrec,',
*    &  'ipair,semi,ecc2,ecc,kw,r,h,gamma=',nksrec,ipair,semi,ecc2,ecc,
*    &    kstar(n+ipair),r(ipair),h(ipair),gamma(ipair)
*         ICM = N + IPAIR
*         I1 = 2*IPAIR - 1
*         I2 = 2*IPAIR
*         RI = SQRT((X(1,ICM) - RDENS(1))**2 +
*    &              (X(2,ICM) - RDENS(2))**2 +
*    &              (X(3,ICM) - RDENS(3))**2)
*         VI = SQRT(XDOT(1,ICM)**2+XDOT(2,ICM)**2+XDOT(3,ICM)**2)
*         EB = -0.5*BODY(I1)*BODY(I2)/SEMI
*         PD = TWOPI*SEMI*SQRT(DABS(SEMI)/BODY(ICM))*TSTAR*365.24D6
*         WRITE (6,16)  TIME+TOFF,NAME(I1),NAME(I2),
*    &         NAME(ICM),KSTAR(I1),KSTAR(I2),KSTAR(ICM),
*    &         IPAIR,DTAU(IPAIR),BODY(I1),BODY(I2),
*    &         R(IPAIR),ECC,SEMI,EB,PD,H(IPAIR),GAMMA(IPAIR),
*    &         STEP(ICM),LIST(1,I1),LIST(1,ICM),
*    &         BODY(I1)*ZMBAR,BODY(I2)*ZMBAR,
*    &         RADIUS(I1)*SU,RADIUS(I2)*SU,R(IPAIR)*SU,RI,VI,ITER
* 16      FORMAT (/,' END KSRECT TIME[NB]',1P,E17.10,' NM1,2,S=',
*    &         3I10,' KW1,2,S=',3I4,' IPAIR',I9,' DTAU',E11.3,
*    &         ' M1,2[NB]',2E11.3,' R12[NB]',E11.3,
*    &         ' e,a,eb[NB]=',2E12.4,E11.3,' P[d]=',E11.3,' H',E11.3,
*    &         ' GAMMA',1P,E11.3,' STEP(ICM)',E11.3,' NPERT',I5,
*    &         ' NB(ICM)',I5,' M1,2[*]',2E11.3,' RAD1,2,S[*]',3E11.3,
*    &         ' RI,VI[NB]=',2E11.3,' ITER=',I5)
*     end if
*     END IF

   50 RETURN
*
      END
