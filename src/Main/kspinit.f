      subroutine kspinit
*     
*     
*     Generate perturber list for primordial binaries during initialization
*     ------------------------------------------------------

      include 'common6.h'
      INCLUDE 'omp_lib.h'
      REAL*8  UI(4),VI(4)
*     
      TIME0 = TIME
!!$omp parallel do 
!!$omp& private(ipair,icm,i1,i2,eb,semi,tk,vi,ui,tp,imod,K)
      DO IPAIR = 1, NPAIRS
         ICM = IPAIR + N
         I2 = 2*IPAIR
         I1 = I2 - 1

*     Form perturber list.
         CALL KSLIST(IPAIR)
*     
*     Transform any unperturbed hard binary to apocentre and set time-step.
         EB = H(IPAIR)*BODY(I1)*BODY(I2)/BODY(ICM)
         SEMI = -0.5*BODY(ICM)/H(IPAIR)
         IF (LIST(1,I1).EQ.0.AND.EB.LT.EBH) THEN
            TK = TWOPI*SEMI*SQRT(SEMI/BODY(ICM))
*     Note TIME is not commensurate after KSPERI (cf. CHTERM & STEPS).
            DO 55 K = 1,4
               UI(K) = U(K,IPAIR)
               VI(K) = UDOT(K,IPAIR)
 55         CONTINUE
*     Determine pericentre time (TP < 0 if TDOT2 < 0) and add TK/2.
            CALL TPERI(SEMI,UI,VI,BODY(ICM),TP)
*     Note: apocentre to apocentre gives almost zero step.
            STEP(I1) = 0.5*MIN(TK,STEP(ICM)) - TP
*     Transform KS variables to peri and by pi/2 to apocentre (skip apo).
            IF (ABS(TDOT2(IPAIR)).GT.1.0E-12.OR.R(IPAIR).LT.SEMI) THEN
               TIME = TIME0
               CALL KSPERI(IPAIR)
               CALL KSAPO(IPAIR)
*     Reset TIME to quantized value (small > 0 or < 0 possible initially).
            ELSE IF (TDOT2(IPAIR).GT.0.0) THEN
               TDOT2(IPAIR) = -1.0E-20
            END IF
         END IF
*     
*     Estimate an appropriate KS slow-down index for G < GMIN.
         IMOD = 1
         IF (LIST(1,I1).EQ.0.AND.SEMI.GT.0.0) THEN
            TK = TWOPI*SEMI*SQRT(SEMI/BODY(ICM))
            IF (KZ(26).GT.0.AND.STEP(ICM).GT.TK) THEN
               IMOD = 1 + INT(LOG(STEP(ICM)/TK)/0.69)
               IMOD = MIN(IMOD,5)
            END IF
         END IF
*
*     Obtain polynomials for perturbed KS motion (standard case & merger).
         TIME = TIME0
         CALL KSPOLY(IPAIR,IMOD)
*
         LIST(2,I2) = -1
*
      END DO
!!$omp end parallel do      
*
*     Check optional output for new KS.
      DO IPAIR = 1, NPAIRS
         ICM = IPAIR + N
         I2 = 2*IPAIR
         I1 = I2 - 1
*
*       Obtain apocentre distance.
      SEMI = -0.5*BODY(ICM)/H(IPAIR)
      ECC2 = (1.0-R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(ICM)*SEMI)
      RAP = SEMI*(1.0 + SQRT(ECC2))
*
      IF (KZ(10).GT.0) THEN
         RIJ2 = (X(1,ICM) - RDENS(1))**2 + (X(2,ICM) - RDENS(2))**2 +
     &          (X(3,ICM) - RDENS(3))**2
         RI = DSQRT(RIJ2)
         VI2 = XDOT(1,ICM)**2+XDOT(2,ICM)**2+XDOT(3,ICM)**2
         PD = TWOPI*SEMI*SQRT(DABS(SEMI)/BODY(ICM))*TSTAR*365.24D6
         if(rank.eq.0)
     &    WRITE (6,60)  TIME+TOFF,NAME(I1),NAME(I2),
     &         NAME(ICM),KSTAR(I1),KSTAR(I2),KSTAR(ICM),
     &         IPAIR,DTAU(IPAIR),BODY(I1),BODY(I2),
     &         R(IPAIR),SQRT(ECC2),SEMI,EB,PD,H(IPAIR),GAMMA(IPAIR),
     &         STEP(ICM),LIST(1,I1),LIST(1,ICM),
     &         BODY(I1)*ZMBAR,BODY(I2)*ZMBAR,
     &         RADIUS(I1)*SU,RADIUS(I2)*SU,R(IPAIR)*SU,RI,DSQRT(VI2)
         call flush(6)
      END IF
      END DO
*
  60      FORMAT (/,' NEW KSREG   TIME[NB]',1P,E17.10,' NM1,2,S=',
     &         3I10,' KW1,2,S=',3I4,' IPAIR',I9,' DTAU',E10.2,
     &         ' M1,2[NB]',2E10.2,' R12[NB]',E10.2,
     &         ' e,a,eb[NB]=',2E12.4,E10.2,' P[d]=',E10.2,' H',E10.2,
     &         ' GAMMA',1P,E10.2,' STEP(ICM)',E10.2,' NPERT',I5,
     &         ' NB(ICM)',I5,' M1,2[*]',2E10.2,' RAD1,2,S[*]',3E10.2,
     &         ' RI,VI[NB]=',2E10.2)
*
      TIME = TIME0
*
      RETURN
*
      END
