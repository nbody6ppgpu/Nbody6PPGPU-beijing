      SUBROUTINE KSRES_OP(IPAIR,I1,I2,Q,RDOT,R2DOT,R3DOT,
     &                                R4DOT,R5DOT,KCASE)
*
*
*       Coordinates & velocities of KS pair for output
*       ------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      REAL*8  Q(3),RDOT(3),R2DOT(3),R3DOT(3),R4DOT(3),R5DOT(3)
      REAL*8  UI(4),V(4),A1(3,4),VDOT(4),V3(4),V4(4),V5(4)
      REAL*8  XL1(3,4),XL2(3,4),XL3(3,4)
*
*       Note: KDUM = 1 for prediction of U & UDOT and = 2 for termination.
      I1 = 2*IPAIR - 1
      I2 = 2*IPAIR
*
*       Ensure appropriate prediction (second pair in binary collision).
      IF (LIST(1,I1).EQ.0.OR.T0(I1).EQ.TIME.OR.TIME.EQ.0.D0) THEN
*       Set current values of regularized coordinates & velocities.
          DO 1 K = 1,4
              UI(K) = U0(K,IPAIR)
              V(K) = UDOT(K,IPAIR)
              VDOT(K) = 2.D0*FU(K,IPAIR)
              V3(K) = 6.D0*FUDOT(K,IPAIR)
              V4(K) = FUDOT2(K,IPAIR)
              V5(K) = FUDOT3(K,IPAIR)
    1     CONTINUE
*         DTU = 0.D0
*       Copy binding energy for routine ENERGY.
          HT = H(IPAIR)
          GO TO 4
      END IF
*
*       Predict U & UDOT to order FUDOT3 using third-order time inversion.
      A2 = 1.0/R(IPAIR)
      A3 = A2*(TIME - T0(I1))
*
*       See whether the time interval should be modified by KSLOW procedure.
      IF (KSLOW(IPAIR).GT.1) THEN
          IMOD = KSLOW(IPAIR)
          A3 = A3/FLOAT(ISLOW(IMOD))
      END IF
*
*       Expand regularized interval to third order.
      A4 = 3.0D0*TDOT2(IPAIR)**2*A2 - TDOT3(IPAIR)
      DTU = ((ONE6*A4*A3 - 0.5D0*TDOT2(IPAIR))*A2*A3 + 1.0)*A3
*       Apply safety test near small pericentre or for unperturbed motion.
      IF (DTU.GT.DTAU(IPAIR)) DTU = 0.8*DTAU(IPAIR)
      IF (DTU.LT.0.0) DTU = 0.0
*
      DTU1 = 0.2D0*DTU
      DTU2 = DTU/24.0D0
      DO 3 K = 1,4
          UI(K) = ((((FUDOT3(K,IPAIR)*DTU1 + FUDOT2(K,IPAIR))*DTU2 +
     &                          FUDOT(K,IPAIR))*DTU + FU(K,IPAIR))*DTU +
     &                                  UDOT(K,IPAIR))*DTU + U0(K,IPAIR)
          V(K) = (((FUDOT3(K,IPAIR)*DTU2 + ONE6*FUDOT2(K,IPAIR))*DTU +
     &              3.0D0*FUDOT(K,IPAIR))*DTU + 2.0D0*FU(K,IPAIR))*DTU +
     &                                                     UDOT(K,IPAIR)
          VDOT(K) = ((FUDOT3(K,IPAIR)*DTU2 + ONE6*FUDOT2(K,IPAIR))*DTU +
     &              3.0D0*FUDOT(K,IPAIR))*DTU + 2.0D0*FU(K,IPAIR)
          V3(K) = (ONE6*FUDOT3(K,IPAIR)*DTU+FUDOT2(K,IPAIR)/2.D0)*DTU +
     &              6.0D0*FUDOT(K,IPAIR)
          V4(K) = FUDOT3(K,IPAIR)*DTU/2.D0 + FUDOT2(K,IPAIR)
          V5(K) = FUDOT3(K,IPAIR)
    3 CONTINUE
*
*       Predict current binding energy per unit mass for routine ADJUST.
      HT = (((HDOT4(IPAIR)*DTU2 + ONE6*HDOT3(IPAIR))*DTU +
     &             0.5D0*HDOT2(IPAIR))*DTU + HDOT(IPAIR))*DTU + H(IPAIR)
*
*       Form current transformation matrix and two-body separation.
    4 CALL MATRIX(UI,XL1)
      CALL MATRIX(V,XL2)
      CALL MATRIX(VDOT,XL3)
*
      RI = UI(1)**2 + UI(2)**2 + UI(3)**2 + UI(4)**2
      RINV = 2.D0/RI
      RINV2 = RINV/RI
      RINV3 = RINV2/RI
      RINV4 = RINV3/RI
      RINV5 = RINV4/RI
*
*       Obtain relative coordinates & velocities.
      DO 7 J = 1,3
          Q(J) = 0.0D0
          RDOT(J) = 0.0D0
          R2DOT(J) = 0.0D0
          R3DOT(J) = 0.0D0
          R4DOT(J) = 0.0D0
          R5DOT(J) = 0.0D0
          DO 6 K = 1,4
              Q(J) = Q(J) + XL1(J,K)*UI(K)
              RDOT(J) = RDOT(J)+RINV*XL1(J,K)*V(K)
              R2DOT(J) = R2DOT(J)+RINV2*(XL2(J,K)*V(K)+XL1(J,K)*VDOT(K))
              R3DOT(J) = R3DOT(J)+RINV3*(3.D0*XL2(J,K)*VDOT(K)+
     &                                  XL1(J,K)*V3(K))
              R4DOT(J) = R4DOT(J)+RINV4*(3.D0*XL3(J,K)*VDOT(K)+
     &                             4.D0*XL2(J,K)*V3(K)+XL1(J,K)*V4(K))
              R5DOT(J) = R5DOT(J)+RINV5*(10.D0*XL3(J,K)*V3(K)+
     &                             5.D0*XL2(J,K)*V4(K)+XL1(J,K)*V5(K))
    6     CONTINUE
    7 CONTINUE
*
*     if(rank.eq.0)
*    &WRITE(555,555)IPAIR,I1,I2,TIME,T0(I1),Q(1),RDOT(1),R2DOT(1),
*    &      R3DOT(1),R4DOT(1),R5DOT(1),XL1(1,1),XL2(1,1),XL3(1,1),
*    &      UI(1),V(1),VDOT(1),V3(1),V4(1),V5(1),DTU
*555  FORMAT(1X,' IP,1,2,T,T0=',3I8,1P,2D13.5,' Q R1-R5=',6D13.5,
*    &    ' XL1-3=',3D13.5,' U,V,VDOT,V3-5=',6D13.5,' DTU=',D13.5)
      IF(KCASE.EQ.1) THEN
         I = N + IPAIR

*     Set coordinates & velocities (XDOT for KDUM = 1) of KS components.
         DO 15 K = 1,3
            X(K,I1) = X(K,I) + BODY(I2)*Q(K)/BODY(I)
            X(K,I2) = X(K,I) - BODY(I1)*Q(K)/BODY(I)
            XDOT(K,I1) = XDOT(K,I) + BODY(I2)*RDOT(K)/BODY(I)
            XDOT(K,I2) = XDOT(K,I) - BODY(I1)*RDOT(K)/BODY(I)
 15      CONTINUE
      END IF
*
      RETURN
*
      END
