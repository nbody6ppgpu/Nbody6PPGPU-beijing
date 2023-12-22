       SUBROUTINE ACCRETION(I)
*
*
*       Accretion of stars onto central massive black hole.
*       Based on "star_destr_ext.c"
*       -------------------------
*       by Taras Panamarev
*
      INCLUDE 'common6.h'
      REAL*8 R_TIDAL, E_kin, E_pot_bh, E_pot, E_corr_tmp,
     &       EPS, EPS_bh, EPS_bh2, RSB, RKB2, RKS2, R2, RR, NNB,
     &       X_IJ, Y_IJ, Z_IJ, VX_IJ, VY_IJ, VZ_IJ, V, TMP,
     &       A_EXT_I(3), ADOT_EXT_I(3), TEMP1, EPS_SS, EPS_SS2,
     &       DTR,DTR2HALF,DT_TMP,DT2HALF, DT3OVER6,RKS,XP(3),EPOT1
     &       ,t_diss_on, raccr, semi
*       Declarations for call of sev... routines
      real*4 b_rs1,b_rs2,b_l1,b_l2,b_te1,b_te2,b_rc1,b_rc2,b_m1,b_m2,
     &       b_mc1,b_mc2,b_ecc, b_a, b_p,
     &       s_rs, s_l, s_te, s_rc, s_m, s_mc
*       Declarations for call of ksres_op routine
      REAL*8  Q(3),RDOT(3),R2DOT(3),R3DOT(3),R4DOT(3),R5DOT(3)
*
      INTEGER ii, L, KW
*      LOGICAL is_active
*
*     Time when accretion gets 'switched on' 
      t_diss_on = 0.125 
*
      kw=-1
*     Update total time:
      TTOT = TIME + TOFF
*     
      if (ttot.lt.t_diss_on) goto 10     
*
      EPS_bh2 = EPS_bh**2
*
       E_pot = 0.0D0
       E_kin = 0.0D0
       E_pot_bh = 0.0D0
*
       R_TIDAL = 0.22*RACCR
*
      iaccr = 0
*	   
      K = I
*    Change x to x0, because x is predicted value
      RSB = SQRT(X0(1,K)**2+X0(2,K)**2+X0(3,K)**2)
*
      IF (RSB.LT.R_TIDAL.and.RSB.ne.0.0) THEN
*
            IF (K.GT.N) THEN 
*             binary is being accreted:
             KSPAIR = K - N
*
             call KSRES_OP(KSPAIR,J1,J2,Q,RDOT,R2DOT,R3DOT,R4DOT,R5DOT)

             if (rank.eq.0) PRINT*, 'A BINARY IS ACCRETED', TTOT,
     &        K-N,NAME(N), name(j1), name(j2)

*     Stellar evolution of two components
            IF (KZ(12).GT.0) THEN
               call sev_one_star_new(J1,kstar(j1),B_RS1,B_L1,
     &              B_TE1,B_MC1,B_RC1,B_M1)
               call sev_one_star_new(J2,kstar(j2),B_RS2,B_L2,
     &              B_TE2,B_MC2,B_RC2,B_M2)
            END IF
*     Binary parameters
            SEMI = -0.5*BODY(K)/H(KSPAIR)
            B_ECC = REAL(SQRT((1.0D0 - R(KSPAIR)/SEMI)**2+
     &           TDOT2(KSPAIR)**2/(BODY(K)*SEMI)))
            B_P = REAL(DAYS*SEMI*SQRT(ABS(SEMI)/BODY(I)))
            B_A = REAL(SEMI*RAU)
*     CONTINUE
         WRITE (28, 96) TTOT,CMBH,NAME(J1),NAME(J2),NAME(K),
     &    KSTAR(J1),KSTAR(J2),KSTAR(K),
     &    BODY(J1),BODY(J2),(X0(J,K), J=1,3),(X0DOT(J,K), J=1,3),
     &    B_A,B_ECC,B_P,B_RS1, B_RS2, B_L1,B_L2,B_TE1,B_TE2
   96 FORMAT(' Bin ',1P,2E16.6,3I10,3I4,17E18.6)
         CALL FLUSH(28)
              IPHASE = 2
*
*            Check for rare case of merged binary component
             IF (NAME(K).LT.0) IPHASE = 7
*      
             GOTO 10
*         
            ENDIF
*       
         X_IJ = X0(1,K)
         Y_IJ = X0(2,K)
         Z_IJ = X0(3,K)
*         
         VX_IJ = X0DOT(1,K)
         VY_IJ = X0DOT(2,K)
         VZ_IJ = X0DOT(3,K)
*
         V = SQRT(VX_IJ**2+VY_IJ**2+VZ_IJ**2)         
*
         R2 = X_IJ**2+Y_IJ**2+Z_IJ**2
         RR = SQRT(R2)
*         
         RV_IJ = VX_IJ*X_IJ+VY_IJ*Y_IJ+VZ_IJ*Z_IJ
         TEMP1 = CMBH/(RR*R2)
*
*	 External force from BH acting on a particle		 
*         
         A_EXT_I(1) = -TeMP1*X_IJ
         A_EXT_I(2) = -TeMP1*Y_IJ
         A_EXT_I(3) = -TeMP1*Z_IJ
*
         ADOT_EXT_I(1) = -TeMP1*(VX_IJ-3.0*RV_IJ*X_IJ/R2)
         ADOT_EXT_I(2) = -TeMP1*(VY_IJ-3.0*RV_IJ*Y_IJ/R2)
         ADOT_EXT_I(3) = -TeMP1*(VZ_IJ-3.0*RV_IJ*Z_IJ/R2)
*            print *, "before sev_one_star:", s_rs, s_l, s_te, kw
                     call sev_one_star_new(K,kw,S_RS,S_L,
     &                    S_TE,S_MC,S_RC,S_M)
*            print *, "after sev_one_star:", s_rs, s_l, s_te, kw
         WRITE (28, 99) TTOT,CMBH,NAME(K),KSTAR(K),BODY(K),
     &      (X0(J,K),J=1,3),(X0DOT(J,K),J=1,3),S_RS,S_L,S_TE
   99 FORMAT(' Sin ',1P,2E16.6,I10,I4,10E18.6)
         CALL FLUSH(28)
                     
* !!!     PUT HERE STARDISC ON/OFF CRITERIUM
*
         E_KIN = 0.5*BODY(K)*(X0DOT(1,K)**2+X0DOT(2,K)**2+X0DOT(3,K)**2)
         E_POT_BH = -CMBH*BODY(K)/SQRT(RSB**2+EPS_bh2)
         NNB = LIST(1,K)
*      EPOT1 = 0.0D0
*      DO 20 J = 1,N
*      IF (J.EQ.K .OR. J.EQ.2*IPAIR-1 .OR. J.EQ.2*IPAIR .OR.
*     *    BODY(J).EQ.0.0D0 .OR. BODY(K).EQ.0.0D0)  GO TO 20
*          A1 = X(1,K) - X(1,J)
*          A2 = X(2,K) - X(2,J)
*          A3 = X(3,K) - X(3,J)
*      A4 = BODY(J)*BODY(K)/DSQRT(A1*A1 + A2*A2 + A3*A3)
*      EPOT1 = EPOT1 - A4
*
         TMP = 0.0D0
*
* !!! The body of do loop below needs to be corrected due to reg and ireg
* !!!      timestep - Done!
* !!!		
        DO 20 ii = IFIRST, NTOT
           IF(ii.NE.K) THEN
              DT = TIME - T0(ii)
              DT1 = 1.5D0*DT
              DT2 = 2.0D0*DT
              X(1:3,ii) = ((FDOT(1:3,ii)*DT + F(1:3,ii))*DT + 
     &                     X0DOT(1:3,ii))*DT + X0(1:3,ii)
              XDOT(1:3,ii) = (FDOT(1:3,ii)*DT1 + F(1:3,ii))*DT2 + 
     &                     X0DOT(1:3,ii)
              RKB2 = X(1,ii)**2+X(2,ii)**2+X(3,ii)**2
              TMP = TMP - BODY(ii)/SQRT(RKB2+EPS_bh2)
              RKS2 = (X0(1,K)-X(1,ii))**2+(X0(2,K)-X(2,ii))**2+
     &               (X0(3,K)-X(3,ii))**2
              TMP = TMP + BODY(ii)/SQRT(RKS2+EPS_ss2)
           END IF
 20     CONTINUE 
*
*       !Accrete the mass of a star:
        CMBH = CMBH + BODY(K)
*       NUM_ACCR = NUM_ACCR+1
*
*		!Create zero-mass particle and put the accreted star very far away from the cluster:
*        E_POT = EPOT1
       E_POT = -BODY(K)*TMP
        E_CORR_TMP = E_KIN + E_POT + E_POT_BH
        E_CORR = E_CORR + E_CORR_TMP
*
       if(rank.eq.0)
     &  PRINT*,'ACCR: RSB,IPH,K,NAME,TIME,EPOT,EKIN,EBH,ECORR,',
     & RSB,IPHASE,K, NAME(K), TTOT, E_POT, E_KIN, E_POT_BH, E_CORR  
* 
       BODY(K) = 0.0
*
*       t0(k) = time
*        T0(K) = TADJ + DTADJ
*        STEP(K) = 1.0D+06
*        STEPR(K) = 1.0D+06
*
        DO 30 L = 1,3
              X0(L,K) = 10000.0*RSCALE*X(L,K)/RR
              X(L,K) = X0(L,K)
              X0DOT(L,K) = SQRT(0.004*ZMASS/RSCALE)*XDOT(L,K)/V
              F(L,K) = 0.0D0
              FDOT(L,K) = 0.0D0
              D0(L,K) = 0.0
              D1(L,K) = 0.0
              D2(L,K) = 0.0
              D3(L,K) = 0.0
              D0R(L,K) = 0.0
              D1R(L,K) = 0.0
              D2R(L,K) = 0.0
              D3R(L,K) = 0.0
 30    CONTINUE
*
*     remove from NXTLST
      call remove_tlist(K,step,dtk)
        step(k)=2.d0*dtk(1)
*     add into GHOST LIST
        call add_tlist(K,STEP,DTK)
        kstar(k) = 15
        tev(k)=1.E6

*     !Remove accreted star from neighbour lists
*	NNB = LIST(1,K)		
        CALL NBREM(K,1,NNB)
        LIST(1,K) = 0
*
*		Set IPHASE < 0 to ensure updating of time-step sequence
*
        IPHASE = -1
        IACCR = 1
*
       if(rank.eq.0) PRINT*, 'BHDATA :',  TTOT, CMBH
*      50 FORMAT ('BHDATA: ',F8.2,E15.6)

*
*      GOTO 10
*
* 40   CONTINUE
*

*  ! end if r_tidal:
       END IF
*      !End DO K=IFIRST,NTOT
10      CONTINUE
*
   
*
      RETURN
*
      END
*----------------------------------------------------------------------------
      subroutine sev_one_star_new(I,KW,RM,LUM,TE,MC,RCC,M1)
*
*
*     Get stellar evolution parameters for one star
*     ---------------------------------------------
*
      include 'common6.h'
      INTEGER I,KW
      REAL*8  LUMS(10),TSCLS(20),GB(10)
      REAL*8  RM8,LUM8,MC8,RCC8,M18,M0
      REAL*4  RM,LUM,TE,MC,RCC,M1
      REAL*8 RSCALE_OUT,MSCALE_OUT,VSCALE_OUT,RAU_OUT,TSCALE_OUT

      IF (KZ(19).EQ.0.AND.KZ(12).EQ.-1) THEN
         RSCALE_OUT=1.0
         MSCALE_OUT=1.0
         VSCALE_OUT=1.0
         RAU_OUT=1.0
         TSCALE_OUT=1.0
      else
         RSCALE_OUT=RBAR
         MSCALE_OUT=ZMBAR
         VSCALE_OUT=VSTAR
         RAU_OUT=RAU
         TSCALE_OUT=TSTAR
      END IF

*
      KW = KSTAR(I)
      AGE = MAX(TIME,TEV0(I))*TSCALE_OUT - EPOCH(I)
      M0 = BODY0(I)*MSCALE_OUT
      M18 = BODY(I)*MSCALE_OUT
      CALL STAR(KW,M0,M18,TM,TN,TSCLS,LUMS,GB,ZPARS)
      CALL HRDIAG(M0,AGE,M18,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &     RM8,LUM8,KW,MC8,RCC8,ME,RE,K2)
*     Temperature (Solar temperature got from Williams, D. R. (1 July 2013). "Sun Fact Sheet". NASA. Retrieved 12 August 2013.)
      TE = REAL(5778*(LUM8/(RM8*RM8))**0.25)

      RM = REAL(RM8)
      LUM = REAL(LUM8)
      MC = REAL(MC8)
      RCC = REAL(RCC8)
      M1 = REAL(M18)

      RETURN

      END

