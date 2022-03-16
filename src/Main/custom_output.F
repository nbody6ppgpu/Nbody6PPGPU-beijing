      subroutine custom_output(NXTLEN,NXTLST,NXTLIMIT,NGHOSTS)
*
*
*     Custom output with HDF5/binary format  (Long Wang)
*
*     Completion and upgrade Aug. 2021 (Rainer Spurzem)
*     Added FDOT, FDSCALE
*-----------------------------------------------------------------------------------
*     IF(KZ(19).EQ.0.AND.KZ(12).EQ.-1) The output is not scaled
*     single: NAME, M[M*], X(1:3)[pc], V(1:3)[km/s], POT[NB],
*             RS[R*], L[L*], Teff[K], MCORE[M*], RSCORE[R*], K*

*     binary: NAME1, NAME2, NAME(ICM), M1[M*], M2[M*], XCM(1:3)[pc], 
*             VCM(1:3)[km/s], XREL(1:3)[AU], VREL(1:3)[km/s], POT[NB],
*             semi[AU], ecc, P[days], Gamma,
*             RS1[R*], RS2[R*], L1[L*], L2[L*], Teff1[K], Teff2[K],
*             MCORE1[M*], MCORE2[M*], RSCORE1[R*], RSCORE2[R*],
*             K*1, K*2, K*(ICM)

*     merger: NAME1,NAME2,NAME3,NAME(INCM), M1[M*], M2[M*], M3[*],
*             XCM(1:3)[pc], VCM(1:3)[km/s], XREL0(1:3)[AU], VREL0(1:3)[km/s],
*             XREL1(1:3)[AU], VREL1(1:3)[km/s],POT[NB],
*             semi0[AU], ecc0, P0[days], semi1[AU], ecc1, P1[days],
*             RS1[R*], RS2[R*], RS3[R*], L1[L*], L2[L*], L3[L*], Teff1[K], Teff2[K], Teff3[K],
*             MCORE1[M*], MCORE2[M*], MCORE3[M*], RSCORE1[R*], RSCORE2[R*], RSCORE3[R*]
*             K*1, K*2, K*3, K*(INCM), 
*     (binary parameter: 0 means inner binary, 1 means outer component )
*     (single parameter: 1/2, means inner binary components, 3 means outer component )
*     (merger parameter: 0 inner binary, 1 outer binary)
*     (POT: potential, sum of - G m_j/r_j)
*-----------------------------------------------------------------------------------
*
      include 'common6.h'
      include 'timing.h'
      include 'galaxy.h'
      include 'output.h'
      include 'scale_out.h'
      PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2,GALIMIT=0.1D0)
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CMC(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
      REAL*8  M,MASS,MC,MIJ,MKK
      REAL*8 TMASK(NMAX),VIREL(3),XIREL(3)
      REAL*8 FREL(3),FDREL(3),D2REL(3),D3REL(3)
      INTEGER NLOOP
*
      INTEGER NXTLST(NMAX),IWBINC
      INTEGER N_SINGLE,N_BINARY,N_MERGER, IMLST(4,30),NIMLST
      REAL*8 MMLST(30),MTMP,M3
      REAL*8 F2DOT(3),F3DOT(3),F2DOT2(3),F3DOT2(3),ONE24,ONE120
      SAVE TMASK
      DATA TMASK/NMAX*-999.0D0/
      PARAMETER  (ONE24=1.0/24.0D0,ONE120=1.0/120.0D0)
*      
      TTOT = TIME+TOFF
      IWBINC = NTOT
      N_SINGLE = 0
      N_BINARY = 0
      N_MERGER = 0
*
*     Try no scaling output if KZ(46)
      IF (KZ(19).EQ.0.AND.KZ(12).EQ.-1) THEN
         RSCALE_OUT=1.0
         MSCALE_OUT=1.0
         DSCALE_OUT=1.0
         VSCALE_OUT=1.0
         RAU_OUT=1.0
         TSCALE_OUT=1.0
         DAYS_OUT=1.0
         FSCALE_OUT=1.0
         FDSCALE_OUT=1.0
         POTSCALE_OUT=1.0
         G_OUT=1.0
         D2SCALE_OUT=1.0
         D3SCALE_OUT=1.0
      else
         RSCALE_OUT=RBAR
         MSCALE_OUT=ZMBAR
         VSCALE_OUT=VSTAR
         DSCALE_OUT=ZMBAR/RBAR**3
         RAU_OUT=RAU
         TSCALE_OUT=TSTAR
         DAYS_OUT=DAYS
         G_OUT=4.302E-3
         POTSCALE_OUT=G_OUT*MSCALE_OUT/RSCALE_OUT
         FSCALE_OUT=G_OUT*(MSCALE_OUT/RSCALE_OUT)**2
         FDSCALE_OUT=FSCALE_OUT/TSCALE_OUT
         D2SCALE_OUT=FDSCALE_OUT/TSCALE_OUT
         D3SCALE_OUT=D2SCALE_OUT/TSCALE_OUT
      END IF

*     choose full or active particles output
      NLOOP = NXTLEN
      IF(KZ(46).EQ.3.OR.KZ(46).EQ.4.OR.TTOT.EQ.0.0D0) THEN
         NLOOP = NXTLIMIT
      END IF
 
*     if(rank.eq.0)print*,' IN ttot,nloop,nxtlen,nxtlimit=',
*    &  ttot,nloop,nxtlen,nxtlimit     
      DO 1 J = 1, NLOOP
         I = NXTLST(J)
*
         S = TIME - T0(I)
         S1 = 1.5*S
         S2 = 2.0*S
*     Note X may become corrected value X0 for zero interval.
         X(1,I) = ((FDOT(1,I)*S + F(1,I))*S + X0DOT(1,I))*S + X0(1,I)
         X(2,I) = ((FDOT(2,I)*S + F(2,I))*S + X0DOT(2,I))*S + X0(2,I)
         X(3,I) = ((FDOT(3,I)*S + F(3,I))*S + X0DOT(3,I))*S + X0(3,I)
         XDOT(1,I) = (FDOT(1,I)*S1 + F(1,I))*S2 + X0DOT(1,I)
         XDOT(2,I) = (FDOT(2,I)*S1 + F(2,I))*S2 + X0DOT(2,I)
         XDOT(3,I) = (FDOT(3,I)*S1 + F(3,I))*S2 + X0DOT(3,I)
*
*      Store F2DOT, F3DOT ready for Taylor series with factor as F,FDOT
         DTR = TIME - T0R(I)
      DO 30 K = 1,3
          F2DOT(K) = ONE24*(D3R(K,I)*DTR+D2R(K,I) + D3(K,I)*S+D2(K,I))
          F3DOT(K) = ONE120*(D3R(K,I) + D3(K,I))
 30   CONTINUE
*     if(rank.eq.0)print*,' i,time,dtr,s,d2,d3=',
*    &  i,time,dtr,s,d2(1,i),d3(1,i)
*     if(rank.eq.0)print*,' i,time,dtr,s,d2r,d3r=',
*    &  i,time,dtr,s,d2r(1,i),d3r(1,i)
*     if(rank.eq.0)print*,' i,time,dtr,s,f2,f3=',
*    &  i,time,dtr,s,f2dot(1),f3dot(1)
*
*     --K.S. binary case--
         IF(I.GT.N.AND.NAME(I).GT.0) THEN
               IPAIR = I - N
*     Resolve unperturbed binary
            CALL KSRES_OP(IPAIR,J1,J2,XIREL,VIREL,
     &        FREL,FDREL,D2REL,D3REL,0)
*     Name
            if (GAMMA(IPAIR).LT.GALIMIT.AND.H(IPAIR).LT.0) THEN
               N_BINARY = N_BINARY + 1
               IK = N_BINARY
               NB_NAMC(IK) = NAME(I)
               NB_NAM1(IK) = NAME(J1)
               NB_NAM2(IK) = NAME(J2)
*     Dynamic information
               B_M1(IK) = REAL(BODY(J1)*MSCALE_OUT)
               B_M2(IK) = REAL(BODY(J2)*MSCALE_OUT)
               B_XC1(IK) = REAL(X(1,I)*RSCALE_OUT) 
               B_XC2(IK) = REAL(X(2,I)*RSCALE_OUT) 
               B_XC3(IK) = REAL(X(3,I)*RSCALE_OUT) 
               B_VC1(IK) = REAL(XDOT(1,I)*VSCALE_OUT) 
               B_VC2(IK) = REAL(XDOT(2,I)*VSCALE_OUT) 
               B_VC3(IK) = REAL(XDOT(3,I)*VSCALE_OUT) 
               B_XR1(IK) = REAL(XIREL(1)*RSCALE_OUT) 
               B_XR2(IK) = REAL(XIREL(2)*RSCALE_OUT) 
               B_XR3(IK) = REAL(XIREL(3)*RSCALE_OUT) 
               B_VR1(IK) = REAL(VIREL(1)*VSCALE_OUT) 
               B_VR2(IK) = REAL(VIREL(2)*VSCALE_OUT) 
               B_VR3(IK) = REAL(VIREL(3)*VSCALE_OUT) 
               B_FC1(IK) = REAL(F(1,I)*FSCALE_OUT)
               B_FC2(IK) = REAL(F(2,I)*FSCALE_OUT)
               B_FC3(IK) = REAL(F(3,I)*FSCALE_OUT)
               B_FDC1(IK) = REAL(FDOT(1,I)*FDSCALE_OUT)
               B_FDC2(IK) = REAL(FDOT(2,I)*FDSCALE_OUT)
               B_FDC3(IK) = REAL(FDOT(3,I)*FDSCALE_OUT)
               B_D21(IK) = REAL(F2DOT(1)*D2SCALE_OUT)
               B_D22(IK) = REAL(F2DOT(2)*D2SCALE_OUT)
               B_D23(IK) = REAL(F2DOT(3)*D2SCALE_OUT)
               B_D31(IK) = REAL(F3DOT(1)*D3SCALE_OUT)
               B_D32(IK) = REAL(F3DOT(2)*D3SCALE_OUT)
               B_D33(IK) = REAL(F3DOT(3)*D3SCALE_OUT)
               B_XR1(IK) = REAL(XIREL(1)*RSCALE_OUT)
               B_XR2(IK) = REAL(XIREL(2)*RSCALE_OUT)
               B_XR3(IK) = REAL(XIREL(3)*RSCALE_OUT)
               B_VR1(IK) = REAL(VIREL(1)*VSCALE_OUT)
               B_VR2(IK) = REAL(VIREL(2)*VSCALE_OUT)
               B_VR3(IK) = REAL(VIREL(3)*VSCALE_OUT)
               B_FR1(IK) = REAL(FREL(1)*FSCALE_OUT)
               B_FR2(IK) = REAL(FREL(2)*FSCALE_OUT)
               B_FR3(IK) = REAL(FREL(3)*FSCALE_OUT)
               B_FDR1(IK) = REAL(FDREL(1)*FDSCALE_OUT)
               B_FDR2(IK) = REAL(FDREL(2)*FDSCALE_OUT)
               B_FDR3(IK) = REAL(FDREL(3)*FDSCALE_OUT)
               B_D21R(IK) = REAL(D2REL(1)*D2SCALE_OUT)
               B_D22R(IK) = REAL(D2REL(2)*D2SCALE_OUT)
               B_D23R(IK) = REAL(D2REL(3)*D2SCALE_OUT)
               B_D31R(IK) = REAL(D3REL(1)*D3SCALE_OUT)
               B_D32R(IK) = REAL(D3REL(2)*D3SCALE_OUT)
               B_D33R(IK) = REAL(D3REL(3)*D3SCALE_OUT)
               B_POT(IK) = REAL(-PHIDBL(I)*POTSCALE_OUT)
               B_G(IK)   = REAL(GAMMA(IPAIR))
               B_STEP(IK) = REAL(STEP(I)*TSCALE_OUT)
               B_STEPR(IK) = REAL(STEPR(I)*TSCALE_OUT)
               B_T0(IK) = REAL(T0(I)*TSCALE_OUT)
               B_T0R(IK) = REAL(T0R(I)*TSCALE_OUT)
*     Stellar evolution of two components
               IF (KZ(12).GT.0) THEN
                  call sev_one_star(J1,NB_KW1(IK),B_RS1(IK),B_L1(IK),
     &                 B_TE1(IK),B_MC1(IK),B_RC1(IK),B_M1(IK))
                  call sev_one_star(J2,NB_KW2(IK),B_RS2(IK),B_L2(IK),
     &                 B_TE2(IK),B_MC2(IK),B_RC2(IK),B_M2(IK))
               END IF
               NB_KWC(IK) = KSTAR(I)
*     Binary parameters
               SEMI = -0.5*BODY(I)/H(IPAIR)
               B_ECC(IK) = REAL(SQRT((1.0D0 - R(IPAIR)/SEMI)**2+
     &              TDOT2(IPAIR)**2/(BODY(I)*SEMI)))
               B_P(IK) = REAL(DAYS_OUT*SEMI*SQRT(ABS(SEMI)/BODY(I)))
               B_A(IK) = REAL(SEMI*RAU_OUT)
            END IF
*     If perturbation is too strong (Gamma>GALIMIT) or the binding energy is positive, the K.S. may not be a binary.
*     Save binary members and c.m. as single particles in any case (R.Sp. Aug 2021).
               TMASK(J1) = TTOT
               TMASK(J2) = TTOT
               TMASK(I)  = TTOT
               N_SINGLE = N_SINGLE + 1
               IK = N_SINGLE
               call custom_store_single(J1,IK)
*     as of Oct. 2021 no correction here, use relative quantities (R.Sp.)
*     correction for the force and potential / D2 / D3 missing R.Sp.
*              DR3I = BODY(J2)/R(IPAIR)**3
*              S_F1(IK) = S_F1(IK) - REAL(XIREL(1)*DR3I*FSCALE_OUT)
*              S_F2(IK) = S_F2(IK) - REAL(XIREL(2)*DR3I*FSCALE_OUT)
*              S_F3(IK) = S_F3(IK) - REAL(XIREL(3)*DR3I*FSCALE_OUT)
*              S_POT(IK) = S_POT(IK) -
*    &              REAL(BODY(J2)/R(IPAIR)*POTSCALE_OUT)
*
               N_SINGLE = N_SINGLE + 1
               IK = N_SINGLE
               call custom_store_single(J2,IK)
*     correction for the force and potential
*              DR3I = -BODY(J1)/R(IPAIR)**3
*              S_F1(IK) = S_F1(IK) - REAL(XIREL(1)*DR3I*FSCALE_OUT) 
*              S_F2(IK) = S_F2(IK) - REAL(XIREL(2)*DR3I*FSCALE_OUT) 
*              S_F3(IK) = S_F3(IK) - REAL(XIREL(3)*DR3I*FSCALE_OUT) 
*              S_POT(IK) = S_POT(IK) -
*    &              REAL(BODY(J1)/R(IPAIR)*POTSCALE_OUT)
*     add also c.m. (R.Sp. Aug 2021).
*              N_SINGLE = N_SINGLE + 1
*              IK = N_SINGLE
*              call custom_store_single(I,IK)
*             
*     CONTINUE
            GO TO 1
         END IF
*     
*     --Merger case--
         IF(I.GT.N.AND.NAME(I).LT.0) THEN
            NIMLST = 0
            NAMEI = NAME(I)
            ICM = I
            IPAIR = ICM - N
            I2 = 2*IPAIR
            I1 = I2 - 1
*     Check possible quartet [[B,S],S] and quintet [[B,S],B], save two components seperately
 4          IF (NAMEI.LT.-2*NZERO) THEN
*     Merger index
               IM = 0
*     Ghost index
               IG = 0
*     Notice, here assume merger number is small
               DO K = 1, NMERGE
                  IF (NAMEM(K).EQ.NAMEI) THEN
                     IM = K
                     GO TO 10
                  END IF
               END DO
*     error messsage
               print*,'Error!: Cannot find merger index for name =',
     &              NAMEI, ', NMERGER =',NMERGE
               call flush(6)
               call abort()
*
 10            DO K = NXTLIMIT+1,NXTLIMIT+NGHOSTS
                  IF (NAME(NXTLST(K)).EQ.NAMEG(IM)) THEN
                     IG = NXTLST(K)
                     GO TO 20
                  END IF
               END DO
*     error messsage
               print*,'Error!: Cannot find ghost of merger IM=',
     &              IM, ', NAMEG(IM) =',NAMEG(IM)
               call flush(6)
               call abort()
*
 20            CONTINUE
               CALL KSRES_OP(IPAIR,J1,J2,XIREL,VIREL,
     &        FREL,FDREL,D2REL,D3REL,1)
*     Outer component
               IF (NAME(IG).LT.0) THEN
                  NIMLST = NIMLST + 1
                  IMLST(1,NIMLST) = J2
                  IMLST(2,NIMLST) = IG
                  IMLST(3,NIMLST) = NAME(IG)
                  MMLST(NIMLST) = CM(4,IM)
               ELSE IF(NAME(IG).GT.NZERO) THEN
                  N_BINARY = N_BINARY + 1
                  JPAIR = IG - N
           CALL KSRES_OP(IPAIR,J1,J2,XIREL,VIREL,
     &        FREL,FDREL,D2REL,D3REL,0)
                  call custom_store_binary(N_BINARY,J2,NAME(IG),
     &                JP1,JP2,CM(3,IM),CM(4,IM),H(JPAIR),R(JPAIR),
     &                TDOT2(JPAIR),XIREL,VIREL,FREL,FDREL,D2REL,D3REL,
     &                KSTAR(IG),GAMMA(JPAIR))
               ELSE
                  N_SINGLE = N_SINGLE + 1
                  IK = N_SINGLE
*     Name
                  NS_NAM(IK) = NAME(IG)
*     Dynamics
                  S_M(IK)  = REAL(BODY(J2)*MSCALE_OUT)
                  S_RSN(IK) = REAL(RS(J2)*RSCALE_OUT)
                  S_X1(IK) = REAL(X(1,J2)*RSCALE_OUT) 
                  S_X2(IK) = REAL(X(2,J2)*RSCALE_OUT) 
                  S_X3(IK) = REAL(X(3,J2)*RSCALE_OUT) 
                  S_V1(IK) = REAL(XDOT(1,J2)*VSCALE_OUT) 
                  S_V2(IK) = REAL(XDOT(2,J2)*VSCALE_OUT) 
                  S_V3(IK) = REAL(XDOT(3,J2)*VSCALE_OUT) 
                  S_F1(IK) = REAL(F(1,I)*FSCALE_OUT) 
                  S_F2(IK) = REAL(F(2,I)*FSCALE_OUT) 
                  S_F3(IK) = REAL(F(3,I)*FSCALE_OUT) 
                  S_FD1(IK) = REAL(FDOT(1,I)*FDSCALE_OUT) 
                  S_FD2(IK) = REAL(FDOT(2,I)*FDSCALE_OUT) 
                  S_FD3(IK) = REAL(FDOT(3,I)*FDSCALE_OUT) 
                  S_D21(IK) = REAL(F2DOT(1)*D2SCALE_OUT)
                  S_D22(IK) = REAL(F2DOT(2)*D2SCALE_OUT)
                  S_D23(IK) = REAL(F2DOT(3)*D2SCALE_OUT)
                  S_D31(IK) = REAL(F3DOT(1)*D3SCALE_OUT)
                  S_D32(IK) = REAL(F3DOT(2)*D3SCALE_OUT)
                  S_D33(IK) = REAL(F3DOT(3)*D3SCALE_OUT)
                  S_POT(IK) = REAL(-PHIDBL(I)*POTSCALE_OUT)
                  S_STEP(IK) = REAL(STEP(I)*TSCALE_OUT)
                  S_STEPR(IK) = REAL(STEPR(I)*TSCALE_OUT)
                  S_T0(IK) = REAL(T0(I)*TSCALE_OUT)
                  S_T0R(IK) = REAL(T0R(I)*TSCALE_OUT)

*     Stellar evolution of single star
                  IF (KZ(12).GT.0) THEN
                     MTMP = BODY(IG)
                     BODY(IG) = BODY(J2)
                     call sev_one_star(IG,NS_KW(IK),S_RS(IK),S_L(IK),
     &                    S_TE(IK),S_MC(IK),S_RC(IK),S_M(IK))
                     BODY(IG) = MTMP
                  END IF
               END IF
*     Inner merger
               NIMLST = NIMLST + 1
               IMLST(1,NIMLST) = J1
               IMLST(2,NIMLST) = -IM
               IMLST(3,NIMLST) = NAMEI + 2*NZERO
               MMLST(NIMLST) = CM(2,IM)
               GO TO 8
            END IF
*     For normal merger
            N_MERGER = N_MERGER + 1
            IK = N_MERGER
            DO K = 1, NMERGE
               IF (NAMEM(K).EQ.NAMEI) THEN
                  IM = K
                  GO TO 6
               END IF
            END DO
*     error messsage
            print*,'Error!: Cannot find merger index for name =',
     &           NAMEI, ', NMERGER =',NMERGE
            call flush(6)
            call abort()
*     Check B-B system, save as two binaries
 6          IF (NAME(I2).GT.NZERO) THEN
               DO K = NXTLIMIT+1,NXTLIMIT+NGHOSTS
                  IF (NAME(NXTLST(K)).EQ.NAMEG(IM)) THEN
                     IG = NXTLST(K)
                     GO TO 17
                  END IF
               END DO
*     error messsage
               print*,'Error!: Cannot find ghost of merger IM=',
     &              IM, ', NAMEG(IM) =',NAMEG(IM)
               call flush(6)
               call abort()
 17            CONTINUE
           CALL KSRES_OP(IPAIR,J1,J2,XIREL,VIREL,
     &        FREL,FDREL,D2REL,D3REL,1)
*     Inner binary
               N_BINARY = N_BINARY + 1
               RIP = SQRT(XREL(1,IM)**2+XREL(2,IM)**2+XREL(3,IM)**2)
               RDOT = XREL(1,IM)*VREL(1,IM) + XREL(2,IM)*VREL(2,IM) + 
     &              XREL(3,IM)*VREL(3,IM)
               call custom_store_binary(N_BINARY,J1,NAMEM(IM),J1,IG,
     &              CM(1,IM),CM(2,IM),HM(IM),RIP,RDOT,
     &              XREL(1,IM),VREL(1,IM),FREL,FDREL,D2REL,D3REL,
     &              KSTARM(IM),GALIMIT)
*     Outer binary
               N_BINARY = N_BINARY + 1
               JPAIR = IG - N
           CALL KSRES_OP(JPAIR,JP1,JP2,XIREL,VIREL,
     &        FREL,FDREL,D2REL,D3REL,0)
               call custom_store_binary(N_BINARY,J2,NAME(I2),JP1,JP2,
     &             CM(3,IM),CM(4,IM),H(JPAIR),R(JPAIR),TDOT2(JPAIR),
     &             XIREL,VIREL,FREL,FDREL,D2REL,D3REL,KSTAR(I2),GALIMIT)
               GO TO 8
            END IF
*     Resolve merger
            M_XC1(IK) = REAL(X(1,ICM)*RSCALE_OUT) 
            M_XC2(IK) = REAL(X(2,ICM)*RSCALE_OUT) 
            M_XC3(IK) = REAL(X(3,ICM)*RSCALE_OUT) 
            M_VC1(IK) = REAL(XDOT(1,ICM)*VSCALE_OUT) 
            M_VC2(IK) = REAL(XDOT(2,ICM)*VSCALE_OUT) 
            M_VC3(IK) = REAL(XDOT(3,ICM)*VSCALE_OUT)
            M_FC1(IK) = REAL(F(1,ICM)*FSCALE_OUT)
            M_FC2(IK) = REAL(F(2,ICM)*FSCALE_OUT)
            M_FC3(IK) = REAL(F(3,ICM)*FSCALE_OUT)
            M_POT(IK) = REAL(-PHIDBL(ICM)*POTSCALE_OUT)

*     Check whether ipair is avaiable
            IF(IPAIR.GT.0) THEN
           CALL KSRES_OP(IPAIR,J1,J2,XIREL,VIREL,
     &        FREL,FDREL,D2REL,D3REL,0)
               M_XR11(IK) = REAL(XIREL(1)*RSCALE_OUT)
               M_XR12(IK) = REAL(XIREL(2)*RSCALE_OUT)
               M_XR13(IK) = REAL(XIREL(3)*RSCALE_OUT)
               M_VR11(IK) = REAL(VIREL(1)*VSCALE_OUT)
               M_VR12(IK) = REAL(VIREL(2)*VSCALE_OUT)
               M_VR13(IK) = REAL(VIREL(3)*VSCALE_OUT)
            ELSE
               IM2 = -IPAIR
               M_XR11(IK) = REAL(XREL(1,IM2)*RSCALE_OUT) 
               M_XR12(IK) = REAL(XREL(2,IM2)*RSCALE_OUT) 
               M_XR13(IK) = REAL(XREL(3,IM2)*RSCALE_OUT) 
               M_VR11(IK) = REAL(VREL(1,IM2)*VSCALE_OUT) 
               M_VR12(IK) = REAL(VREL(2,IM2)*VSCALE_OUT) 
               M_VR13(IK) = REAL(VREL(3,IM2)*VSCALE_OUT) 
            END IF
            M_XR01(IK) = REAL(XREL(1,IM)*RSCALE_OUT) 
            M_XR02(IK) = REAL(XREL(2,IM)*RSCALE_OUT) 
            M_XR03(IK) = REAL(XREL(3,IM)*RSCALE_OUT) 
            M_VR01(IK) = REAL(VREL(1,IM)*VSCALE_OUT) 
            M_VR02(IK) = REAL(VREL(2,IM)*VSCALE_OUT) 
            M_VR03(IK) = REAL(VREL(3,IM)*VSCALE_OUT) 
            M_M1(IK) = REAL(CM(1,IM)*MSCALE_OUT) 
            M_M2(IK) = REAL(CM(2,IM)*MSCALE_OUT) 
            M_M3(IK) = REAL(BODY(J2)*MSCALE_OUT)
            
            NM_NAMC(IK) = NAMEI
            NM_NAM1(IK) = NAME(J1)
            NM_NAM2(IK) = NAMEG(IM)
            NM_NAM3(IK) = NAME(J2)

            SEMI = 0.D0
            IF (IPAIR.GT.0) THEN
               SEMI = -0.5*BODY(ICM)/H(IPAIR)
               M_ECC1(IK) = REAL(SQRT((1.0D0 - R(IPAIR)/SEMI)**2 + 
     &              TDOT2(IPAIR)**2/(BODY(ICM)*SEMI)))
            ELSE
               IM2 = -IPAIR
               SEMI = -0.5*BODY(ICM)/HM(IM2)
               RJ = SQRT(XREL(1,IM2)**2+XREL(2,IM2)**2+XREL(3,IM2)**2)
C               RDOT = XREL(1,IM2)*VREL(1,IM2) + XREL(2,IM2)*VREL(2,IM2)
C     &              + XREL(3,IM2)*VREL(3,IM2)
               RDOT = 2.0*(UM(1,IM2)*UMDOT(1,IM2)+UM(2,IM2)*UMDOT(2,IM2)
     &              +UM(3,IM2)*UMDOT(3,IM2)+UM(4,IM2)*UMDOT(4,IM2))
               M_ECC1(IK) = REAL(SQRT((1.0D0 - RJ/SEMI)**2 + 
     &              RDOT**2/(BODY(ICM)*SEMI)))
            END IF
            M_P1(IK) = REAL(DAYS_OUT*SEMI*SQRT(ABS(SEMI)/BODY(ICM)))
            M_A1(IK) = REAL(SEMI*RAU_OUT)
            SEMI0 = -0.5*BODY(J1)/HM(IM)
            RJ = SQRT(XREL(1,IM)**2 + XREL(2,IM)**2 + XREL(3,IM)**2)
C            RDOT = XREL(1,IM)*VREL(1,IM) + XREL(2,IM)*VREL(2,IM)
C     &           + XREL(3,IM)*VREL(3,IM)
            RDOT = 2.0*(UM(1,IM)*UMDOT(1,IM) + UM(2,IM)*UMDOT(2,IM) +
     &           UM(3,IM)*UMDOT(3,IM) + UM(4,IM)*UMDOT(4,IM))
            M_ECC0(IK) = REAL(SQRT((1.0D0 - RJ/SEMI0)**2 +
     &           RDOT**2/(BODY(J1)*SEMI0)))
            M_P0(IK) = REAL(DAYS_OUT*SEMI0*SQRT(ABS(SEMI0)/BODY(J1)))
            M_A0(IK) = REAL(SEMI0*RAU_OUT)
            
*     stellar evolution
            IF (KZ(12).GT.0) THEN
               DO K = NXTLIMIT+1,NXTLIMIT+NGHOSTS
                  IF (NAME(NXTLST(K)).EQ.NAMEG(IM)) THEN
                     IG = NXTLST(K)
                     GO TO 7
                  END IF
               END DO
*     error messsage
               print*,'Error!: Cannot find ghost of merger IM=',
     &              IM, ', NAMEG(IM) =',NAMEG(IM)
               call flush(6)
               call abort()
*     first star
 7             MTMP = BODY(J1)
               BODY(J1) = CM(1,IM)
               call sev_one_star(J1,NM_KW1(IK),M_RS1(IK),M_L1(IK),
     &              M_TE1(IK),M_MC1(IK),M_RC1(IK),M_M1(IK))
               BODY(J1) = MTMP
*     second star
               MTMP = BODY(IG)
               BODY(IG) = CM(2,IM)
               call sev_one_star(IG,NM_KW2(IK),M_RS2(IK),M_L2(IK),
     &              M_TE2(IK),M_MC2(IK),M_RC2(IK),M_M2(IK))
               BODY(IG) = MTMP
*     third star
               call sev_one_star(J2,NM_KW3(IK),M_RS3(IK),M_L3(IK),
     &              M_TE3(IK),M_MC3(IK),M_RC3(IK),M_M3(IK))
               NM_KWC(IK) = KSTAR(ICM)
            END IF
 8          IF(NIMLST.GT.0) THEN
*     ICM (M,X,XD), IPAIR(XIREL,VIREL), NAMEI, M3
               ICM = IMLST(1,NIMLST)
               IPAIR = IMLST(2,NIMLST)
               NAMEI = IMLST(3,NIMLST)
               M3 = MMLST(NIMLST)
               NIMLST = NIMLST - 1
               GO TO 4
            END IF
            GO TO 1
         END IF
*
*     --Chain case--
         IF(NAME(I).EQ.0) THEN
            DO L = 1,NCH
               N_SINGLE = N_SINGLE + 1
               IK = N_SINGLE
               NS_NAM(IK) = NAMEC(L)
               S_M(IK)  = REAL(BODYC(L)*MSCALE_OUT)
               S_RSN(IK) = REAL(RS(L)*RSCALE_OUT)
               S_X1(IK) = REAL((XCH(3*L-2)+X(1,ICH))*RSCALE_OUT)
               S_X2(IK) = REAL((XCH(3*L-1)+X(2,ICH))*RSCALE_OUT)
               S_X3(IK) = REAL((XCH(3*L)+X(3,ICH))*RSCALE_OUT)
               S_V1(IK) = REAL((VCH(3*L-2)+XDOT(1,ICH))*VSCALE_OUT)
               S_V2(IK) = REAL((VCH(3*L-1)+XDOT(2,ICH))*VSCALE_OUT)
               S_V3(IK) = REAL((VCH(3*L)+XDOT(3,ICH))*VSCALE_OUT)
               S_F1(IK) = REAL(F(1,I)*FSCALE_OUT) 
               S_F2(IK) = REAL(F(2,I)*FSCALE_OUT) 
               S_F3(IK) = REAL(F(3,I)*FSCALE_OUT) 
               S_POT(IK) = REAL(-PHIDBL(I)*POTSCALE_OUT)
*     stellar evolution         
               IF(KZ(12).GT.0) THEN
                  IF (L.EQ.1) THEN
                     MTMP = BODY(I)
                     BODY(I) = BODYC(L)
                     call sev_one_star(I,NS_KW(IK),S_RS(IK),S_L(IK),
     &                    S_TE(IK),S_MC(IK),S_RC(IK),S_M(IK))
                     BODY(I) = MTMP
                  ELSE
                     IG = 0
                     DO K = NXTLIMIT+1,NXTLIMIT+NGHOSTS
                        IF(NAME(NXTLST(K)).EQ.NAMEC(L)) THEN
                           IG = NXTLST(K)
                           GO TO 11
                        END IF
                     END DO
*     error message
                     print*,'Error!: Chain ghost name =',NAMEC(L),
     &                    'not found.'
                     call flush(6)
                     call abort()
 11                  MTMP = BODY(IG)
                     BODY(IG) = BODYC(L)
                     call sev_one_star(IG,NS_KW(IK),S_RS(IK),S_L(IK),
     &                    S_TE(IK),S_MC(IK),S_RC(IK),S_M(IK))
                     BODY(IG) = MTMP
                  END IF
               END IF
            END DO
*     continue
            GO TO 1
         END IF
*
*     Nearest neighbor index
         I2 = IMINR(I)
         IF(I2.LE.0) I2 = I + 1
 12      IF(BODY(I2).EQ.0.D0.AND.I2.LE.N) THEN
            I2 = I2 + 1
            GO TO 12
         END IF
            
*     --Single star case--
*     escape the registered wide binary index or second component is multiple system
         IF(TMASK(I).EQ.TTOT.OR.TMASK(I2).EQ.TTOT
     &        .OR.I2.GT.N.OR.NAME(I2).LE.0) THEN
            IF(TMASK(I).NE.TTOT) THEN
               TMASK(I) = TTOT
               N_SINGLE = N_SINGLE + 1
*           if(rank.eq.0.and.abs(ttot-1.d0).lt.1.e-04)
*    &    write(6,*)' s ',i,name(i),tmask(i),n_single
               call custom_store_single(I,N_SINGLE)
            END IF
*     CONTINUE
            GO TO 1
*
         ELSE
*     --Wide binary case--
            ttime = TIME
            call JPRED_int(I2,ttime)
            RIJ2 = 0.D0
            VIREL2 = 0.D0
            RDOT = 0.D0
            S2 = TIME - T0(I2)
            DTR2 = TIME - T0R(I2)
            DO K = 1,3
               XIREL(K) = X(K,I) - X(K,I2)
               VIREL(K) = XDOT(K,I) - XDOT(K,I2)
               RIJ2 = RIJ2 + XIREL(K)**2
               VIREL2 = VIREL2 + VIREL(K)**2
               RDOT = RDOT + XIREL(K)*VIREL(K)
*      Store F2DOT2, F3DOT2 for I2 ready for Taylor series with factor as F,FDOT
      F2DOT2(K) = ONE24*(D3R(K,I2)*DTR2+D2R(K,I2)+D3(K,I2)*S2+D2(K,I2))
      F3DOT2(K) = ONE120*(D3R(K,I2) + D3(K,I2))
            END DO
            RIJ = SQRT(RIJ2)
            ZMB = BODY(I) + BODY(I2)
            SEMI = 2.0/RIJ - VIREL2/ZMB
            SEMI = 1.0/SEMI
*     Check whether satisfy wide binary criterion
            IF (SEMI.GT.0.0.AND.SEMI.LT.0.1*RSCALE) THEN
               IWBINC = IWBINC + 1
               N_BINARY = N_BINARY + 1
               IK = N_BINARY
*     set time flag to avoid duplicate case
               TMASK(I) = TTOT
               TMASK(I2) = TTOT
*           if(rank.eq.0.and.abs(ttot-1.d0).lt.1.e-04)
*    &    write(6,*)' wb ',i,name(i),tmask(i),n_binary,i2,name(i2),
*    &                tmask(i2)
*     Name
               NB_NAMC(IK) = IWBINC
               NB_NAM1(IK) = NAME(I)
               NB_NAM2(IK) = NAME(I2)
*     Dynamic information
               B_M1(IK) = REAL(BODY(I)*MSCALE_OUT)
               B_M2(IK) = REAL(BODY(I2)*MSCALE_OUT)

               B_XC1(IK) = REAL((BODY(I)*X(1,I)+BODY(I2)*X(1,I2))
     &              /ZMB*RSCALE_OUT)
               B_XC2(IK) = REAL((BODY(I)*X(2,I)+BODY(I2)*X(2,I2))
     &              /ZMB*RSCALE_OUT)
               B_XC3(IK) = REAL((BODY(I)*X(3,I)+BODY(I2)*X(3,I2))
     &              /ZMB*RSCALE_OUT)
               B_VC1(IK) = REAL((BODY(I)*XDOT(1,I)+BODY(I2)*XDOT(1,I2))
     &              /ZMB*VSCALE_OUT)
               B_VC2(IK) = REAL((BODY(I)*XDOT(2,I)+BODY(I2)*XDOT(3,I2))
     &              /ZMB*VSCALE_OUT)
               B_VC3(IK) = REAL((BODY(I)*XDOT(3,I)+BODY(I2)*XDOT(3,I2))
     &              /ZMB*VSCALE_OUT)
               B_FC1(IK) = REAL((BODY(I)*F(1,I)+BODY(I2)*F(1,I2))
     &              /ZMB*FSCALE_OUT)
               B_FC2(IK) = REAL((BODY(I)*F(2,I)+BODY(I2)*F(2,I2))
     &              /ZMB*FSCALE_OUT)
               B_FC3(IK) = REAL((BODY(I)*FDOT(3,I)+BODY(I2)*FDOT(3,I2))
     &              /ZMB*FSCALE_OUT)
               B_FDC1(IK) = REAL((BODY(I)*FDOT(1,I)+BODY(I2)*FDOT(1,I2))
     &              /ZMB*FDSCALE_OUT)
               B_FDC2(IK) = REAL((BODY(I)*FDOT(2,I)+BODY(I2)*FDOT(2,I2))
     &              /ZMB*FDSCALE_OUT)
               B_FDC3(IK) = REAL((BODY(I)*FDOT(3,I)+BODY(I2)*FDOT(3,I2))
     &              /ZMB*FDSCALE_OUT)
               B_D21(IK) = REAL((BODY(I)*F2DOT(1)+BODY(I2)*F2DOT2(1))
     &              /ZMB*D2SCALE_OUT)
               B_D22(IK) = REAL((BODY(I)*F2DOT(2)+BODY(I2)*F2DOT2(2))
     &              /ZMB*D2SCALE_OUT)
               B_D23(IK) = REAL((BODY(I)*F2DOT(3)+BODY(I2)*F2DOT2(3))
     &              /ZMB*D2SCALE_OUT)
               B_D31(IK) = REAL((BODY(I)*F3DOT(1)+BODY(I2)*F3DOT2(1))
     &              /ZMB*D3SCALE_OUT)
               B_D32(IK) = REAL((BODY(I)*F3DOT(2)+BODY(I2)*F3DOT2(2))
     &              /ZMB*D3SCALE_OUT)
               B_D33(IK) = REAL((BODY(I)*F3DOT(3)+BODY(I2)*F3DOT2(3))
     &              /ZMB*D3SCALE_OUT)
               B_XR1(IK) = REAL(XIREL(1)*RSCALE_OUT) 
               B_XR2(IK) = REAL(XIREL(2)*RSCALE_OUT) 
               B_XR3(IK) = REAL(XIREL(3)*RSCALE_OUT) 
               B_VR1(IK) = REAL(VIREL(1)*VSCALE_OUT) 
               B_VR2(IK) = REAL(VIREL(2)*VSCALE_OUT) 
               B_VR3(IK) = REAL(VIREL(3)*VSCALE_OUT) 
               B_FR1(IK) = REAL(FREL(1)*FSCALE_OUT)
               B_FR2(IK) = REAL(FREL(2)*FSCALE_OUT)
               B_FR3(IK) = REAL(FREL(3)*FSCALE_OUT)
               B_FDR1(IK) = REAL(FDREL(1)*FDSCALE_OUT)
               B_FDR2(IK) = REAL(FDREL(2)*FDSCALE_OUT)
               B_FDR3(IK) = REAL(FDREL(3)*FDSCALE_OUT)
               B_D21R(IK) = REAL(D2REL(1)*D2SCALE_OUT)
               B_D22R(IK) = REAL(D2REL(2)*D2SCALE_OUT)
               B_D23R(IK) = REAL(D2REL(3)*D2SCALE_OUT)
               B_D31R(IK) = REAL(D3REL(1)*D3SCALE_OUT)
               B_D32R(IK) = REAL(D3REL(2)*D3SCALE_OUT)
               B_D33R(IK) = REAL(D3REL(3)*D3SCALE_OUT)
               B_POT(IK) = REAL(-PHIDBL(I)*POTSCALE_OUT)
*     Assume small perturbation (R.Sp. Oct. 2021)
               B_G(IK) = 0.D0
*     Stellar evolution of two components
               IF (KZ(12).GT.0) THEN
                  call sev_one_star(I,NB_KW1(IK),B_RS1(IK),B_L1(IK),
     &                 B_TE1(IK),B_MC1(IK),B_RC1(IK),B_M1(IK))
                  call sev_one_star(I2,NB_KW2(IK),B_RS2(IK),B_L2(IK),
     &                 B_TE2(IK),B_MC2(IK),B_RC2(IK),B_M2(IK))
               END IF
               NB_KWC(IK) = -1
*     Binary parameters
               B_ECC(IK) = REAL(SQRT((1.0D0 - RIJ/SEMI)**2 + 
     &              RDOT**2/(ZMB*SEMI)))
               BODYI = (B_M1(IK) + B_M2(IK))/MSCALE_OUT
               B_P(IK) = REAL(DAYS_OUT*SEMI*SQRT(ABS(SEMI)/BODYI))
               B_A(IK) = REAL(SEMI*RAU_OUT)
*     Treat also as single star (R.Sp. added Aug. 2021).
               TMASK(I) = TTOT
               TMASK(I2) = TTOT
*          if(rank.eq.0.and.abs(ttot-1.d0).lt.1.e-04)
*    &    write(6,*)' ws ',i,name(i),tmask(i),n_single
               N_SINGLE = N_SINGLE + 1
               call custom_store_single(I,N_SINGLE)
*          if(rank.eq.0.and.abs(ttot-1.d0).lt.1.e-04)
*    &    write(6,*)' ws ',i,name(i),tmask(i),n_single
               N_SINGLE = N_SINGLE + 1
               call custom_store_single(I2,N_SINGLE)
               GO TO 1
            ELSE
*     Treat as single star
               TMASK(I) = TTOT
*           if(rank.eq.0.and.abs(ttot-1.d0).lt.1.e-04)
*    &    write(6,*)' s2 ',i,name(i),tmask(i),n_single
               N_SINGLE = N_SINGLE + 1
               call custom_store_single(I,N_SINGLE)
*     continue
               GO TO 1
            END IF
         END IF
*
 1    CONTINUE

*     Prepare global quantities for writing
          S_AS(1) = REAL(TTOT)
          S_AS(2) = REAL(FLOAT(NPAIRS))
          S_AS(3) = REAL(RBAR)
          S_AS(4) = REAL(ZMBAR)
          S_AS(5) = REAL(FLOAT(N))
          S_AS(6) = REAL(TSTAR)
          S_AS(7) = REAL(RDENS(1))
          S_AS(8) = REAL(RDENS(2))
          S_AS(9) = REAL(RDENS(3))
          S_AS(10) = REAL(TTOT/TCR0)
          S_AS(11) = REAL(TSCALE)
          S_AS(12) = REAL(VSTAR)
          S_AS(13) = REAL(RC)
          S_AS(14) = REAL(NC)
          S_AS(15) = REAL(VC)
          S_AS(16) = REAL(RHOM)
          S_AS(17) = REAL(CMAX)
          S_AS(18) = REAL(RSCALE)
          S_AS(19) = REAL(RSMIN)
          S_AS(20) = REAL(DMIN1)
          S_AS(21) = REAL(RG(1))
          S_AS(22) = REAL(RG(2))
          S_AS(23) = REAL(RG(3))
          S_AS(24) = REAL(VG(1))
          S_AS(25) = REAL(VG(2))
          S_AS(26) = REAL(VG(3))
          S_AS(27) = REAL(TIDAL(1))
          S_AS(28) = REAL(TIDAL(2))
          S_AS(29) = REAL(TIDAL(3))
          S_AS(30) = REAL(TIDAL(4))
          S_AS(31) = REAL(GMG)
          S_AS(32) = REAL(OMEGA)
          S_AS(33) = REAL(DISK)
          S_AS(34) = REAL(A)
          S_AS(35) = REAL(B)
          S_AS(36) = REAL(ZMET)
          S_AS(37:56) = REAL(ZPARS(1:20))
          S_AS(57) = REAL(ETAI)
          S_AS(58) = REAL(ETAR)
          S_AS(59) = REAL(ETAU)
          S_AS(60) = REAL(ECLOSE)
          S_AS(61) = REAL(DTMIN)
          S_AS(62) = REAL(RMIN)
          S_AS(63) = REAL(GMIN)
          S_AS(64) = REAL(GMAX)
          S_AS(65) = REAL(SMAX)
          S_AS(66) = REAL(FLOAT(NNBOPT))
          S_AS(67) = REAL(EPOCH0)
          S_AS(68) = REAL(FLOAT(N_SINGLE))
          S_AS(69) = REAL(FLOAT(N_BINARY))
          S_AS(70) = REAL(FLOAT(N_MERGER))
*
*     output, note N_SINGLE is now N_STAR (R.Sp. Oct. 2021).
      if(rank.eq.0) then
          print*,'CUSTOM OUTPUT: TTOT',TTOT,'N_STAR',
     &     N_SINGLE,'N_BINARY',N_BINARY,'N_MERGER',N_MERGER,
     &     ' Output reduced 0/1=',KZ(48)
          call output_single(TTOT,N_SINGLE,KZ(12),KZ(46),KZ(48))
          call output_binary(TTOT,N_BINARY,KZ(12),KZ(46),KZ(48))
          call output_merger(TTOT,N_MERGER,KZ(12),KZ(46),KZ(48))
      end if
*
      if (NXTLEN.eq.NXTLIMIT.or.TTOT.eq.0.D0) then
         ZSM0 = ZSMASS0*MSCALE_OUT
         ZBM0 = ZBMASS0*MSCALE_OUT
*        call global_output(ZSM0,ZBM0,RC,
*    &        N_SINGLE,N_BINARY,N_MERGER,KZ(7),KZ(19))
      END IF
      
      RETURN

      END
      
*----------------------------------------------------------------------------
      subroutine sev_one_star(I,KW,RM,LUM,TE,MC,RCC,M1)
*
*
*     Get stellar evolution parameters for one star
*     ---------------------------------------------
*
      include 'common6.h'
      include 'scale_out.h'
      INTEGER I,KW
      REAL*8  LUMS(10),TSCLS(20),GB(10)
      REAL*8  RM8,LUM8,MC8,RCC8,M18,M0
      REAL*4  RM,LUM,TE,MC,RCC,M1

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

*----------------------------------------------------------------------------
      subroutine custom_store_single(I,IK)
*
*
*     Save one single star in vectors
*     ---------------------------------
*     
      include 'common6.h'
      include 'galaxy.h'
      include 'output_single.h'
      include 'scale_out.h'
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
      REAL*8 F2DOT(3),F3DOT(3),ONE24,ONE120,S,DTR
      PARAMETER  (ONE24=1.0/24.0D0,ONE120=1.0/120.0D0)
C*     single
          S = TIME - T0(I)
          DTR = TIME - T0R(I)
*      Store F2DOT, F3DOT ready for Taylor series with factor as F,FDOT
      DO 30 K = 1,3
          F2DOT(K) = ONE24*(D3R(K,I)*DTR+D2R(K,I) + D3(K,I)*S+D2(K,I))
          F3DOT(K) = ONE120*(D3R(K,I) + D3(K,I))
 30   CONTINUE
*
*     Name and Label for binary members / c.m. (R.Sp. Aug. 2021).
      NS_NAM(IK) = NAME(I)
      NS_LAB(IK) = 0
      IF(I.LT.IFIRST)NS_LAB(IK) = -1
      IF(I.GT.N)NS_LAB(IK) = 1
*     Dynamics
      S_M(IK)  = REAL(BODY(I)*MSCALE_OUT)
      S_RSN(IK) = REAL(RS(I)*RSCALE_OUT)
      S_X1(IK) = REAL(X(1,I)*RSCALE_OUT) 
      S_X2(IK) = REAL(X(2,I)*RSCALE_OUT) 
      S_X3(IK) = REAL(X(3,I)*RSCALE_OUT) 
      S_V1(IK) = REAL(XDOT(1,I)*VSCALE_OUT) 
      S_V2(IK) = REAL(XDOT(2,I)*VSCALE_OUT) 
      S_V3(IK) = REAL(XDOT(3,I)*VSCALE_OUT) 
      S_F1(IK) = REAL(F(1,I)*FSCALE_OUT) 
      S_F2(IK) = REAL(F(2,I)*FSCALE_OUT) 
      S_F3(IK) = REAL(F(3,I)*FSCALE_OUT) 
      S_POT(IK) = REAL(-PHIDBL(I)*POTSCALE_OUT)
      S_FD1(IK) = REAL(FDOT(1,I)*FDSCALE_OUT)
      S_FD2(IK) = REAL(FDOT(2,I)*FDSCALE_OUT)
      S_FD3(IK) = REAL(FDOT(3,I)*FDSCALE_OUT)
      S_D21(IK) = REAL(F2DOT(1)*D2SCALE_OUT)
      S_D22(IK) = REAL(F2DOT(2)*D2SCALE_OUT)
      S_D23(IK) = REAL(F2DOT(3)*D2SCALE_OUT)
      S_D31(IK) = REAL(F3DOT(1)*D3SCALE_OUT)
      S_D32(IK) = REAL(F3DOT(2)*D3SCALE_OUT)
      S_D33(IK) = REAL(F3DOT(3)*D3SCALE_OUT)
      S_STEP(IK) = REAL(STEP(I))
      S_STEPR(IK) = REAL(STEPR(I))
      S_T0(IK) = REAL(T0(I))
      S_T0R(IK) = REAL(T0R(I))

*     if(rank.eq.0)print*,' css i,ik,d2,d3,d2sc,d3sc=',
*    &  i,ik,s_d21(1),s_d31(1),d2scale_out,d3scale_out
*     if(rank.eq.0)print*,' css i,ik,f2,f3=',
*    &  i,ik,f2dot(1),f3dot(1)

*     Stellar evolution of single star
      IF (KZ(12).GT.0) THEN
         call sev_one_star(I,NS_KW(IK),S_RS(IK),S_L(IK),
     &        S_TE(IK),S_MC(IK),S_RC(IK),S_M(IK))
      END IF

      RETURN

      END

*----------------------------------------------------------------------------
      subroutine custom_store_binary(IK,I,NCM,J1,J2,M1,M2,HIP,RIP,RDOT,
     &     XIREL,VIREL,FREL,FDREL,D2REL,D3REL,KC,GI)
*
*
*     Save special binary system
*     ---------------------------------
*     
      include 'common6.h'
      include 'output_binary.h'
      include 'scale_out.h'
      COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)

C*     Binary
      REAL*8 XIREL(3),VIREL(3),FREL(3),FDREL(3),D2REL(3),D3REL(3)
      REAL*8 M1, M2, MTMP,MCM,RDOT,GI
      REAL*8 F2DOT(3),F3DOT(3),ONE24,ONE120,S,DTR
      PARAMETER  (ONE24=1.0/24.0D0,ONE120=1.0/120.0D0)
*     Name
      NB_NAMC(IK) = NCM
      NB_NAM1(IK) = NAME(J1)
      NB_NAM2(IK) = NAME(J2)
*
          S = TIME - T0(I)
          DTR = TIME - T0R(I)
*      Store F2DOT, F3DOT ready for Taylor series with factor as F,FDOT
      DO 30 K = 1,3
          F2DOT(K) = ONE24*(D3R(K,I)*DTR+D2R(K,I) + D3(K,I)*S+D2(K,I))
          F3DOT(K) = ONE120*(D3R(K,I) + D3(K,I))
 30   CONTINUE
*     Dynamic information
      B_M1(IK) = REAL(M1*MSCALE_OUT) 
      B_M2(IK) = REAL(M2*MSCALE_OUT) 
      B_XC1(IK) = REAL(X(1,I)*RSCALE_OUT) 
      B_XC2(IK) = REAL(X(2,I)*RSCALE_OUT) 
      B_XC3(IK) = REAL(X(3,I)*RSCALE_OUT) 
      B_VC1(IK) = REAL(XDOT(1,I)*VSCALE_OUT) 
      B_VC2(IK) = REAL(XDOT(2,I)*VSCALE_OUT) 
      B_VC3(IK) = REAL(XDOT(3,I)*VSCALE_OUT) 
      B_FC1(IK) = REAL(F(1,I)*FSCALE_OUT)
      B_FC2(IK) = REAL(F(2,I)*FSCALE_OUT)
      B_FC3(IK) = REAL(F(3,I)*FSCALE_OUT)
      B_FDC1(IK) = REAL(FDOT(1,I)*FDSCALE_OUT)
      B_FDC2(IK) = REAL(FDOT(2,I)*FDSCALE_OUT)
      B_FDC3(IK) = REAL(FDOT(3,I)*FDSCALE_OUT)
      B_D21(IK) = REAL(F2DOT(1)*D2SCALE_OUT)
      B_D22(IK) = REAL(F2DOT(2)*D2SCALE_OUT)
      B_D23(IK) = REAL(F2DOT(3)*D2SCALE_OUT)
      B_D31(IK) = REAL(F3DOT(1)*D3SCALE_OUT)
      B_D32(IK) = REAL(F3DOT(2)*D3SCALE_OUT)
      B_D33(IK) = REAL(F3DOT(3)*D3SCALE_OUT)
      B_XR1(IK) = REAL(XIREL(1)*RSCALE_OUT) 
      B_XR2(IK) = REAL(XIREL(2)*RSCALE_OUT) 
      B_XR3(IK) = REAL(XIREL(3)*RSCALE_OUT) 
      B_VR1(IK) = REAL(VIREL(1)*VSCALE_OUT) 
      B_VR2(IK) = REAL(VIREL(2)*VSCALE_OUT) 
      B_VR3(IK) = REAL(VIREL(3)*VSCALE_OUT) 
      B_FR1(IK) = REAL(FREL(1)*FSCALE_OUT)
      B_FR2(IK) = REAL(FREL(2)*FSCALE_OUT)
      B_FR3(IK) = REAL(FREL(3)*FSCALE_OUT)
      B_FDR1(IK) = REAL(FDREL(1)*FDSCALE_OUT)
      B_FDR2(IK) = REAL(FDREL(2)*FDSCALE_OUT)
      B_FDR3(IK) = REAL(FDREL(3)*FDSCALE_OUT)
      B_D21R(IK) = REAL(D2REL(1)*D2SCALE_OUT)
      B_D22R(IK) = REAL(D2REL(2)*D2SCALE_OUT)
      B_D23R(IK) = REAL(D2REL(3)*D2SCALE_OUT)
      B_D31R(IK) = REAL(D3REL(1)*D3SCALE_OUT)
      B_D32R(IK) = REAL(D3REL(2)*D3SCALE_OUT)
      B_D33R(IK) = REAL(D3REL(3)*D3SCALE_OUT)
      B_POT(IK) = REAL(-PHIDBL(I)*POTSCALE_OUT)
      B_G(IK)   = REAL(GI)
      B_STEP(IK) = REAL(STEP(I))
      B_STEPR(IK) = REAL(STEPR(I))
      B_T0(IK) = REAL(T0(I))
      B_T0R(IK) = REAL(T0R(I))

*     Stellar evolution of two components
      IF (KZ(12).GT.0) THEN
         MTMP = BODY(J1)
         BODY(J1) = M1
         call sev_one_star(J1,NB_KW1(IK),B_RS1(IK),B_L1(IK),
     &        B_TE1(IK),B_MC1(IK),B_RC1(IK),B_M1(IK))
         BODY(J1) = MTMP
         MTMP = BODY(J2)
         BODY(J2) = M2
         call sev_one_star(J2,NB_KW2(IK),B_RS2(IK),B_L2(IK),
     &        B_TE2(IK),B_MC2(IK),B_RC2(IK),B_M2(IK))
         BODY(J2) = MTMP
      END IF
      NB_KWC(IK) = KC
*     Binary parameters
      MCM = M1 + M2
      SEMI = -0.5*MCM/HIP
      B_ECC(IK) = REAL(SQRT((1.0D0 - RIP/SEMI)**2 + 
     &     RDOT**2/(MCM*SEMI)))
      B_P(IK) = REAL(DAYS_OUT*SEMI*SQRT(ABS(SEMI)/MCM))
      B_A(IK) = REAL(SEMI*RAU_OUT)

      RETURN

      END
