       SUBROUTINE DRAGFORCE(I,FIRR,FD,XI,VI)
*
*
*       Drag force & first derivative.
*       -------------------------
*
      INCLUDE 'common6.h'
      REAL*8  FIRR(3),FD(3),FDRAG(3),FDDRAG(3),XI(3),VI(3),
     &        R2DISC,R2DENS,RDISC,DDENS,RHODOT,VXRE,VYRE,VZRE,VRE,
     &        CKFDRAG,DVXRE,DVYRE,DVZRE,DVRE,KFDOT,EXCENTR1(NMAX),
     &        AORBIT(NMAX),BORBIT(NMAX),CORBIT(NMAX),INCLIN1(NMAX),
     &        DISCMOM,STARMOM,ANGLE1,INCLIN2(NMAX),RR0_S,
     &        RDOT,Rzero,sigma,R2_DELTA,R2zero,H2Rzero,S,BETA_S,
     &        RRzero, RR2zero, RDISC_3over,RDISC_5over, SQRT_MBH,
     &        alphax, ZETA, RRCRIT_ZETA,R_ZETA,R0_ZETA,
     &        RCRIT_ZETA,R_,Z2, R_ZETA_R, RRCRIT, VXDISC,VYDISC,
     &        r_lim2, h_lim2,
     &        Nscale,st_cross_sec,VRE_ov,
     &        VX2,DVX,VDX,VY2,DVY,VDY,VZ2,DVZ,VDZ,
     &        Rstar2,Qd,Rd,DRAG,CONST,BHLDYN,init_rad2,
     &        VESC4,Lambda,M_cl,Rstar2_avg,M_avg,
     &        ADDRAG(3),ABHL(3),FGEO(3),FOLD,c4,
*** New vector based variables
     &        VREL(3),VRR(3),DVR(3),DDVR(3),DVREL(3),VINV3,VINV4,
*** New variables to define Dyn drag
     &        Omega,c_sound,c_sound_ov,sunm,sunr,
     &        R_min,R_max,R_min_ov,
*** New variables for rotation curve LS
     &        RR,TOTALMASS,VDISC2,HDISC,HDISC2,PREFACTOR,
     &        DSIGMA,SQRT_MTOT,GRAD_P
*     &        ,xold(3), vold(3), vv_new, vv_old, SQRT_VD
*
*         GF: New variables added (from radius to Lambda)
*
      INTEGER I, K, N0
*
*      LOGICAL LPR
*
*     Limiting factors for the accretion disk
      r_lim2 = 0.5*0.5
      h_lim2 = 5.0d-03**2
*
      ttot = time + toff
*
*      t_diss_on = 0.125
      if (ttot.lt.t_diss_on) goto 10

*      if (ttot.eq.0.0) then
*      LPR=NAME(I).EQ.250
*
* GF Aug 2022: New definition of constants using DATA/

      Rzero=0.22
      S=4.0
      BETA_S=0.7
      alphax = 0.75
*
*      R_CRIT = 0.0257314
      R2zero=Rzero*Rzero
      H2Rzero=HZ*HZ*R2zero
*
      N0=4000
      Qd=5.0D0
      sunm = 1.D0/ZMBAR !sunmass in nbody units
      sunr = 1.D0/SU !sun radius in nbody units       
      
      Rstar2_avg=2.7097296576790755*(sunr**2)
      M_avg= 1.D0*sunm

      Nscale=(float(N0)/float(N))*LOG(0.4*N)/LOG(0.4*N0)
*      Qzero= 0.01

      Q_DRAG = Qzero*Nscale 
*
*              Auxiliary Calculations
      R2DISC=XI(1)**2+XI(2)**2
      RR=DSQRT(XI(1)**2+XI(2)**2+XI(3)**2) 
      Z2=XI(3)**2
*
      IF (r2disc.lt.r_lim2.and.z2.lt.h_lim2) then
*
      RDISC=DSQRT(R2DISC)
      R_=DSQRT(R2DISC+Z2)

      sigma=(2.0-ALPHA)*Mdisc/(TWOPI*dsqrt(TWOPI)*HZ*Rzero)
      RDOT=(XI(1)*VI(1)+XI(2)*VI(2))/RDISC
*	
*      New density profile 
*
      IF (RDISC.LT.R_CRIT) THEN 
*
      RRCRIT=RDISC/R_CRIT
      ZETA=1.0
      RRCRIT_ZETA=RRCRIT
      R_ZETA=1.0/R2DISC
      R0_ZETA=1.0/R2zero
      RCRIT_ZETA=1.0/(R_CRIT*R_CRIT)
      PREFACTOR=2.0
*
      ELSE
*
      ZETA=0.0
      RRCRIT_ZETA=1.0
      R_ZETA=1.0
      R0_ZETA=1.0
      RCRIT_ZETA=1.0
      PREFACTOR=3.0
*      
      END IF
*
      R_ZETA_R=R_ZETA/R_
*
      RRzero=RDISC/Rzero
      RR2zero=R2DISC/R2zero
      RR0_S=RRzero**S
      R2_DELTA=RRzero**(-alphax)*EXP(-BETA_S*RR0_S)/RRCRIT_ZETA
*
*   Disk density
      DDENS=sigma*R2_DELTA*DEXP(-Z2/(2.0*H2Rzero*RRCRIT_ZETA**2)) 
*
*   The density derivative
        RHODOT=(-RDOT)*(ZETA+alphax+BETA_S*S*RR0_S)/RDISC
     &       -(VI(3)*XI(3)*R_ZETA-ZETA*(R_ZETA_R)*Z2*RDOT)
     &        /(H2Rzero*RCRIT_ZETA)
*
      RDISC_3over=RDISC**(1.5)
      RDISC_5over=RDISC**(2.5)
      SQRT_MBH=DSQRT(CMBH)
*
*   LS,2022: New rotation curve:
*
*   Enable different additional terms:
*     Old rotation curve      
      IF (irot_opt.eq.1) THEN
      	VXDISC=-SQRT_MBH*XI(2)/RDISC_3over
      	VYDISC=SQRT_MBH*XI(1)/RDISC_3over
	VDISC2=CMBH/RDISC
	TOTALMASS=CMBH
	XSLOPE=0.
*
*     Only with mass profile
      ELSE IF (irot_opt.eq.2) THEN
*   	Calculate enclosed mass LS
      	call bisection(RDISC, I)
*
      	TOTALMASS=CMBH+XINTERMASS(I)
      	SQRT_MTOT=DSQRT(TOTALMASS)
*
	VXDISC=-SQRT_MTOT*XI(2)/RDISC_3over
      	VYDISC=SQRT_MTOT*XI(1)/RDISC_3over
        VDISC2=TOTALMASS/RDISC
*
      ELSE IF (irot_opt.eq.9) THEN
*   	Calculate enclosed mass LS
      	call bisection(RDISC,I)
*
 	TOTALMASS=CMBH+XINTERMASS(I)
        HDISC = HZ*RZero*RRCRIT_ZETA
        HDISC2 = HDISC*HDISC
        DSIGMA = -1/RDISC*alphax
* 	Ignore cut-off - with cut-off see below:
*        DSIGMA = -1/RDISC*(alphax+BETA_S*S*RRzero**S)
	GRAD_P = (1+ HDISC2/RDISC*(
     &         XSLOPE(I)/TOTALMASS
     &         +DSIGMA
     &         -PREFACTOR/RDISC))
        VDISC2=TOTALMASS*GRAD_P/RDISC
        SQRT_VD=DSQRT(ABS(VDISC2/R2DISC))
*   New Velocity of the disc
        VXDISC=-XI(2)*SQRT_VD
        VYDISC=XI(1)*SQRT_VD
      ENDIF
*
*
*   Old Velocity of the disc
*      VXDISC=-SQRT_MBH*XI(2)/RDISC_3over
*      VYDISC=SQRT_MBH*XI(1)/RDISC_3over


* (2)  Compute Angular velocity and gas sound speed

        Omega= SQRT_MBH/SQRT(RDISC_3over)
        c_sound= HZ/Omega

*
*   Relative velocity
*      VXRE=VXDISC-VI(1)
*      VYRE=VYDISC-VI(2)
*      VZRE=-VI(3)
*      VRE=DSQRT(VXRE*VXRE + VYRE*VYRE + VZRE*VZRE)

** GF,2021: Here follows an experimental re-def of relative velocity as
* vectors

      VREL(1)=VXDISC-VI(1)
      VREL(2)=VYDISC-VI(2)
      VREL(3)=-VI(3)
      VRE=DSQRT(VREL(1)*VREL(1)+VREL(2)*VREL(2)+VREL(3)*VREL(3))
*
*
* in the case of BH, VESC= speed of light c
*c is in dimensionless units
      c4 = (1.89D+04)**4
      VINV3=1.D0/VRE**3
      VINV4=1.D0/VRE**4

* GF 2022 Calculate Lambda based on the disk density

*      IF (DDENS.lt.1.0D-10) THEN
         Lambda= 10.0D0
*      ELSE
*         c_sound_ov = 1.0/c_sound
*         Lambda= VRE*c_sound_ov
*      ENDIF

*GF 2022: Set option to have scaling (111) with FGEO, FBHL, FDYN
*scaling (11) for 1 to 1 simulations of sBHs (FBHL & FDYN)
* option (1) for geometric drag only

      IF (idragf.eq.111) THEN

         DRAG=Q_DRAG*DDENS*(TWOPI/2.D0*(Rzero**2)/ZMASS)*RADIUS(I)**2
     &   /(Rstar2_avg)*M_avg/BODY(I)

*         CONST=VESC4*(1.D0 + LOG(Lambda))/Qd

           IF (KSTAR(I).lt.14) THEN !Condition to set c as vesc with sBH
                  VESC4= (2.D0*BODY(I)/RADIUS(I))**2
           ELSE IF (KSTAR(I).eq.14) THEN
                VESC4=c4
           ENDIF

* LS 2023 - calculation with vesc4 moved after definition of vesc4
         CONST=VESC4*(1.D0 + LOG(Lambda))/Qd
         
         BHLDYN= CONST*VINV4
         FOLD= DRAG*VRE
         CKFDRAG=FOLD*(1.D0+BHLDYN)

         FGEO(1:3)=FOLD*VREL(1:3) !needed for overall derivative

      ELSE IF (idragf.eq.11) THEN

         DRAG= TWOPI*2.D0*BODY(I)**2*LOG(Lambda) ! Formula from Ostriker99
         CKFDRAG=DRAG*VINV3*DDENS

      ELSE IF (idragf.eq.1) THEN

         DRAG=Q_DRAG*DDENS*(TWOPI/2.D0*(Rzero**2)/ZMASS)*RADIUS(I)**2
     &   /(Rstar2_avg)*M_avg/BODY(I)
         CKFDRAG=DRAG*VRE


      ENDIF


*    GF:Drag components 

      FDRAG(1:3)=CKFDRAG*VREL(1:3)
*
*
*     The dissipation energy:
*      EDISS1=EDISS1+BODY(I)*(FDRAG(1)*VI(1)+FDRAG(2)*VI(2)+
*     &       FDRAG(3)*VI(3))*STEP(I) 
*
      ediss1=ediss1+0.5*body(i)*((a_drag(1,i)+fdrag(1))*(x(1,i)-x0(1,i))
     &                        +(a_drag(2,i)+fdrag(2))*(x(2,i)-x0(2,i)) 
     &                        +(a_drag(3,i)+fdrag(3))*(x(3,i)-x0(3,i)) )

* 
*
*      do ii=1,3
*         adrag(ii) = fdrag(ii)
*      end do
*      
        DTR = TIME - T0R(I)
*
*
** GF 2022: Re-def of DVRE as vectors
** LS,2022: New rotation curve - irot_opt=1 is Gaias formulation:
*
*     Old rotation curve      
      IF (irot_opt.eq.1) THEN
*
        DVREL(1)=-SQRT_MBH*VI(2)/RDISC_3over
     &      +SQRT_MBH*1.5*XI(2)*RDOT/RDISC_5over
     &      -F(1,I)
        DVREL(2)=-SQRT_MBH*VI(1)/RDISC_3over
     &      +SQRT_MBH*1.5*XI(1)*RDOT/RDISC_5over
     &      -F(2,I)
*
*     Only with mass profile
      ELSE IF (irot_opt.eq.2) THEN
*
        DVREL(1)=-SQRT_MTOT*VI(2)/RDISC_3over
     &      +1.5*SQRT_MTOT*XI(2)*RDOT/RDISC_5over
     &	    +0.5*RDOT*XSLOPE(I)/SQRT_MTOT*VXDISC
     &      -F(1,I)
        DVREL(2)=-SQRT_MTOT*VI(1)/RDISC_3over
     &      +SQRT_MTOT*1.5*XI(1)*RDOT/RDISC_5over
     &	    +0.5*RDOT*XSLOPE(I)/SQRT_MTOT*VYDISC
     &      -F(2,I)
*
      ELSE IF (irot_opt.eq.9) THEN
*
        DVREL(1)= VXDISC/XI(2)*VI(2)
     &      -1.5*VXDISC*RDOT/RDISC
     &	    +0.5*RDOT*XSLOPE(I)/TOTALMASS*VXDISC
     &      +0.5*RDOT*VXDISC/GRAD_P*(2/R2DISC
     &	    -1/R2DISC)*(GRAD_P-1)*RDISC
     &      +0.5*RDOT*VXDISC/GRAD_P*HDISC2/RDISC*(
*     &         XSLOPE(I)*XSLOPE(I)/TOTALMASS/TOTALMASS
     &	    -DSIGMA/RDISC
     &	    -PREFACTOR/R2DISC)
     &	    -F(1,I)
        DVREL(2)= VYDISC/XI(1)*VI(1)
     &      -1.5*VYDISC*RDOT/RDISC
     &	    +0.5*XSLOPE(I)/TOTALMASS*VYDISC
     &      +0.5*RDOT*VYDISC/GRAD_P*(2/R2DISC*ZETA
     &	    -1/R2DISC)*(GRAD_P-1)*RDISC
     &      +0.5*RDOT*VYDISC/GRAD_P*HDISC2/RDISC*(
*     &         XSLOPE(I)*XSLOPE(I)/TOTALMASS/TOTALMASS
     &	    -DSIGMA/RDISC
     &	    -PREFACTOR/R2DISC)
     &      -F(2,I)
*	WRITE(*,*)"DVREL",DVREL(1),DVREL(2)
      ENDIF
*
*
        DVREL(3)=-F(3,I)

        DVRE=(VREL(1)*DVREL(1)+VREL(2)*DVREL(2)+VREL(3)*DVREL(3))/VRE
*
* Old version:
*      DVXRE=-SQRT_MBH*VI(2)/RDISC_3over
*     &      +SQRT_MBH*1.5*XI(2)*RDOT/RDISC_5over
*     &      -F(1,I)
*      DVYRE=-SQRT_MBH*VI(1)/RDISC_3over
*     &      +SQRT_MBH*1.5*XI(1)*RDOT/RDISC_5over
*     &      -F(2,I)
*      DVZRE=-F(3,I)
*     
*      DVRE=(VXRE*DVXRE+VYRE*DVYRE+VZRE*DVZRE)/VRE
*
*      KFDOT=DRAG*M_avg/BODY(I)
*   

**************
* GF 2022 Compute drag force derivative 
             IF (idragf.eq.111) THEN !Scaled version
          
              VRR(1:3) = VREL(1:3)*VINV3
              DVR(1:3) = -3.0 *VREL(1:3)*VINV3/VRE
              DDVR(1:3)= VINV3*DVREL(1:3)

              ADDRAG(1:3)=DRAG*VRE*(RHODOT*VREL(1:3) !FGEO deriv.
     &        +DVRE/VRE*VREL(1:3)+DVREL(1:3))

              ABHL(1:3)= CONST*DDENS*(RHODOT*VRR(1:3) !FBHL & FDYN deriv.
     &         +DVR(1:3)+DDVR(1:3))
              
              FDDRAG(1:3)=ADDRAG(1:3)*(1.0+BHLDYN)  !Overall derivative
     &          +FGEO(1:3)*ABHL(1:3)

             ELSE IF (idragf.eq.11) THEN !No scaling derivative
                 
              VRR(1:3) = VREL(1:3)*VINV3
              DVR(1:3) = -3.0 *VREL(1:3)*VINV3/VRE
              DDVR(1:3)= VINV3*DVREL(1:3)

              FDDRAG(1:3)= DRAG*DDENS*(RHODOT*VRR(1:3)+DVR(1:3)
     &          +DDVR(1:3)) !This is only including FBHL & FDYN on sBHs

             ELSE IF (idragf.eq.1) THEN

              FDDRAG(1:3)=DRAG*VRE*(RHODOT*VREL(1:3) !Fgeo deriv.
     &        +DVRE/VRE*VREL(1:3)+DVREL(1:3))
                 
            ENDIF              
*
*    
*      fdx = dsqrt(fdrag(1)**2+fdrag(2)**2+fdrag(3)**2)
*      fddx = dsqrt(fddrag(1)**2+fddrag(2)**2+fddrag(3)**2)
*      stepdd = fdx/fddx
*      IF(stepdd.lt.step(i))print*,' i,dt,dtd=',i,step(i),stepdd
*
        do k=1,3
           a_drag(k,i) = fdrag(k)
        end do
*
** LS,2023: Save rot. velocity and rel. vel. in vrot.97
*
      IF (mod(TIME,2.).eq.0) THEN
        WRITE (97, 101) ttot,RR,VDISC2,XINTERMASS(I),VRE,RDISC,I,
     &          VREL(1),VREL(2),VREL(3), a_drag(1, I),
     &          a_drag(2, I), a_drag(3, I)
 101    FORMAT (e15.6, e15.6, e15.6, e15.6, e15.6, e15.6, I15.0,
     &          e15.6, e15.6, e15.6,e18.8, e20.8, e18.11)
        CALL FLUSH(27)
      ENDIF
*
*
        ELSE 
        do k=1,3
           a_drag(k,i) = 0.0
           fddrag(k) = 0.0
        end do
* end if (r_lim2 & h_lim2):
      END IF
*
*              Total Force Acting on a Star
      DO K=1,3
         FIRR(K)=FIRR(K)+A_DRAG(K,I)
         FD(K)=FD(K)+FDDRAG(K)
*
*   !!  Try experimental ediss calculation:
*     We need to add drag component to the total force and it's
*     derivative in order to use XVPRED
*         F(K,I) = F(K,I) + FDRAG(K)
*         FDOT(K,I) = FDOT(K,I) + FDDRAG(K)  
*         xold(k) = x(k,i)
*         vold(k) = xdot(k,i)
      END DO 
*
*   Predict coordinates and velocities
*   experimental feature
*      if (kdr.eq.0) then
*      call xvpred(i,nnb)
*      vv_old2 = vold(1)**2+vold(2)**2+vold(3)**2
*      vv_new2 = xdot(1,i)**2+xdot(2,i)**2+xdot(3,i)**2
*      ediss1 = ediss1 + 0.5*body(i)*(vv_new2-vv_old2)      
*         F(K,I) = F(K,I) - FDRAG(K)
*         FDOT(K,I) = FDOT(K,I) - FDDRAG(K)
*      end if
*
*
  10    RETURN
*
      END
