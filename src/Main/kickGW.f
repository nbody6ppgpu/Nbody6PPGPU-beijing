      SUBROUTINE KICKGW(I1,I2,kicktype)
*
*     Written on July 2021 by Manuel Arca Sedda 
*     (the kick.F routine has been used as a basis for this)
*
*     Velocity kick for merging WD, neutron stars or black holes.
*     ---------------------------------------------------
*

      Include 'kspars.h'
      INCLUDE 'common6.h'

      REAL*8 RAN2
      REAL*8 VK(4)
      REAL*8 mbh1,mbh2,abh1,abh2
      
      REAL*8 mas1,mas2,mas3,Qr, acosa, acosb, acosg
      REAL*8 mfinal, sfinal

*      SPIN(I1) = 0.7
*      SPIN(I2) = 0.7

      Qr = BODY(I1)/BODY(I2)
      IF(Qr .GT. 1.0)THEN
         Qr = 1.0 / Qr
      ENDIF
      

      IF(kicktype .EQ. -2)THEN
         VKICK = 0.0
      ELSE IF(kicktype .EQ. -1)THEN
         VKICK = 1.E4         
      ELSE IF(kicktype .EQ. 0)THEN
****  THIS APPROXIMATE THE VKICK - Q PLANE IF SPIN = 0.01 **********         
         mas1 = 212.18D0
         mas2 = 86.7878D0
         mas3 = 0.769689D0         
         VKICK = mas1*Qr**1.2D0/(mas2*Qr**1.2D0+2.0D0) 
         VKICK = VKICK * exp(-mas3*Qr**4.2D0)
         VKICK = 10.D0 ** VKICK
         mfinal= 0.985 * (BODY0(I1)+BODY0(I2)) * ZMBAR
         sfinal= 0.2        
      ELSE IF(kicktype .EQ. 1)THEN
****** adimensional spins:      SPIN(I1), SPIN(I2)         

         !SPIN(I1) = 0.5
         !SPIN(I2) = 0.5

         acosa = -1. + 2.*RAN2(IDUM1)
         acosb = -1. + 2.*RAN2(IDUM1) 
         acosg = -1. + 2.*RAN2(IDUM1)

         VKICK = 0.0
         mfinal= 0.0
         sfinal= 0.0

         mbh1 = BODY(I1)*ZMBAR
         mbh2 = BODY(I2)*ZMBAR
    
         IF(ASPN(I1).LE.0.0.OR.ASPN(I1).GT.1.0)THEN     
           ASPN(I1) = 0.5
           SPIN(I1) = 0.5
         ENDIF
         IF(ASPN(I2).LE.0.0.OR.ASPN(I2).GT.1.0)THEN
           ASPN(I2) = 0.5
           SPIN(I2) = 0.5
         ENDIF   

         abh1 = ASPN(I1)
         abh2 = ASPN(I2)
         
         CALL GWKICKS(mbh1, mbh2 ,abh1 ,abh2, Qr,
     &              acosa, acosb, acosg, VKICK, mfinal, sfinal)
         
         ASPN(I1) = sfinal
         ASPN(I2) = sfinal
         SPIN(I1) = sfinal
         SPIN(I2) = sfinal

      ENDIF
      
      !mfinal = mfinal / ZMBAR
      !BODY(I1) = 0.0
      !BODY(I2) = mfinal
      
      

      THETA = RAN2(IDUM1)*TWOPI
      PHI   = (RAN2(IDUM1)-0.5D0)*TWOPI/2.D0
      SPHI  = SIN(PHI)
      X1    = ASIN(SPHI)
      CPHI  = COS(X1)
      VK(1) = COS(THETA)*CPHI*VKICK
      VK(2) = SIN(THETA)*CPHI*VKICK
      VK(3) = SPHI*VKICK
      VK(4) = VKICK
      
      IF(rank.eq.0) THEN
         
         WRITE (6,310)  TTOT*TSTAR, NAME(I1), NAME(I2), KSTAR(I1), 
     &        KSTAR(I2), BODY(I1)*ZMBAR, BODY(I2)*ZMBAR,
     &        ASPN(I1),ASPN(I2), Qr,
     &        VKICK, VK(1),VK(2),VK(3)
 310       FORMAT (' GW recoil:  TIME[Myr] NAME(I1) NAME(I2) K*(I1) ',
     &        'K*(I2) M(I1)[M*] M(I2)[M*] SPIN(I1) SPIN(I2) Q', 
     &        ' VGW[km/s] (V,vx,vy,vz)',
     &        1P,E14.3,0P,2I10,2I4,5F9.3,
     &        1P,4E14.3)
         
         NBKICK = NBKICK + 1
         
      ENDIF

      
            
      I = I2           
      VI2 = 0.0
      VF2 = 0.0
*      CALL JPRED(I,TIME,TIME)
      DO 10 K = 1,3
          VI2 = VI2 + XDOT(K,I)**2
          XDOT(K,I) = XDOT(K,I) + VK(K)/VSTAR
          X0DOT(K,I) = XDOT(K,I)
          VF2 = VF2 + XDOT(K,I)**2
 10   CONTINUE
*     Modify energy loss due to increased velocity of single particle.
      DETMP = - 0.5*BODY(I)*(VF2 - VI2)

      I = I1            
      VI2 = 0.0
      VF2 = 0.0
*     CALL JPRED(I,TIME,TIME)
      DO 11 K = 1,3
          VI2 = VI2 + XDOT(K,I)**2
          XDOT(K,I) = XDOT(K,I) + VK(K)/VSTAR
          X0DOT(K,I) = XDOT(K,I)
          VF2 = VF2 + XDOT(K,I)**2
 11   CONTINUE
*     Modify energy loss due to increased velocity of single particle.
      DETMP = DETMP - 0.5*BODY(I)*(VF2 - VI2)


      ECDOT = ECDOT + DETMP
      
*     ks MPI communication
*** removed the following call by M.A.S. Manuel Arca Sedda 12 Nov 2021
*     call ksparmpi(K_store,K_real8,K_ECDOT,0,0,DETMP)


      NKICK = NKICK + 1

      
      RETURN
*
      END



      

      SUBROUTINE GWKICKS(m1, m2, a1, a2, q, cosa, cosb, cosg,
     &  vkick,mfina,sfina)
      
******* CALCULATION OF KICKS ACCORDING TO CAMPANELLI, LOUSTO, ZLOCHOWER, ETC.. SEE ARCA SEDDA ET AL (2020) ***************
      REAL *8 m1, m2, a1, a2, q, cosa, cosb, cosg
      REAL *8 sina,sinb,sing
      REAL *8 a2par,a2per1,a2per2
      REAL *8 a1par,a1per1,a1per2
      REAL *8 Deltapar,Deltaper1,Deltaper2,Deltaper
      REAL *8 dir1,pm,dir2,dir
      REAL *8 SEpar, PHI, PHI1, A, B, H, XI
      REAL *8 V11, VA, VB, VC, ETA, VM, VPER
      REAL *8 VTER,SPPE,VPAR,vkiper1,vkiper2,vkipar,vkick
      REAL *8 rnd
      REAL *8 mfina, sfina
      
      PARAMETER(M_PI=3.1415)
      
      sina = sin(acos(cosa))
      sinb = sin(acos(cosb))
      sing = sin(acos(cosg))
      
      
      a2par = a2 * cosg
      a2per1= a2 * sing 
      a2per2= 0.0
      
      a1par = a1 * cosb
      a1per1= a1 * sinb*cosa
      a1per2= a1 * sinb*sina
      
      Deltapar  = (m1+m2)*(m1+m2) / (1+q) * (a2par  - q*a1par)
      Deltaper1 = (m1+m2)*(m1+m2) / (1+q) * (a2per1 - q*a1per1)
      Deltaper2 = (m1+m2)*(m1+m2) / (1+q) * (a2per2 - q*a1per2)
      
      Deltaper  = sqrt(Deltaper1*Deltaper1 + Deltaper2*Deltaper2)

      CALL RANDOM_NUMBER(rnd)
      dir1 = rnd
      CALL RANDOM_NUMBER(rnd)
      pm   = 1.-2.*rnd
      dir2 = abs(pm)/pm * sqrt(1.-dir1*dir1)
      dir  = sqrt(dir1*dir1 + dir2*dir2)
      
      
      SEpar = 2.*(a2par + q*q*a1par)/((1.+q)*(1.+q))
      
      
      PHI = (Deltaper1*dir1 + Deltaper2*dir2)/(Deltaper * dir)
      CALL RANDOM_NUMBER(rnd)
      PHI1= 2.*M_PI*rnd
      A = 1.2E4 ! km/s  
      B = -0.93 ! adim
      H = 6.9E3 !km/s
      XI= 145.*M_PI/180. !degrees 
      V11 = 3677.76 ! km/s
      VA  = 2481.21 ! km/s
      VB  = 1792.45 ! km/s
      VC  = 1506.52 ! km/s
      
      ETA= q/((1.+q)*(1.+q)) !asimmetric mass ratio
      
      
      VM   = A*ETA*ETA*sqrt(1.-4.*ETA)*(1.+B*ETA)      
      
      VPER = H*ETA*ETA/(1.+q) * (a2par - q*a1par)
      
      VTER = V11 + VA*SEpar + VB*SEpar*SEpar + VC*SEpar*SEpar*SEpar
      SPPE = sqrt((a2per1 - q*a1per1)*(a2per1 - q*a1per1)
     &       + (a2per2 - q*a1per2)*(a2per2 - q*a1per2))*cos(PHI-PHI1)
      VPAR = 16.*ETA*ETA/(1.+q) * VTER * SPPE
      
      vkiper1 = VM + VPER*cos(XI)
      vkiper2 = VPER*sin(XI)
      vkipar  = VPAR
      
      vkick = sqrt(vkiper1*vkiper1 + vkiper2*vkiper2 + vkipar*vkipar)

      mfina = 0.0
      sfina = 0.0
      
      CALL bbh_final_mass_non_precessing_UIB2016(m1,m2,a1,a2,mfina)
      CALL bbh_final_spin_non_precessing_UIB2016(m1,m2,
     &                               a1par,a2par,sfina)


      
      END SUBROUTINE

      
!      FUNCTION pow(A,B)
!      REAL *8 power
!      power = A**B     
!      END FUNCTION POW
      


      SUBROUTINE bbh_UIBfits_setup(mbh1, mbh2,
     &     chibh1, chibh2, out)
      REAL *8 m1, m2, chi1, chi2, m, msq, m1sq, m2sq
      REAL *8 eta, eta2, eta3, eta4
      REAL *8 S1, S2,Stot, Shat, Shat2, Shat3, Shat4
      REAL *8 chidiff, chidiff2, sqrt2, sqrt3, sqrt1m4eta
      REAL *8 out(15)
      REAL *8 mbh1,mbh2,chibh1,chibh2
      

      m1 = mbh1
      m2 = mbh2
      chi1=chibh1
      chi2=chibh2

      IF(m2.GT.m1)THEN
         m2 = mbh1
         chi2 = chibh1
         m1 = mbh2
         chi1 = chibh2
      ENDIF

      m    = m1+m2
      
      msq  = m*m
      m1sq = m1*m1
      m2sq = m2*m2
      
      eta  = m1*m2/msq
      if(eta>0.25)eta = min(eta,0.25d0)
      if(eta<0.0) eta = 0.0
      
      eta2 = eta*eta
      eta3 = eta2*eta
      eta4 = eta2*eta2
      
!     spin variables (in m = 1 units)
      S1    = chi1*m1sq/msq     !# spin angular momentum 1
      S2    = chi2*m2sq/msq     !# spin angular momentum 2
      Stot  = S1+S2             !# total spin
      Shat  = (chi1*m1sq+chi2*m2sq)/(m1sq+m2sq) !# effective spin, = msq*Stot/(m1sq+m2sq)
      Shat2 = Shat*Shat
      Shat3 = Shat2*Shat
      Shat4 = Shat2*Shat2
      
      chidiff  = chi1 - chi2
      chidiff2 = chidiff*chidiff
      
!    # typical squareroots and functions of eta
      sqrt2 = 2.**0.5
      sqrt3 = 3.**0.5
      
      sqrt1m4eta = (1. - 4.*eta**0.5)
      
      out(1) = m
      out(2) = eta
      out(3) = eta2
      out(4) = eta3
      out(5) = eta4
      out(6) = Stot
      out(7) = Shat
      out(8) = Shat2
      out(9) = Shat3
      out(10) = Shat4
      out(11) = chidiff
      out(12) = chidiff2
      out(13) = sqrt2
      out(14) = sqrt3
      out(15) = sqrt1m4eta
      

      END SUBROUTINE




      SUBROUTINE bbh_final_mass_non_precessing_UIB2016(m1, m2,
     &                                     chi1, chi2,mfin)

!## double Functions::bbh_final_mass_non_precessing_UIB2016(double m1, double m2, double chi1, double chi2, string version){


      REAL *8 out(15)
      REAL *8 m1,m2,chi1,chi2,mfina
      REAL *8 m, eta, eta2, eta3, eta4
      REAL *8 Stot, Shat, Shat2, Shat3, Shat4
      REAL *8 chidiff, chidiff2, sqrt2, sqrt3, sqrt1m4eta 
      REAL *8 b10, b20, b30, b50
      REAL *8 a2, a3, a4, b1, b2, b3, b5
      REAL *8 f20, f30, f50, f10, f21
      REAL *8 d10, d11, d20, d30, d31, f11, f31, f51
      REAL *8 Erad, Mf, mfin
      REAL *8 Erad1, Erad2, Erad3, Erad4, Erad5, Erad6, Erad7, Erad8
      REAL *8 Erad_den, Erad_num
 

      out(1) = 0.0
      out(2) = 0.0
      out(3) = 0.0
      out(4) = 0.0
      out(5) = 0.0
      out(6) = 0.0
      out(7) = 0.0
      out(8) = 0.0
      out(9) = 0.0
      out(10) =0.0
      out(11) =0.0
      out(12) =0.0
      out(13) =0.0
      out(14) =0.0
      out(15) =0.0
      
      CALL bbh_UIBfits_setup(m1, m2, chi1, chi2, out)

      m          = out(1)
      eta        = out(2)
      eta2       = out(3)
      eta3       = out(4)
      eta4       = out(5)
      Stot       = out(6)
      Shat       = out(7)
      Shat2      = out(8)
      Shat3      = out(9)
      Shat4      = out(10)
      chidiff    = out(11) 
      chidiff2   = out(12)
      sqrt2      = out(13)
      sqrt3      = out(14)
      sqrt1m4eta = out(15)

     
     
!# rational-function Pade coefficients (exact) from Eq. (22) of 1611.00332v2
      b10 = 0.346;
      b20 = 0.211;
      b30 = 0.128;
      b50 = -0.212;
!     # fit coefficients from Tables VII-X of 1611.00332v2
!     # values at increased numerical precision copied from
!     # https://git.ligo.org/uib-papers/finalstate2016/blob/master/LALInference/EradUIB2016v2_pyform_coeffs.txt
!     # git commit f490774d3593adff5bb09ae26b7efc6deab76a42
      a2 = 0.5609904135313374;
      a3 = -0.84667563764404;
      a4 = 3.145145224278187;
      b1 = -0.2091189048177395;
      b2 = -0.19709136361080587;
      b3 = -0.1588185739358418;
      b5 = 2.9852925538232014;
      f20 = 4.271313308472851;
      f30 = 31.08987570280556;
      f50 = 1.5673498395263061;
      f10 = 1.8083565298668276;
      f21 = 0.;
      d10 = -0.09803730445895877;
      d11 = -3.2283713377939134;
      d20 = 0.01118530335431078;
      d30 = -0.01978238971523653;
      d31 = -4.91667749015812;
      f11 = 15.738082204419655;
      f31 = -243.6299258830685;
      f51 = -0.5808669012986468;
      

      Erad1=((1. - 2.0/3.0*sqrt2)*eta + a2*eta2 + a3*eta3 + a4*eta4)
      Erad2=b10*b1*Shat*(f10 + f11*eta + (16. - 16.*f10 - 4.*f11)*eta2)
      Erad3=b20*b2*Shat2*(f20 + f21*eta + (16. - 16.*f20 - 4.*f21)*eta2)
      Erad4=b30*b3*Shat3*(f30 + f31*eta + (16. - 16.*f30 - 4.*f31)*eta2)
      Erad5=b50*b5*Shat*(f50 + f51*eta + (16. - 16.*f50 - 4.*f51)*eta2)
      Erad6=d10*sqrt1m4eta*eta2*(1. + d11*eta)*chidiff
      Erad7=d30*Shat*sqrt1m4eta*eta*(1. + d31*eta)*chidiff
      Erad8=d20*eta3*chidiff2
      
      Erad_num = (Erad1 * ( 1.d0 + Erad2 + Erad3 + Erad4 ))
      Erad_den = (1. + Erad5)
      
      Erad = Erad_num / Erad_den + Erad6 + Erad7 + Erad8 

      
!     # Convert to actual final mass
      Mf = m*(1.-Erad);
      mfin = Mf;



      END SUBROUTINE


      SUBROUTINE bbh_final_spin_non_precessing_UIB2016(m1, m2,
     &                                              chi1, chi2,sfin)
  
!     # double Functions::bbh_final_spin_non_precessing_UIB2016(double m1, double m2, double chi1, double chi2, string version){
  
      REAL *8 out(15)

      REAL *8 m, eta, eta2, eta3, eta4, m1, m2, chi1, chi2
      REAL *8 sfina, Stot, Shat, Shat2, Shat3, Shat4
      REAL *8 chidiff, chidiff2, sqrt2, sqrt3, sqrt1m4eta

      REAL *8 sfin
      REAL *8 a20, a30, a50, b10, b20, b30, b50

      REAL *8 a2, a3, a5, b1, b2, b3, b5
      REAL *8 f21, f31, f50, f11, f52
      REAL *8 d10, d11, d20, d30, d31, f12, f22, f32, f51
      REAL *8 Lorb, chif
      REAL *8 Lorb_num,Lorb_den,Lorb4_pre,Lorb4_ins
      REAL *8 Lorb1,Lorb2,Lorb3,Lorb4,Lorb5,Lorb6,Lorb7
       

      out(1) = 0.0
      out(2) = 0.0
      out(3) = 0.0
      out(4) = 0.0
      out(5) = 0.0
      out(6) = 0.0
      out(7) = 0.0
      out(8) = 0.0
      out(9) = 0.0
      out(10) =0.0
      out(11) =0.0
      out(12) =0.0
      out(13) =0.0
      out(14) =0.0
      out(15) =0.0
      

      CALL bbh_UIBfits_setup(m1, m2, chi1, chi2, out)

      m          = out(1);
      eta        = out(2);
      eta2       = out(3);
      eta3       = out(4);
      eta4       = out(5);
      Stot       = out(6);
      Shat       = out(7);
      Shat2      = out(8);
      Shat3      = out(9);
      Shat4      = out(10);
      chidiff    = out(11); 
      chidiff2   = out(12); 
      sqrt2      = out(13);
      sqrt3      = out(14);
      sqrt1m4eta = out(15);
      
!     # rational-function Pade coefficients (exact) from Eqs. (7) and (8) of 1611.00332v2
      a20 = 5.24;
      a30 = 1.3;
      a50 = 2.88;
      b10 = -0.194;
      b20 = 0.0851;
      b30 = 0.00954;
      b50 = -0.579;
!     # fit coefficients from Tables I-IV of 1611.00332v2
!     # values at increased numerical precision copied from
!     # https://git.ligo.org/uib-papers/finalstate2016/blob/master/LALInference/FinalSpinUIB2016v2_pyform_coeffs.txt
!     # git commit f490774d3593adff5bb09ae26b7efc6deab76a42
      a2 = 3.8326341618708577;
      a3 = -9.487364155598392;
      a5 = 2.5134875145648374;
      b1 = 1.0009563702914628;
      b2 = 0.7877509372255369;
      b3 = 0.6540138407185817;
      b5 = 0.8396665722805308;
      f21 = 8.77367320110712;
      f31 = 22.830033250479833;
      f50 = 1.8804718791591157;
      f11 = 4.409160174224525;
      f52 = 0.;
      d10 = 0.3223660562764661;
      d11 = 9.332575956437443;
      d20 = -0.059808322561702126;
      d30 = 2.3170397514509933;
      d31 = -3.2624649875884852;
      f12 = 0.5118334706832706;
      f22 = -32.060648277652994;
      f32 = -153.83722669033995;
      f51 = -4.770246856212403;
      
      Lorb_num = (2.*sqrt3*eta + a20*a2*eta2 + a30*a3*eta3)
      Lorb_den = (1. + a50*a5*eta)
      
      Lorb1=b10*b1*Shat*(f11*eta+f12*eta2+(64.-16.*f11-4.*f12)*eta3)
      Lorb2=b20*b2*Shat2*(f21*eta+f22*eta2+(64.-16.*f21-4.*f22)*eta3)
      Lorb3=b30*b3*Shat3*(f31*eta+f32*eta2+(64.-16.*f31-4.*f32)*eta3)
      Lorb4_pre = b50*b5*Shat
      Lorb4_ins = (64.-64.*f50-16.*f51-4.*f52)
      Lorb4=Lorb4_pre*(f50+f51*eta+f52*eta2+Lorb4_ins*eta3)
      Lorb5=d10*sqrt1m4eta*eta2*(1. + d11*eta)*chidiff
      Lorb6=d30*Shat*sqrt1m4eta*eta3*(1. + d31*eta)*chidiff
      Lorb7=d20*eta3*chidiff2
      
      Lorb = Lorb_num/Lorb_den
      Lorb = Lorb + (Lorb1+Lorb2+Lorb3)/(1.0+Lorb4)
      Lorb = Lorb +Lorb5+Lorb6+Lorb7
      
      
!     # Convert to actual final spin      
      chif = Lorb + Stot;
      
      sfin = chif;
      
      
      
      END SUBROUTINE
