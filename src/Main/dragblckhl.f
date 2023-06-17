        SUBROUTINE DRAGBLCKHL(FIRR,FD,XI,VI)
*
*
*       Black hole force & first derivative.
*       -------------------------
*       Based on Ch. Omarov's 2002 routine and phi-grape.c
*       by Taras Panamarev
*
*
      INCLUDE 'common6.h'
      REAL*8   XI(3),VI(3),FIRR(3),FD(3),
     &        FBLACKH(3),FDBLACKH(3)
     &        ,x_ij,y_ij,z_ij,vx_ij,vy_ij,vz_ij,temp1,rr,r2,rv_ij
     &        ,epsBH 
*
*     XI and VI are corrected positions and velocities 
      x_ij = XI(1)
      y_ij = XI(2)
      z_ij = XI(3)
*
      vx_ij = VI(1)
      vy_ij = VI(2)
      vz_ij = VI(3)
*
*     Black hole softening (take from imput file):
      epsBH = EPS_bh 
      r2 = x_ij**2+y_ij**2+z_ij**2+epsBH**2
      rr = dsqrt(r2)
      rv_ij = x_ij*vx_ij+y_ij*vy_ij+z_ij*vz_ij
      temp1 = cmbh/rr
*      pot_ext(I) = temp1
      temp1 = temp1/r2
*
*     Force from the black Hole:
      FBLACKH(1) = -temp1*x_ij
      FBLACKH(2) = -temp1*y_ij
      FBLACKH(3) = -temp1*z_ij
*    First derivative of the force
      FDBLACKH(1) = -temp1 * (vx_ij - 3.0*rv_ij * x_ij/r2)
      FDBLACKH(2) = -temp1 * (vy_ij - 3.0*rv_ij * y_ij/r2)
      FDBLACKH(3) = -temp1 * (vz_ij - 3.0*rv_ij * z_ij/r2)
*     
*     Total Force Acting on a Star
*
       DO K=1,3
         FIRR(K)=FIRR(K)+FBLACKH(K)
         FD(K)=FD(K)+FDBLACKH(K)
       END DO 
*
      RETURN
*
      END
      
