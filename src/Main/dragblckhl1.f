        SUBROUTINE DRAGBLCKHL1(I)
*
*
*       Initialization of Black hole force & first derivative.
*       -------------------------
*       Based on Ch. Omarov's 2002 routine and phi-grape.c
*       by Taras Panamarev
*
      INCLUDE 'common6.h'
      REAL*8  FBLACKH(3),FDBLACKH(3), firr(3,i), fd(3,i)
     &        ,x_ij,y_ij,z_ij,vx_ij,vy_ij,vz_ij,rr,r2,rv_ij,temp1
     &       ,epsBH 
*
      x_ij = X(1,I)
      y_ij = X(2,I)
      z_ij = X(3,I)
*
      vx_ij = XDOT(1,I)
      vy_ij = XDOT(2,I)
      vz_ij = XDOT(3,I)
*
*     Black hole softening (take from input file)
      epsBH = EPS_bh 
      r2 = x_ij**2+y_ij**2+z_ij**2+epsBH**2
      rr = dsqrt(r2)
      rv_ij = x_ij*vx_ij + y_ij*vy_ij + z_ij*vz_ij
*
      temp1 = CMBH/rr
*      pot_ext(I) = temp1
*
      temp1 = temp1/r2
*
*    Force from the Black Hole:
      FBLACKH(1) = -temp1*x_ij
      FBLACKH(2) = -temp1*y_ij
      FBLACKH(3) = -temp1*z_ij
*    First derivative of the BH force
      FDBLACKH(1) = -temp1 * (vx_ij - 3.0*rv_ij * x_ij/r2)
      FDBLACKH(2) = -temp1 * (vy_ij - 3.0*rv_ij * y_ij/r2)
      FDBLACKH(3) = -temp1 * (vz_ij - 3.0*rv_ij * z_ij/r2)
*     
*      Add BH force to the total force
*
       DO K=1,3
         FI(K,I)=FI(K,I)+FBLACKH(K)
         D1(K,I)=D1(K,I)+FDBLACKH(K)
       END DO 
*
      RETURN
*
      END
