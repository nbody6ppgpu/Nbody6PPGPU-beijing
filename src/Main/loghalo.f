      SUBROUTINE LOGHALO_INIT(I,R_SC,Q_SC,V_zero)
*
*
*     Initialization. External potential: logarithmic dark matter halo (BT2008, eq 2.71a)  
*     force and first derivatives  
*     ---------------------------
*     by Taras Panamarev
*

      INCLUDE 'common6.h'
      REAL*8  fhalo(3),fdhalo(3),r_sc,q_sc,v_zero
     &        ,x_ij,y_ij,z_ij,vx_ij,vy_ij,vz_ij
     &        ,tmp,z_tmp,z2_tmp,r_tmp,r2_tmp
*     &        ,r_sc, q_sc, v_zero 
*
*     XI and VI are corrected positions and velocities 
      x_ij = X(1,i)
      y_ij = X(2,i)
      z_ij = X(3,i)
*
      vx_ij = xdot(1,i)
      vy_ij = xdot(2,i)
      vz_ij = xdot(3,i)
*
      r_tmp = r_sc**2 + x_ij**2 + y_ij**2 + z_ij**2/q_sc**2
      
      tmp = -v_zero**2/r_tmp

      fhalo(1) = x_ij * tmp
      fhalo(2) = y_ij * tmp
      fhalo(3) = z_ij/q_sc**2 * tmp

      tmp = tmp/r_tmp

      fdhalo(1) = tmp * (vx_ij*r_tmp -
     &            2*x_ij*(x_ij*vx_ij+y_ij*vy_ij+z_ij*vz_ij/q_sc**2))

      fdhalo(2) = tmp* (vy_ij*r_tmp -
     &            2*y_ij*(x_ij*vx_ij+y_ij*vy_ij+z_ij*vz_ij/q_sc**2))

      fdhalo(3) = tmp* (vz_ij*r_tmp -
     & 2*z_ij*(x_ij*vx_ij+y_ij*vy_ij+z_ij*vz_ij/q_sc**2))/q_sc**2

*     Total Force Acting on a Star
*
       DO K=1,3
         FI(K,i)=FI(K,i)+fhalo(K)
         D1(K,i)=D1(K,i)+fdhalo(K)
       END DO 
*
      RETURN
*
      END
 
*---------------------------------------------------------------------------------------------------------------
      SUBROUTINE LOGHALO(FIRR,FD,XI,VI,R_SC,Q_SC,V_zero)
*
*
*     External potential: logarithmic dark matter halo (BT2008, eq 2.71a)  
*     force and first derivatives  
*     ---------------------------
*     by Taras Panamarev
*

      INCLUDE 'common6.h'
      REAL*8   XI(3),VI(3),FIRR(3),FD(3),
     &        fhalo(3),fdhalo(3)
     &        ,x_ij,y_ij,z_ij,vx_ij,vy_ij,vz_ij
     &        ,tmp,z_tmp,z2_tmp,r_tmp,r2_tmp
     &        ,r_sc, q_sc, v_zero 
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
      r_tmp = r_sc**2 + x_ij**2 + y_ij**2 + z_ij**2/q_sc**2
      
      tmp = -v_zero**2/r_tmp

      fhalo(1) = x_ij * tmp
      fhalo(2) = y_ij * tmp
      fhalo(3) = z_ij/q_sc**2 * tmp

      tmp = tmp/r_tmp

      fdhalo(1) = tmp * (vx_ij*r_tmp -
     &            2*x_ij*(x_ij*vx_ij+y_ij*vy_ij+z_ij*vz_ij/q_sc**2))

      fdhalo(2) = tmp* (vy_ij*r_tmp -
     &            2*y_ij*(x_ij*vx_ij+y_ij*vy_ij+z_ij*vz_ij/q_sc**2))

      fdhalo(3) = tmp* (vz_ij*r_tmp -
     & 2*z_ij*(x_ij*vx_ij+y_ij*vy_ij+z_ij*vz_ij/q_sc**2))/q_sc**2

*     Total Force Acting on a Star
*
       DO K=1,3
         FIRR(K)=FIRR(K)+fhalo(K)
         FD(K)=FD(K)+fdhalo(K)
       END DO 
*
      RETURN
*
      END
 


