      SUBROUTINE POTGAL_INIT(I,A_SCALE,B_SCALE,M_SCALE)
*
*     Initialization of the xternal potential: buldge,disk or halo (Miyamoto & Nagai 1975,
*     BT 2008, eq 2.69a)  
*     force and first derivatives  (as in phi-GRAPE.c)
*     ---------------------------
*     by Taras Panamarev
*
      INCLUDE 'common6.h'
      REAL*8  fext(3),fdext(3)
     &        ,x_ij,y_ij,z_ij,vx_ij,vy_ij,vz_ij
     &        ,tmp,z_tmp,z2_tmp,r_tmp,r2_tmp
     &        ,a_scale, b_scale, m_scale 
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
      z2_tmp = z_ij**2+b_scale**2
      z_tmp = sqrt(z2_tmp)

      r2_tmp = x_ij**2 + y_ij**2 + (z_tmp + a_scale)**2
      r_tmp = sqrt(r2_tmp)

*      pot_ext = m_scale/r_tmp
      tmp = m_scale/(r2_tmp*r_tmp)

      fext(1) = -tmp*x_ij
      fext(2) = -tmp*y_ij
      fext(3) = -tmp*z_ij * (z_tmp + a_scale)/z_tmp

      tmp = m_scale / (z_tmp*r2_tmp*r2_tmp*r_tmp)

      fdext(1) = tmp * (- vx_ij*z_tmp*r2_tmp
     &                    + 3.0*x_ij*vx_ij*x_ij*z_tmp
     &                    + 3.0*x_ij*vy_ij*y_ij*z_tmp
     &                    + 3.0*x_ij*vz_ij*z_ij*(z_tmp + a_scale)**2 )

      fdext(2) = tmp * (- vy_ij*z_tmp*r2_tmp
     &                    + 3.0*y_ij*vx_ij*x_ij*z_tmp
     &                    + 3.0*y_ij*vy_ij*y_ij*z_tmp
     &                    + 3.0*y_ij*vz_ij*z_ij*(z_tmp + a_scale)**2 )

      fdext(3) = tmp * (- vz_ij*(z_tmp + a_scale)*
     &           ( x_ij**2*( z2_tmp*z_tmp + a_scale*b_scale**2 ) +
     &             y_ij**2*( z2_tmp*z_tmp + a_scale*b_scale**2 ) -
     &           ( 2.0*a_scale*(z_ij**2 - b_scale**2)*z_tmp +
     &             2.0*z_ij**4 + b_scale**2*z_ij**2 -
     &                 b_scale**2*( a_scale**2 + b_scale**2) ) ) 
     &       + 3.0*vx_ij*x_ij*z_ij*z2_tmp*(z_tmp + a_scale)
     &       + 3.0*vy_ij*y_ij*z_ij*z2_tmp*(z_tmp + a_scale) ) / z2_tmp 
*
*     Total Force Acting on a Star
*
       DO K=1,3
         FI(K,i)=FI(K,i)+fext(K)
         D1(K,i)=d1(K,i)+fdext(K)
       END DO 
*
      RETURN
*
      END
 

*----------------------------------------------------------------------------------------------------------------
*
      SUBROUTINE POTGAL(FIRR,FD,XI,VI,A_SCALE,B_SCALE,M_SCALE)
*
*
*     External potential: buldge,disk or halo (Miyamoto & Nagai 1975,
*     BT 2008, eq 2.69a)  
*     force and first derivatives  (as in phi-GRAPE.c)
*     ---------------------------
*     by Taras Panamarev
*

      INCLUDE 'common6.h'
      REAL*8   XI(3),VI(3),FIRR(3),FD(3),
     &        fext(3),fdext(3)
     &        ,x_ij,y_ij,z_ij,vx_ij,vy_ij,vz_ij
     &        ,tmp,z_tmp,z2_tmp,r_tmp,r2_tmp
     &        ,a_scale, b_scale, m_scale 
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
      z2_tmp = z_ij**2+b_scale**2
      z_tmp = sqrt(z2_tmp)

      r2_tmp = x_ij**2 + y_ij**2 + (z_tmp + a_scale)**2
      r_tmp = sqrt(r2_tmp)

*      pot_ext = m_scale/r_tmp
      tmp = m_scale/(r2_tmp*r_tmp)

      fext(1) = -tmp*x_ij
      fext(2) = -tmp*y_ij
      fext(3) = -tmp*z_ij * (z_tmp + a_scale)/z_tmp

      tmp = m_scale / (z_tmp*r2_tmp*r2_tmp*r_tmp)

      fdext(1) = tmp * (- vx_ij*z_tmp*r2_tmp
     &                    + 3.0*x_ij*vx_ij*x_ij*z_tmp
     &                    + 3.0*x_ij*vy_ij*y_ij*z_tmp
     &                    + 3.0*x_ij*vz_ij*z_ij*(z_tmp + a_scale)**2 )

      fdext(2) = tmp * (- vy_ij*z_tmp*r2_tmp
     &                    + 3.0*y_ij*vx_ij*x_ij*z_tmp
     &                    + 3.0*y_ij*vy_ij*y_ij*z_tmp
     &                    + 3.0*y_ij*vz_ij*z_ij*(z_tmp + a_scale)**2 )

      fdext(3) = tmp * (- vz_ij*(z_tmp + a_scale)*
     &           ( x_ij**2*( z2_tmp*z_tmp + a_scale*b_scale**2 ) +
     &             y_ij**2*( z2_tmp*z_tmp + a_scale*b_scale**2 ) -
     &           ( 2.0*a_scale*(z_ij**2 - b_scale**2)*z_tmp +
     &             2.0*z_ij**4 + b_scale**2*z_ij**2 -
     &                 b_scale**2*( a_scale**2 + b_scale**2) ) ) 
     &       + 3.0*vx_ij*x_ij*z_ij*z2_tmp*(z_tmp + a_scale)
     &       + 3.0*vy_ij*y_ij*z_ij*z2_tmp*(z_tmp + a_scale) ) / z2_tmp 
*
*     Total Force Acting on a Star
*
       DO K=1,3
         FIRR(K)=FIRR(K)+fext(K)
         FD(K)=FD(K)+fdext(K)
       END DO 
*
      RETURN
*
      END
 


