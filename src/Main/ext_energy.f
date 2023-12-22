       SUBROUTINE EXT_ENERGY
*
*
*       Potential energy of the central massive blackhole
*        -------------------------------------------------
*
       INCLUDE 'common6.h'
*
       real*8 eblckhl, e_bulge,e_disk,e_halo,e_loghalo
     &        cmbh, eps_bh, r_tmp1
*     
*
       eblckhl = 0.0d0
       e_bulge = 0.0
       e_disk = 0.0
       e_halo = 0.0
       e_loghalo = 0.0
       e_ext = 0.0

      do i = ifirst, ntot
*       
       R2_tmp1 = x(1,i)**2+x(2,i)**2
*
      if (cmbh.gt.0.0)  
     &  eblckhl = eblckhl-cmbh*body(i)/DSQRT(X(1,I)**2+
     &        X(2,I)**2+X(3,I)**2+EPS_bh)

      if (m_bulge.gt.0.0) 
     &     e_bulge = e_bulge -body(i)*m_bulge/
     &     dsqrt(r2_tmp1+(a_bulge+dsqrt(x(3,i)**2+b_bulge**2))**2) 

       if (m_disk.gt.0.0) 
     &     e_disk = e_disk -body(i)*m_disk/
     &     dsqrt(r2_tmp1+(a_disk+dsqrt(x(3,i)**2+b_disk**2))**2) 

       if (m_halo.gt.0.0) 
     &     e_halo = e_halo -body(i)*m_halo/
     &     dsqrt(r2_tmp1+(a_halo+dsqrt(x(3,i)**2+b_halo**2))**2) 

       if (r_scale.gt.0.0)
     &     e_loghalo = e_loghalo + 0.5*body(i)*v_nod**2*
     &     log(r_scale**2+r2_tmp1+x(3,i)**2/q_scale)
     
      enddo
*
      e_ext = e_bulge + e_disk + e_halo + e_loghalo
*
      return
*
        end

