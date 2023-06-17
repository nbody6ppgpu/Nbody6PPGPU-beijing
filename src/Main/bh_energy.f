       SUBROUTINE BH_ENERGY
*
*
*       Potential energy of the central massive blackhole
*        -------------------------------------------------
*
       INCLUDE 'common6.h'
*
       eblckhl = 0.0d0

        do i = ifirst, ntot

       eblckhl = eblckhl-cmbh*body(i)/DSQRT(X(1,I)**2+
     &        X(2,I)**2+X(3,I)**2+EPS_bh)
       
       enddo
*
      return
*
        end
