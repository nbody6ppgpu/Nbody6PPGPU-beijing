      SUBROUTINE INIT_DRAG(I1,I2)
*
*     Initialization of drag force
*     -----------------------------
*
      include 'common6.h'
*
*
      do i = i1,i2
*         
         do k=1,3
            a_drag(k,i) = 0.0
         enddo
*
       enddo

      return
      END 
