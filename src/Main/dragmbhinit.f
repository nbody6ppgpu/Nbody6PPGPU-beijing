      SUBROUTINE DRAGMBHINIT(I1,I2)
*
*
*       Force & first derivative.
*       -------------------------
*
      INCLUDE 'common6.h'
*     Call Drag Force
      DO 90 I=I1,I2
      CALL DRAGBLCKHL1(I)
*      if (Qzero.gt.0.0) CALL DRAGFORCE1(I)
*      DO 95 K = 1,3
*
*        F(K,I)=FI(K,I)+FR(K,I)
*        FDOT(K,I)=D1(K,I)+D1R(K,I)
*
*        D0(K,I) = FI(K,I)
*        D0R(K,I) = FR(K,I)
*        FIDOT(K,I)=D1(K,I)
*        FRDOT(K,I)=D1R(K,I)
*   95 CONTINUE
   90 CONTINUE
*
      RETURN     
      END
