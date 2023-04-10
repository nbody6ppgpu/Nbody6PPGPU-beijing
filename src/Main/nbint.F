       SUBROUTINE NBINT(NPART_LEN,NXTLST_LOC,FIRR,FD,IMINR_LOC)
*
*
*       Regular neighbor integration for GPU version
*     --------------------
*
      include 'common6.h'
      include 'tlist.h'
      include 'omp_lib.h'
      
      REAL*8  FIRR(3,NMAX),FD(3,NMAX)
      REAL*8  RIJMIN,XI(3),XIDOT(3)
      INTEGER II,I,NNB,L,J,NXTLST_LOC(NPART_LEN),IMINR_LOC(NMAX)
      REAL*8  A1,A2,A3,DV(3),RIJ2,DR3I,DRDV,DRDP

!$omp parallel do 
!$omp& private(II,I,NNB,L,J,A1,A2,A3,DV,RIJ2,DR2I,DR3I,DRDV,
!$omp&  DRDP,RIJMIN)
      DO II = 1, NPART_LEN
         I = NXTLST_LOC(II)
         ttime = time
         call jpred_int(I,ttime)
*     Initialize arrays and use neighbor list from GPU
         NNB = LIST(1,I)         
         FIRR(1:3,II) = 0.0D0
         FD(1:3,II) = 0.0D0
         IF(NNB.EQ.0) IMINR_LOC(I) = -1
         RIJMIN = 1.E20
*     Perform only extended neighbour loop if GPU has been used.
         DO L = 2, NNB + 1
            J = LIST(L,I)
            ttime = time
            call jpred_int(J,ttime)
            A1 = X(1,J) - X(1,I)
            A2 = X(2,J) - X(2,I)
            A3 = X(3,J) - X(3,I)
*     Predicted coordinates avoids spurious force differences.
            DV(1) = XDOT(1,J) - XDOT(1,I)
            DV(2) = XDOT(2,J) - XDOT(2,I)
            DV(3) = XDOT(3,J) - XDOT(3,I)
            RIJ2 = A1*A1 + A2*A2 + A3*A3
            IF(RIJ2.LT.RIJMIN) THEN
               RIJMIN = RIJ2
               IMINR_LOC(II) = J
            END IF
            DR2I = 1.0/RIJ2
            DR3I = BODY(J)*DR2I*SQRT(DR2I)
            DRDV = A1*DV(1) + A2*DV(2) + A3*DV(3)
            DRDP = 3.0*DRDV*DR2I
*     
            FIRR(1,II) = FIRR(1,II) + A1*DR3I
            FIRR(2,II) = FIRR(2,II) + A2*DR3I
            FIRR(3,II) = FIRR(3,II) + A3*DR3I
            FD(1,II) = FD(1,II) + (DV(1) - A1*DRDP)*DR3I
            FD(2,II) = FD(2,II) + (DV(2) - A2*DRDP)*DR3I
            FD(3,II) = FD(3,II) + (DV(3) - A3*DRDP)*DR3I
         END DO
      END DO
!$omp end parallel do

      RETURN

      END
      
