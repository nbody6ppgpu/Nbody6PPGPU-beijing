      subroutine rmesc_tlist(I,B_FLAG)
*
*
*     Remove particle in NXTLST for escape
*
*     B_FLAG: binary remove flag

      include 'params.h'
      include 'mpi_base.h'
      include 'tlist.h'
      INTEGER L,J,LK,K
      LOGICAL RM_FLAG,B_FLAG

*     Remove flag.
      RM_FLAG = .false.

*     Step level tracer
      K = NDTMAX
      L = 1
*     Here avoid use DO structure since NXTLIMIT is not constant
 1    IF(L.LE.NXTLIMIT) THEN
         J = NXTLST(L)
 10      IF(L.GT.NDTK(K)) THEN
            K = K - 1
            GO TO 10
         END IF
         IF(J.GT.I) THEN
            NXTLST(L) = NXTLST(L) - 1
         ELSE IF(J.EQ.I) THEN
*     decrease nxtlst ending point
*           print*,'REMOVE I',I,'NXTLIMIT',NXTLIMIT,'L',L,'K',K,
*    &           'NDTMIN',NDTMIN,'NDTK(K)',NDTK(K),'NDTMAX',NDTMAX,
*    &           'B_FLAG',B_FLAG
*     --07/08/14 16:49-lwang-end----------------------------------------*
            IF(NXTLIMIT.LE.1) THEN
               if(rank.eq.0)write(6,*) 'Error: No particle in NXTLST!'
               call flush(6)
               call abort()
            END IF
            NXTLIMIT = NXTLIMIT - 1
*     Reduce all NDTK outside NDTMIN by one
            DO LK = 1,NDTMIN-1
               NDTK(LK) = NDTK(LK) - 1
            END DO
*     Replace removing particle index by the end index of step level K
            IF(L.NE.NDTK(K)) THEN
               NXTLST(L) = NXTLST(NDTK(K))
               IF(NXTLST(L).GT.I) NXTLST(L) = NXTLST(L) - 1
C*     Avoid reduce index in position NXTLIMIT (Wrong)
C            ELSE IF(K.GT.NDTMIN.AND.NXTLST(NDTK(K-1)).GT.I) THEN
C               NXTLST(NDTK(K-1)) = NXTLST(NDTK(K-1)) - 1
            END IF
            NDTK(K) = NDTK(K) - 1
            DO LK = K-1,NDTMIN,-1
*     Shift last index position to beginning of level L
               IF(NDTK(LK+1)+1.NE.LK) THEN
                  NXTLST(NDTK(LK+1)+1) = NXTLST(NDTK(LK))
*     Avoid the missing -1 if the removed one is the ending of level
                  IF(L.EQ.NDTK(LK+1)+1.AND.NXTLST(L).GT.I) THEN
                     NXTLST(L) = NXTLST(L)-1
                  END IF
               END IF 
*     Reduce step level L position by one
               NDTK(LK) = NDTK(LK) - 1
            END DO
            RM_FLAG = .true.
         END IF
         L = L + 1
         GO TO 1
      END IF

*     For ghost list
      IF(NGHOSTS.GT.0) THEN
         IF(RM_FLAG) NXTLST(NXTLIMIT+1) = NXTLST(NXTLIMIT+1+NGHOSTS)
         L = NXTLIMIT+1
 2       IF(NXTLST(L).EQ.I) THEN
            IF(RM_FLAG) THEN
       if(rank.eq.0)
     &      write(6,*) 'Error: Particle I',I,'exist in both NXTLST ',
     &              'and Ghost list! L',L
               call flush(6)
               call abort()
            END IF
            NXTLST(L) = NXTLST(NXTLIMIT+NGHOSTS)
            NGHOSTS = NGHOSTS - 1
            RM_FLAG = .true.
         END IF
         IF(NXTLST(L).GT.I) NXTLST(L) = NXTLST(L) - 1
         L = L + 1
         IF(L.LE.NXTLIMIT+NGHOSTS) GO TO 2
      END IF

*     For NLSTDELAY
      IF(NLSTDELAY(1).GT.0) THEN
         L = 2
 3       IF(NLSTDELAY(L).EQ.I) THEN
            IF(RM_FLAG) THEN
               if(rank.eq.0)
     &         write(6,*) 'Error: Particle I',I,'exist in both NXTLST ',
     &              'and NLSTDELAY! L',L
               call flush(6)
               call abort()
            END IF
            NLSTDELAY(L) = NLSTDELAY(NLSTDELAY(1)+1)
            NLSTDELAY(1) = NLSTDELAY(1) - 1
            RM_FLAG = .true.
         END IF
         IF(NLSTDELAY(L).GT.I) NLSTDELAY(L) = NLSTDELAY(L) - 1
         L = L + 1
         IF(L.LE.NLSTDELAY(1)+1) GO TO 3
      END IF

*     Shift all index by 2 due to KS components removed
      IF(B_FLAG) THEN
         DO LK = 1, NXTLIMIT+NGHOSTS
            NXTLST(LK) = NXTLST(LK) - 2
         END DO
         IF(NLSTDELAY(1).GT.0) THEN
            DO LK = 2, NLSTDELAY(1)+1
               NLSTDELAY(LK) = NLSTDELAY(LK) - 2
            END DO
         END IF
      END IF

*     IF(rank.eq.0.)THEN
*        write(6,*)rank,' End rmesc RM_FLAG,B_FLAG,I=',
*    &    RM_FLAG,B_FLAG,I,
*    &   ' nxtlim,ngh,l,k,ndtmin,max,ndtk(min,k,max) ',NXTLIMIT,NGHOSTS,
*    &   L,K,NDTMIN,NDTMAX,NDTK(NDTMIN),NDTK(K),NDTK(NDTMAX)
*        call flush(6)
*     END IF
*
      call shrink_tlist

      RETURN

      END
