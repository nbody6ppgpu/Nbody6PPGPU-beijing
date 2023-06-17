      SUBROUTINE REDUCE_TSTEP_DRAG(I,XI,VI,ttmp)
*
*
*     Reducing timesteps of particles approaching the AD
*     
*
      INCLUDE 'common6.h'
      REAL*8  Z_NEW_DRAG,R_NEW_DRAG2, hz1, dt_tmp,
     &        xi(3), vi(3) 
*
*    Reduce time steps when particle approaches z=0 plane
*
          dt_tmp = time-t0(i)
*
          Z_NEW_DRAG = XI(3)+DT_TMP*VI(3)
          R_NEW_DRAG2 = XI(1)**2+XI(2)**2
*

      if (r_new_drag2.gt.0.22*0.22*1.2*1.2) goto 10

          IF (R_NEW_DRAG2.GT.R_CRIT*R_CRIT) THEN
              hz1 = HZ*0.22
          ELSE
              hz1 = HZ*0.22*(SQRT(R_NEW_DRAG2)/R_CRIT)
          END IF
*
          IF (ABS(XI(3)).GT.hz1*2.2) THEN
*
*           OPEN (UNIT=97,STATUS='UNKNOWN',FORM='FORMATTED',
*     &           FILE='reduced_steps.dat', ACCESS='APPEND')
*                
             IF (Z_NEW_DRAG*XI(3).LT.0.0) THEN
*                call delay_remove_tlist(i,STEP,DTK)
                 ttmp = 0.125*ttmp
*
*             WRITE(97,99)NAME(i),ttmp,0.125*ttmp,TTOT,(XI(K),k=1,3),hz1
*   99        FORMAT (I7, 8(1X, E20.10))      
*             CALL FLUSH(28)
*
                 step(i) = 0.125*step(i)
*                stepr(i) = 0.125*stepr(i)
*                call delay_store_tlist(i)
             END IF
*              
              IF (ABS(XI(3)).LT.hz1*1.8) then 
*                        call delay_remove_tlist(i,step,dtk)  
                         ttmp = 0.5*ttmp
*              WRITE(97,99)NAME(i),ttmp,0.5*ttmp,TTOT,(XI(K),k=1,3),hz1
*                         step(i) = 0.5*step(i)
*                        stepr(i) = 0.5*stepr(i)
*                        call delay_store_tlist(i)
              ENDIF
          END IF  
*
*       Set new time
*         TIMENW(J) = T0(J) + STEP(J)
*
 10   continue
*     
      RETURN
*
      END
*
