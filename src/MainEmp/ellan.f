
C***********************************************************************
C
C
                        SUBROUTINE ellan
C
C
C***********************************************************************
C
C
C     Subroutine to analyze the ellipticity of the system
C     Written by Christian Theis 199x; updated RSp Aug. 2021
C
C=======================================================================

        INCLUDE 'common6.h'
        COMMON/POTDEN/  RHO(NMAX),XNDBL(NMAX),PHIDBL(NMAX)
C   Declaration of local variables.
C   -------------------------------
        DOUBLE PRECISION etoti(1:nmax),
     &     ti(3,3),tiwork(3,3),dwork(3),ework(3),lam(3)
        INTEGER i,j,k,ief,nbound,nstart,nnext,np,indexev(3)
        INTEGER index(1:nmax)

        INTEGER nef,nef1,nnpart
        PARAMETER(nef=17,nef1=nef+1)
        INTEGER npart(nef1),ncum(nef1)
        DOUBLE PRECISION xf(nef1),ba(nef1),ca(nef1),taue(nef1),
     &     xmass(nef1),xmcum(nef1),r3av(nef1),r2av(nef1),zav(nef1),
     &     avmass(nef1),vrav(nef1),vzav(nef1),vrot(nef1),sigr2(nef1),
     &     sigph2(nef1),sigz2(nef1),evec(3,3,nef1),xpo,ypo,zpo
        DOUBLE PRECISION vrad(nmax),vtan(nmax)
*
        PARAMETER (tiny2=1.D-30)
      DATA XF/0.001D0,0.003D0,0.005D0,0.01D0,0.03D0,0.05D0,0.1D0,
     &     0.2D0,0.3D0,0.4D0,0.5D0,0.6D0,0.7D0,0.8D0,0.9D0,0.95D0,
     &     0.99D0,1.0D0/
C=======================================================================
        
C       calculate specific energy of particles
C       --------------------------------------

        DO 100 i=ifirst,ntot
        etoti(i-ifirst+1) = 0.5D0 * (xdot(1,i)**2 + xdot(2,i)**2 +
     &                         xdot(3,i)**2) - phidbl(i)
100     CONTINUE

C      calculate number of bound particles
C      -----------------------------------

        nbound = 0
        DO 150 i=1,ntot-ifirst+1
           IF(etoti(i).LT.0.D0) nbound = nbound + 1
150     CONTINUE

C       sort for particle energy
C       ------------------------

        nnpart = ntot-ifirst+1
        CALL indexx(nnpart,etoti,index)

C       initialize tensor of inertia
C       ----------------------------

        DO 210 i=1,3
           DO 200 k=1,3
              ti(i,k) = 0.D0
200        CONTINUE
210     CONTINUE

C       LOOP over fraction of most bound particles and all particles
C       ------------------------------------------------------------

        nstart   = 1
*
        DO 500 ief=1,nef1
*
        xmass(ief) = 0.d0
        npart(ief) = 0
        r3av(ief) = 0.d0
        r2av(ief) = 0.d0
        zav(ief) = 0.d0
        vrav(ief) = 0.d0
        vrot(ief) = 0.d0
        vzav(ief) = 0.d0

           IF(ief.LE.nef) THEN
C                                  only fraction of bound particles
C                                  --------------------------------
              nnext = NINT(xf(ief) * nbound)
           ELSE
C                                   all particles
C                                   -------------
              nnext = ntot - ifirst + 1
           ENDIF

C-----------------------------------------------------------------
C--      at least two particles are required for ellipticity...
C-----------------------------------------------------------------
           IF(nnext.LT.2) THEN
              ba(ief)  = 999.
              ca(ief)  = 999.
              taue(ief) = 999.
              DO 320 k=1,3
                 DO 310 j=1,3
                    evec(k,j,ief) = 0.
310              CONTINUE
320           CONTINUE

           ELSE

C       calculate tensor of inertia relative to density centre (R.Sp.)
C       ---------------------------
              DO 400 i=nstart,nnext
                 ipo = index(i) + ifirst - 1
                 xpo = x(1,ipo)-rdens(1)
                 ypo = x(2,ipo)-rdens(2)
                 zpo = x(3,ipo)-rdens(3)
*
*       compute data for analogue to Lagrangian radii (RSp Oct 2022)
*
                 xmass(ief) = xmass(ief) + body(ipo)
                 npart(ief) = npart(ief) + 1
*
                 rrad2 = xpo**2 + ypo**2
                 rrad = rrad2 + zpo**2
                 rrad2 = dsqrt(rrad2)
                 rrad = dsqrt(rrad)
                 r2av(ief) = r2av(ief) + body(ipo)*rrad2
                 r3av(ief) = r3av(ief) + body(ipo)*rrad
                 zav(ief) = zav(ief) + body(ipo)*zpo
*
                 vri = 0.d0
                 do k = 1,3
                 vri = vri + xdot(k,ipo)*(x(k,ipo)-rdens(k))/rrad
                 end do
*     X-Y plane projected position square
                 rr12 = x(1,ipo)**2 + x(2,ipo)**2
*     X-Y plane radial velocity value * sqrt(RR12)
                 xr12 = 0.d0
                 do k = 1,2
                     xr12 = xr12 + xdot(k,ipo)*(x(k,ipo)-rdens(k))
                 end do
                 vrad(ipo) = xr12/dsqrt(rr12)
                 vrav(ief) = vrav(ief) + body(ipo)*vrad(ipo)
*
*     Rotational velocity
                 vrot1 = xdot(1,ipo) - xr12/rr12*x(1,ipo)
                 vrot2 = xdot(2,ipo) - xr12/rr12*x(2,ipo)
*     Rotational direction sign
                 xsign = vrot1*x(2,ipo)/dsqrt(rr12) -
     &                   vrot2*x(1,ipo)/dsqrt(rr12)
                 vrotm = dsqrt(vrot1**2+vrot2**2)
*     Cumulate average rotational velocity in shell, mass weighted
                 if(xsign.gt.0.d0) then
                     vtan(ipo) = vrotm
                     vrot(ief) = vrot(ief) + body(ipo)*vrotm
                 else
                     vtan(ipo) = -vrotm
                     vrot(ief) = vrot(ief) - body(ipo)*vrotm
                 end if
*
                 vzav(ief) = vzav(ief) + body(ipo)*xdot(3,ipo)
*
                 ti(1,1) = ti(1,1) + body(ipo) * (ypo**2 + zpo**2)
                 ti(2,2) = ti(2,2) + body(ipo) * (xpo**2 + zpo**2)
                 ti(3,3) = ti(3,3) + body(ipo) * (xpo**2 + ypo**2)
                 ti(1,2) = ti(1,2) - body(ipo) * xpo*ypo
                 ti(1,3) = ti(1,3) - body(ipo) * xpo*zpo
                 ti(2,3) = ti(2,3) - body(ipo) * ypo*zpo
400           CONTINUE
*
              r3av(ief) = r3av(ief)/xmass(ief)
              r2av(ief) = r2av(ief)/xmass(ief)
              zav(ief) = zav(ief)/xmass(ief)
              vrav(ief) = vrav(ief)/xmass(ief)
              vrot(ief) = vrot(ief)/xmass(ief)
              vzav(ief) = vzav(ief)/xmass(ief)
              avmass(ief) = xmass(ief)/dble(npart(ief))
*
              sigr2(ief) = 0.d0
              sigph2(ief) = 0.d0
              sigz2(ief) = 0.d0
              DO 430 i=nstart,nnext
                 ipo = index(i) + ifirst - 1
                 sigr2(ief) = sigr2(ief) + 
     &                body(ipo)*(vrad(ipo)-vrav(ief))**2
                 sigph2(ief) = sigph2(ief) +
     &                body(ipo)*(vtan(ipo)-vrot(ief))**2
                 sigz2(ief) = sigz2(ief) + 
     &                        body(ipo)*(xdot(3,ipo)-vzav(ief))**2
430           CONTINUE

              sigr2(ief) = sigr2(ief)/xmass(ief)
              sigph2(ief) = sigph2(ief)/xmass(ief)
              sigz2(ief) = sigz2(ief)/xmass(ief)
*
C       set off-axis values by symmetry
C       -------------------------------

              ti(2,1) = ti(1,2)
              ti(3,1) = ti(1,3)
              ti(3,2) = ti(2,3)
*
C=======================================================
C       determine eigenvalues and axis of inertia
C=======================================================

C------------------------------------------------------
C--          copy tensor of inertia
C------------------------------------------------------
              DO 420 i=1,3
                 DO 410 k=1,3
                    tiwork(i,k) = ti(i,k) 
410              CONTINUE
420           CONTINUE
              np = 3

C------------------------------------------------------
C--          calculate eigenvalues and eigenvectors
C------------------------------------------------------
              CALL tred2(tiwork,np,np,dwork,ework)
              CALL tqli(dwork,ework,np,np,tiwork)

C--               sort for increasing eigenvalues
              CALL indexx(np,dwork,indexev)
C--               find eigenvectors
              DO 450 i=1,np
                 lam(i) = dwork(indexev(i))
                 DO 440 k=1,np
                    evec(k,i,ief) = tiwork(k,indexev(i))
440              CONTINUE
450           CONTINUE

              xhelp    = lam(3) + lam(2) - lam(1)
              xhelp1   = lam(2) - lam(3) + lam(1)
c              IF(xhelp1.LT.0.D0) THEN
c                 PRINT*,' ellan: xhelp1 < 0',xhelp1,tnow
c                 xhelp1 = 0.D0
c             ENDIF
              ba(ief)  = SQRT(MAX(tiny2,lam(3)-lam(2)+lam(1)) / xhelp)
              ca(ief)  = SQRT(MAX(tiny2,xhelp1) / xhelp) 
              taue(ief) = (ba(ief)-ca(ief)) /MAX(tiny2,(1.D0 - ca(ief)))

              nstart = nnext + 1

           ENDIF

500     CONTINUE

            xmcum(1) = xmass(1)
            ncum(1) = npart(1)
        DO 600 ief=2,nef1
            xmcum(ief) = xmcum(ief-1) + xmass(ief)
            ncum(ief) = ncum(ief-1) + npart(ief)
600     CONTINUE
C==================================================================
C==         OUTPUT of data
C===================================================================
      if(rank.eq.0)then
            WRITE (6,40) (XF(K),K=1,NEF1)
 40         FORMAT (/,11X,'TIME   E/ET:',1P,18(1X,E9.2))
            WRITE (6,401) TTOT,(XMASS(K),K=1,NEF1)
 401        FORMAT (3X,1P,E12.4,' MSHELL:',18(1X,E9.2))
            WRITE (6,402) TTOT,(XMCUM(K),K=1,NEF1)
 402        FORMAT (3X,1P,E12.4,'   MCUM:',18(1X,E9.2))
            WRITE (6,403) TTOT,(NPART(K),K=1,NEF1)
 403        FORMAT (3X,1P,E12.4,'  NPART:',18(1X,I9))
            WRITE (6,404) TTOT,(NCUM(K),K=1,NEF1)
 404        FORMAT (3X,1P,E12.4,'   NCUM:',18(1X,I9))
            WRITE (6,405) TTOT,(AVMASS(K),K=1,NEF1)
 405        FORMAT (3X,1P,E12.4,' AVMASS:',18(1X,E9.2))
            WRITE (6,406) TTOT,(R3AV(K),K=1,NEF1)
 406        FORMAT (3X,1P,E12.4,'   R3AV:',18(1X,E9.2))
            WRITE (6,407) TTOT,(R2AV(K),K=1,NEF1)
 407        FORMAT (3X,1P,E12.4,'   R2AV:',18(1X,E9.2))
            WRITE (6,408) TTOT,(ZAV(K),K=1,NEF1)
 408        FORMAT (3X,1P,E12.4,'    ZAV:',18(1X,E9.2))
            WRITE (6,412) TTOT,(VROT(K),K=1,NEF1)
 412        FORMAT (3X,1P,E12.4,'   VROT:',18(1X,E9.2))
            WRITE (6,413) TTOT,(VRAV(K),K=1,NEF1)
 413        FORMAT (3X,1P,E12.4,'   VRAV:',18(1X,E9.2))
            WRITE (6,414) TTOT,(VZAV(K),K=1,NEF1)
 414        FORMAT (3X,1P,E12.4,'   VZAV:',18(1X,E9.2))
            WRITE (6,415) TTOT,(SIGR2(K),K=1,NEF1)
 415        FORMAT (3X,1P,E12.4,'  SIGR2:',18(1X,E9.2))
            WRITE (6,416) TTOT,(SIGPH2(K),K=1,NEF1)
 416        FORMAT (3X,1P,E12.4,' SIGPH2:',18(1X,E9.2))
            WRITE (6,417) TTOT,(SIGZ2(K),K=1,NEF1)
 417        FORMAT (3X,1P,E12.4,'  SIGZ2:',18(1X,E9.2))
*
            WRITE (6,41) TTOT,(BA(K),K=1,NEF1)
 41         FORMAT (3X,1P,E12.4,'   B/A: ',18(1X,E9.2))
            WRITE (6,42) TTOT,(CA(K),K=1,NEF1)
 42         FORMAT (3X,1P,E12.4,'   C/A: ',18(1X,E9.2))
            WRITE (6,43) TTOT,(TAUE(K),K=1,NEF1)
 43         FORMAT (3X,1P,E12.4,'   TAU: ',18(1X,E9.2))
      END IF
*
        RETURN
        
        END
