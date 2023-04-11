      SUBROUTINE INPUT
*
*
*       Parameter input.
*       ----------------
*
      INCLUDE 'common6.h'
      EXTERNAL VERIFY
      character*1 level
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
*
*
       NAMELIST /ININPUT/ N, NFIX, NCRIT, NRAND, NNBOPT, NRUN, NCOMM,
     &    ETAI, ETAR, RS0, DTADJ, DELTAT, TCRIT,QE, RBAR, ZMBAR,
     &    KZ, DTMIN, RMIN, ETAU, ECLOSE, GMIN, GMAX, SMAX, LEVEL
*      Define Default Values
*       
*
      IF(rank.eq.0)THEN
      READ (NML=ININPUT, IOSTAT=IIC, UNIT=5)
      IF (IIC.NE.0) WRITE (6, '("INPUT: invalid Namelist format")')
*
*       Read & print the main input parameters.
*        READ (5,*)  N, NFIX, NCRIT, NRAND, NNBOPT, NRUN, NCOMM
C Termination time in physical units, TCRITP, read in nbody6.F
*        READ (5,*)  ETAI, ETAR, RS0, DTADJ, DELTAT, TCRIT,
*    &               QE, RBAR, ZMBAR
*        READ (5,*)  (KZ(J),J=1,50)
*        READ (5,*)  DTMIN, RMIN, ETAU, ECLOSE, GMIN, GMAX, SMAX
*       Keep tidal circularization together with GR (RS March 2019)
      KZ273 = 0
      IF(KZ(27).GT.2) THEN
         KZ273 = 3
         KZ(27) = KZ(27) - 2
      END IF
*
      END IF
*     Allow step for HDF5 files to be larger than one. (R.Sp. Sep. 2021)
      HDFMAX = SMAX
*     Check SMAX to make sure it have correct value
      if(rank.eq.0) THEN
         IF(SMAX.GT.1) THEN
            print*, 'Warning! SMAX > 1.0, reduce to 1.0.'
            SMAX = 1.0
         else
            DTNSMAX = 1.0
 1          IF(SMAX/DTNSMAX.LE.0.75) THEN
               DTNSMAX = 0.5D0*DTNSMAX
               IF(DTNSMAX.GT.1E-19) GO TO 1
            END IF
            SMAX = DTNSMAX
         END IF
      END IF
*
#if MPIINIT
      CALL MPI_BCAST(N,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NFIX,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NCRIT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NRAND,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NNBOPT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NRUN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(NCOMM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
*
      CALL MPI_BCAST(KZ(1),50,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(KZ273,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
*
      CALL MPI_BCAST(ETAI,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ETAR,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RS0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DTADJ,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DELTAT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TCRIT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(QE,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RBAR,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ZMBAR,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DTMIN,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RMIN,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      RMIN2 = RMIN**2 
      CALL MPI_BCAST(ETAU,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ECLOSE,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(GMIN,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(GMAX,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(SMAX,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(LEVEL,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif
*
      if(rank.eq.0)then
         WRITE (6,10)
   10    FORMAT (/////,15X,'N  NFIX  NCRIT  NRAND  NNBOPT  NRUN NCOMM')
         WRITE (6,12)  N, NFIX, NCRIT, NRAND, NNBOPT, NRUN, NCOMM
   12    FORMAT (/,I16,I6,2I7,I8,2I6)
*
C New: (Aug.1998, P.Kroupa)
         WRITE(6,15)
   15    FORMAT (//,12X,' ETAI      ETAR      RS0       DTADJ',
     &                  '     DELTAT',
     &                  '     TCRITP    TCRIT     QE', 
     &                  '        RBAR       ZMBAR')
         WRITE (6,20)  ETAI, ETAR, RS0, DTADJ, DELTAT, TCRITP, TCRIT, 
     &              QE, RBAR,
     &              ZMBAR
   20    FORMAT (/,10X,1P10E10.1)
*
         WRITE (6,22)
   22    FORMAT (//,12X,'OPTIONS')
         WRITE (6,24)  (J,J=1,50)
   24    FORMAT (/,9X,50I3)
         WRITE (6,26)  (KZ(J),J=1,50)
   26    FORMAT (/,9X,50I3)
         WRITE (6,28)
   28    FORMAT (//,12X,'DTMIN     RMIN      ETAU      ECLOSE    GMIN',
     &        '      GMAX     SMAX')
         WRITE (6,30)  DTMIN, RMIN, ETAU, ECLOSE, GMIN, GMAX, SMAX
   30    FORMAT (/,9X,1P7E10.1,/)
      end if
      call flush(6)
*        Read stellar evolution and collision parameters.
      CALL READSE(LEVEL)
*
*       Define total particle number & neighbour membership range.
      NTOT = N
      NZERO = N
      NNBMAX = MIN(N/2,LMAX - 50)
      ZNBMIN = MAX(0.01*FLOAT(NNBMAX),1.0)
      ZNBMAX = 0.9*FLOAT(NNBMAX)
*       Save initial ETAI.
      ETA0 = ETAI
      RSMIN = RS0
      RC = RS0
*
*       Perform a simple validation check on main input parameters.
      CALL VERIFY
*
      GPRINT(1) = 0.0
      DELTAS = 0.0
*     Suppress this KZ(4) since the output need special analysis tool
C      IF (KZ(4).GT.0) THEN
C*       Read parameters for binary evolution analysis.
C          K = KZ(4)
C          if(rank.eq.0)then
C          READ (5,*)  DELTAS, ORBITS(1), (GPRINT(J),J=1,K)
C          end if
C*
C#if MPIINIT
C      CALL MPI_BCAST(DELTAS,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
C      CALL MPI_BCAST(ORBITS(1),9,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
C      CALL MPI_BCAST(GPRINT(1),9,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
C#endif
C*
C      if(rank.eq.0)WRITE (6,40)  DELTAS, ORBITS(1), (GPRINT(J),J=1,K)
C   40     FORMAT (//,12X,'DELTAS =',F6.2,'  ORBITS(1) =',F6.2,
C     &                                            '  GPRINT(J) =',9F7.3)
C*       Modify binary output factor by perturbation at different levels.
C          DO 50 L = 2,K
C              ORBITS(L) = ORBITS(1)*(GPRINT(1)/GPRINT(L))**0.3333
C   50     CONTINUE
C      END IF
*
C Old version:
*       Set random number skip for routine DATA.
c      IDUM1 = NRAND
C NEW version (14.08.98, P.Kroupa):
C*       Set random number SEED for routine DATA.
      IDUM1 = -1*NRAND
c+++ Notify others of this change on log file:
C      if(rank.eq.0)then
C      write(6,*)
C      write(6,*)' ****** NOTE: new random number seed initialisation!'
C      write(6,*)' ****** AND new ran2 from new ed. of Press et al.'
C      write(6,*)
C      end if
*
*
*       Save square of c.m. approximation parameter (coupled to GMIN).
      CMSEP2 = GMIN**(-0.666666667)
*
      RETURN
*
      END
*
      SUBROUTINE READSE(LEVEL)
*
*       Parameter input for SSE/BSE/Coll (R.Sp. Mar23)
*       ----------------------------------------------
*
      INCLUDE 'common6.h'
* Input commons for SSE/BSE===RSp.Mar23======
* sse parameters for hrdiag, mlwind
      integer ecflag,wdflag,nsflag,psflag,mdflag,bhflag,kmech,idum
      real*8 mch,mxns0,mxns1,nwind,bwind,flbv,disp,ecsig,
     *       wdsig1,wdsig2,wdkmax,vvfac
      common/sse/ecflag,wdflag,nsflag,psflag,mdflag,bhflag,
     *       kmech,idum,mch,mxns0,mxns1,nwind,bwind,flbv,disp,ecsig,
     *       wdsig1,wdsig2,wdkmax,vvfac
      namelist/insse/ecflag,wdflag,nsflag,psflag,mdflag,bhflag,
     *       kmech,mch,mxns0,mxns1,nwind,bwind,flbv,disp,ecsig,
     *       wdsig1,wdsig2,wdkmax,vvfac
* bse parameters for comenv, kick, roche, mdot, coal, mix
      integer bhspinfl,kicktype
      real*8 lambd1,alphac,xk2,xk3,acc1,acc2,xbeta,xxi,
     *      epsnov,eddfac,gamm1
      common/bse/lambd1,alphac,bhspinfl,kicktype,xk2,xk3,
     *      acc1,acc2,xbeta,xxi,epsnov,eddfac,gamm1
      namelist/inbse/lambd1,alphac,bhspinfl,kicktype,xk2,xk3,
     *      acc1,acc2,xbeta,xxi,epsnov,eddfac,gamm1
      character*1 level
* collision parameters for coal, mix
      real*8 fctorcl
      common/coll/fctorcl
      namelist/incoll/fctorcl
* End Input for SSE/BSE=====================

*       Start SSE/BSE/Coll Parameter Definition
      IF(LEVEL.EQ.'A')THEN
*
          ecflag = 0
          wdflag = 1
          nsflag = 1
          psflag = 0
          mdflag = 1
          bhflag = 0
*
          nwind = 0.5d0   !  nwind = neta  for mlwind
          bwind = 0.d0
          flbv = 1.5d0
          kmech = 1
          disp = 190.d0
          ecsig = 20.d0
          wdsig1 = 2.d0
          wdsig2 = 2.d0
          wdkmax = 6.d0
          vvfac = 0.d0     ! vvfac = vfac for kick
*
          lambdac = 0.5d0  ! lambdac = lambda for comenv
          alphac = 3.d0    ! alphac = alpha1  for comenv
          bhspinfl = 0
          kicktype = -2
*
          xbeta = 0.125d0  ! xbeta = beta in mdot, roche
          xxi = 1.d0       ! xxi = xi in mdot, roche
          epsnov = 1.d-3
          eddfac = 1.d2
          gamm1 = -1.d0
*
      ELSE IF(LEVEL.EQ.'B')THEN
          ecflag = 1
          wdflag = 1
          nsflag = 3
          psflag = 0
          mdflag = 3
          bhflag = 2
*
          nwind = 0.5d0   !  nwind = neta  for mlwind
          bwind = 0.d0
          flbv = 1.5d0
          kmech = 1
          disp = 265.d0
          ecsig = 3.d0
          wdsig1 = 2.d0
          wdsig2 = 2.d0
          wdkmax = 6.d0
          vvfac = 0.d0     ! vvfac = vfac for kick
*
          lambdac = 0.5d0  ! lambdac = lambda for comenv
          alphac = 3.d0    ! alphac = alpha1  for comenv
          bhspinfl = 0
          kicktype = -2
*
          xbeta = 0.125d0  ! xbeta = beta in mdot, roche
          xxi = 1.d0       ! xxi = xi in mdot, roche
          epsnov = 1.d-3
          eddfac = 1.d2
          gamm1 = -1.d0
*
      ELSE IF(LEVEL.EQ.'C')THEN
          ecflag = 1
          wdflag = 1
          nsflag = 4
          psflag = 1
          mdflag = 4
          bhflag = 2
*
          nwind = 0.477d0   !  nwind = neta  for mlwind
          bwind = 0.d0
          flbv = 1.5d0
          kmech = 4
          disp = 265.d0
          ecsig = 3.d0
          wdsig1 = 2.d0
          wdsig2 = 2.d0
          wdkmax = 6.d0
          vvfac = 0.d0     ! vvfac = vfac for kick
*
          lambdac = 0.0d0  ! lambdac = lambda for comenv
          alphac = 1.d0    ! alphac = alpha1  for comenv
          bhspinfl = 0
          kicktype = -2
*
          xbeta = 0.125d0  ! xbeta = beta in mdot, roche
          xxi = 1.d0       ! xxi = xi in mdot, roche
          epsnov = 1.d-3
          eddfac = 1.d2
          gamm1 = -1.d0
*
*      ELSE IF(LEVEL.EQ.'D')THEN
      ELSE IF(LEVEL.EQ.'0')THEN
          PRINT*,' Run without stellar evolution'
      ELSE
          PRINT*,' No SSE/BSE Level defined - STOP '
#if MPIINIT
      CALL MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
      STOP
#endif
      END IF
*       Following Parameters defined independently of Level.
      IF(LEVEL.NE.'0') THEN        
          mch = 1.44
          mxns0 = 1.8d0
          mxns1 = 2.5d0
*
          xk2 = 0.21d0     !  xk2 = k2,k21,k22 argument hrdiag
          xk3 = 0.21d0     !  xk3 = k3 in mdot, roche
          acc1 = 3.920659d8
          acc2 = 1.5d0
*
          fctorcl = 0.5d0
      END IF
*
          if(rank.eq.0.and.level.ne.'0')then
              read(nml=insse,iostat=iic,unit=5)
              if (iic.ne.0)print*,' insse: invalid Namelist format'
              read(nml=inbse,iostat=iic,unit=5)
              if (iic.ne.0)print*,' inbse: invalid Namelist format'
              read(nml=incoll,iostat=iic,unit=5)
              if (iic.ne.0)print*,' incoll: invalid Namelist format'
*
              write(6,"(/,' SSE/BSE/Coll Level = ',A1,/)")level

      write(6,30)ecflag,wdflag,nsflag,psflag,mdflag,bhflag,
     *       kmech,mch,mxns0,mxns1,nwind,bwind,flbv,disp,ecsig,
     *       wdsig1,wdsig2,wdkmax,vvfac
 30   format(10X,' SSE FLAGS: EC,WD,NS,PS,MD,BH=',6I4,' KMECH=',I4,/,
     * 10X,' MASSES: Chandra, NS0, NS1=',1P,3D11.2,/,
     * 10X,' WIND PARAMS: NETA, BWIND, FLBV=',
     * 3D11.2,/,10X,' KICK PARAMS: DISP, ECSIG, WDSIG1,2,',
     * '  WDKMAX, VFAC=',6D11.2,/)
      write(6,40)lambd1,alphac,bhspinfl,kicktype,xk2,xk3,
     *      acc1,acc2,xbeta,xxi,epsnov,eddfac,gamm1
 40   format(10X,' BSE PARAMS: LAMBDA, ALPHA=',1P,2D11.2,' BHSPINFL,',
     * ' KICKTYPE=',2I4,/,10X,' K2,K3=',2D11.2,' ACC1,2, BETA, XI=',
     * 4D11.2,/,10X,' EPSNOV, EDDFAC, GAMM1=',3D11.2,/)
      write(6,50)fctorcl
 50   FORMAT(10X,' COLL PARAMS: FactorCl=',1P,D11.2,/)
*             write(nml=insse,unit=6)
*             write(nml=incoll,unit=6)
          end if
*
*       End SSE/BSE/Coll Parameter Definition
      RETURN
*
      END
