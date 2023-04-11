      SUBROUTINE MODIFY
*
*
*       Parameter modification at restart.
*       ----------------------------------
*
      INCLUDE 'common6.h'
      EXTERNAL VERIFY
      character*1 level
*
       NAMELIST /ININPUT/ N, NFIX, NCRIT, NRAND, NNBOPT, NRUN, NCOMM,
     &    ETAI, ETAR, RS0, DTADJ, DELTAT, TCRIT,QE, RBAR, ZMBAR,
     &    KZ, DTMIN, RMIN, ETAU, ECLOSE, GMIN, GMAX, SMAX, LEVEL
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
*
      RBAR0 = RBAR
      ZMBAR0 = ZMBAR
*
      if(rank.eq.0)then
          READ (NML=ININPUT, IOSTAT=IIC, UNIT=5)
          IF (IIC.NE.0) WRITE (6, '("INPUT: invalid Namelist format")')
*
          IF(KSTART.GT.2)THEN
              PRINT*,' KSTART > 2 obsolete in new namelist input, STOP'
#if MPIINIT
              CALL MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
              STOP
#endif
          END IF
*       Keep tidal circularization together with GR (RS March 2019)
          KZ273 = 0
          IF(KZ(27).GT.2) THEN
             KZ273 = 3
             KZ(27) = KZ(27) - 2
          END IF
*
      end if
*
*     Allow step for HDF5 files to be larger than one. (R.Sp. Sep. 2021)
      HDFMAX = SMAX
*     Check SMAX to make sure it have correct value
      if(rank.eq.0) then
         IF(SMAX.GT.1) THEN
            print*, 'Warning! SMAX > 1.0, reduce to 1.0.'
            SMAX = 1.0
         ELSE
            DTNSMAX = 1.0
 1          IF(SMAX/DTNSMAX.LE.0.75) THEN
               DTNSMAX = 0.5D0*DTNSMAX
               IF(DTNSMAX.GT.1E-19) GO TO 1
            END IF
            SMAX = DTNSMAX
         END IF
      end if
*
      TCRIT = TTOT + TCRIT
*
      if(rank.eq.0)then
         WRITE (6,10)
 10      FORMAT(//,30X,'RESTART PARAMETERS:',/,30X,'==================')
         WRITE (6,11) TTOT/TCR0, TIME/TCR0, TCRIT/TCR0, TTOT, TIME,
     &            TCRIT
 11      FORMAT (/,' START AT TTOT/TIME ',2F16.8,' STOP INTENDED AT ',
     &            F16.8,' TCR0',/,' START AT TTOT/TIME ',2F16.8,
     &            ' STOP INTENDED AT ',F16.8,' NBODY-UNITS ',/)
*
         WRITE (6,12)  N, NFIX, NCRIT, NRAND, NNBOPT, NRUN, NCOMM
   12    FORMAT (/,I16,I6,2I7,I8,2I6)
*       Note that RBAR and ZMBAR should be taken from original data read by mydump.
         WRITE(6,15)
   15    FORMAT (/,12X,' ETAI      ETAR      RS0       DTADJ',
     &                  '     DELTAT',
     &                  '     TCRITP    TCRIT     QE',
     &                  '        RBAR0      ZMBAR0')
         WRITE (6,20)  ETAI, ETAR, RS0, DTADJ, DELTAT, TCRITP, TCRIT,
     &              QE, RBAR0, ZMBAR0
   20    FORMAT (/,10X,1P10E10.1)
*
         WRITE (6,22)
   22    FORMAT (/,12X,'OPTIONS')
         WRITE (6,24)  (J,J=1,50)
   24    FORMAT (/,9X,50I3)
         WRITE (6,26)  (KZ(J),J=1,50)
   26    FORMAT (/,9X,50I3)
         WRITE (6,28)
   28    FORMAT (/,12X,'DTMIN     RMIN      ETAU      ECLOSE    GMIN',
     &        '      GMAX     SMAX')
         WRITE (6,30)  DTMIN, RMIN, ETAU, ECLOSE, GMIN, GMAX, SMAX
   30    FORMAT (/,9X,1P7E10.1)
         WRITE (6,31)
   31    FORMAT(/,30X,'=====================================',/)
*
         RBAR = RBAR0
         ZMBAR = ZMBAR0
      end if
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
#ifdef TT
*** FlorentR - restart with new tensors
      IF (KZ(14).EQ.9) THEN
        CALL TTINIT
      ENDIF
*** FRenaud
#endif
*        Read stellar evolution and collision parameters.
      CALL READSE(LEVEL)
*
*       Perform a simple validation check on main input parameters.
      CALL VERIFY
*
*       Save the new parameters on tape/disc in case a restart is needed.
      CALL MYDUMP(1,1)
*
      RETURN
*
      END
