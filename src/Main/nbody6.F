      PROGRAM NBODY6
*
*             N B O D Y 6++
*             *************
*
*       Regularized AC N-body code with triple & binary collisions.
*       --------------------------------------------------------
*
*       Hermite integration scheme with block-steps (V 4.0.0 April/99).
*       ------------------------------------------------------------------
*
*       Developed by Sverre Aarseth, IOA, Cambridge.
*       ............................................
*       Message Passing Version NBODY6++ for Massively Parallel Systems
*       Developed by Rainer Spurzem, ARI, Heidelberg
*       
*       Hybrid parallelization (GPU, AVX/SSE, OpenMP + MPI) 
*       Developed by Long Wang, KIAA, Peking University
*
      INCLUDE 'common6.h'
      INCLUDE 'timing.h'
      include 'omp_lib.h'
      COMMON/STSTAT/  TINIT,NIR,NIB,NRGL,NKS
#ifdef DEBUG
*     --10/03/14 19:40-lwang-debug--------------------------------------*
***** Note:------------------------------------------------------------**
      COMMON/adt/ adtime,dumptime,dprintt,dtprint,namep
*     --10/03/14 19:41-lwang-end----------------------------------------*
#endif
      EXTERNAL MERGE
      NAMELIST /INNBODY6/ KSTART, TCOMP, TCRTP0,
     &    isernb,iserreg,iserks,adtime,dumptime,dprintt,dtprint,namep
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#define MPIINIT 0
#endif

#if MPIINIT
*       Initialize MPI
      CALL MPI_INIT(ierr)
*     CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,ithread,ierr)
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,group,ierr)
      CALL MPI_GROUP_SIZE(group,isize,ierr)
      CALL MPI_GROUP_RANK(group,rank,ierr)
      ibarcount=0
      write(6,11) rank,isize,group
 11   format('MPI-initial: This is rank=',I6,' size=',I6,' group=',I11)
#endif
*
*       Initialize the timer.
      CALL CPUTIM(ttota)

*       Get threads number
#ifdef OMP
!$omp parallel 
      icore=OMP_get_num_threads()
!$omp end parallel
      PRINT*,'RANK: ',rank,' OpenMP Number of Threads: ',icore
#else
      icore = 1
#endif

#ifdef PARALLEL
      call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif      
      call flush(6)
*
*       Read start/restart indicator & CPU time.
*     19:41-lwang-debug-10/03/14-----------------------------*
      IF(rank.eq.0)THEN
#ifdef DEBUG
      READ (NML=INNBODY6, IOSTAT=IIC, UNIT=5)
      IF (IIC.NE.0) WRITE (6, '("NBODY6D: invalid Namelist format")')
*     READ (5,*)  KSTART, TCOMP, TCRTP0,
*    &    isernb,iserreg,iserks,adtime,dumptime,dprintt,dtprint,namep
#else
      READ (NML=INNBODY6, IOSTAT=IIC, UNIT=5)
      IF (IIC.NE.0) WRITE (6, '("NBODY6: invalid Namelist format")')
*     READ (5,*)  KSTART, TCOMP, TCRTP0,
*    &    isernb,iserreg,iserks
#endif
      END IF
#if MPIINIT
      CALL MPI_BCAST(isernb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iserreg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iserks,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(KSTART,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TCOMP,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(TCRTP0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
*
#ifdef DEBUG
      CALL MPI_BCAST(adtime,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dumptime,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dprintt,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dtprint,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(namep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
*
      isernb = max(isize,isernb*icore)
      iserreg = max(isize,iserreg*icore)
*      iserks = max(isize,iserks*icore)

      IF(rank.eq.0) THEN
        PRINT*,' MPI: iserreg,isernb,iserks,icore=',
     &        iserreg,isernb,iserks,icore
#ifdef DEBUG
        PRINT*,' DEBUG: adtime,dumptime,dprintt,dtprint,namep =',
     &      adtime,dumptime,dprintt,dtprint,namep
#endif
      END IF
#endif
*
      IF (KSTART.EQ.1) THEN
*
          TCRITP = TCRTP0
          if(rank.eq.0)WRITE (6,*) ' KSTART, TCOMP, TCRITP =',
     &                               KSTART, TCOMP, TCRITP
*       Read input parameters, perform initial setup and obtain output.
          CPU = TCOMP
          CALL START
          call cputim(tt7)
          CALL ADJUST
          call cputim(tt8)
          if(rank.eq.0)ttadj = ttadj + (tt8-tt7)*60.
      ELSE
*
*       Read previously saved COMMON variables from tape/disc on unit 1.
*       Backup kstart value before call mydump
          KSTART0 = KSTART
          CALL MYDUMP(0,1)
*       Overwrite old TCRITP from input
          IF(TCRTP0.GT.0.D0) TCRITP = TCRTP0
          if(rank.eq.0)WRITE (6,*) ' KSTART, TCOMP, TCRITP =',
     &                               KSTART, TCOMP, TCRITP
*       Reset kstart to input value
          KSTART = KSTART0
*       
          IF (NDUMP.GE.3) STOP
*       Safety indicator preventing repeated restarts set in routine CHECK.
          CPU = TCOMP
          CPU0 = 0.0 
*       Set IPHASE = -1 for new NLIST in routine INTGRT (Hermite version).
          IPHASE = -1
*
*       Initialize evolution parameters which depend on metallicity.
          IF (KZ(19).GE.3) THEN
              CALL ZCNSTS(ZMET,ZPARS)
          END IF
*
*       Check reading modified restart parameters (KSTART = 3, 4 or 5).
*       Or control output RSp March 2019 (KSTART = 2)
          CALL MODIFY
*
*       Open all other files.
          if(rank.eq.0)CALL FILE_INIT
*
*       Find massive back hole after restart
          IF (KZ(24).EQ.1) call IMBHRESTART

      END IF
*
* (R.Sp.)Set time flag and step number flags for beginning of run
      TINIT = TTOT
      NIR = NSTEPI
      NIB = NSTEPB
      NRGL = NSTEPR
      NKS = NSTEPU
*
      call cputim(tt2)
      if(rank.eq.0)ttinitial = ttinitial + (tt2-ttota)*60.
*       Advance solutions until next output or change of procedure.
    1 CONTINUE
      call cputim(tt1)
*
*     if(rank.eq.0.and.time.lt.0.01)print*,' Bef INTGRT t=',time
      CALL INTGRT
*     if(rank.eq.0.and.time.lt.0.01)print*,' Aft INTGRT t=',time
*
      call cputim(tt2)
      if(rank.eq.0)ttint = ttint + (tt2-tt1)*60.
*
      IF (IPHASE.EQ.1) THEN
*       Prepare new KS regularization.
      call cputim(tt1)
          CALL KSREG
          CALL FLUSH(6)
      call cputim(tt2)
      if(rank.eq.0)ttks = ttks + (tt2-tt1)*60.
      if(rank.eq.0)ttksinit = ttksinit + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.2) THEN
*       Terminate KS regularization.
      call cputim(tt1)
          CALL KSTERM
      call cputim(tt2)
      if(rank.eq.0)ttks = ttks + (tt2-tt1)*60.
      if(rank.eq.0)ttksterm = ttksterm + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.3) THEN
*       Perform energy check & parameter adjustments and print diagnostics.
          call cputim(tt7)
          CALL ADJUST
          call cputim(tt8)
          if(rank.eq.0)ttadj = ttadj + (tt8-tt7)*60.
*
      ELSE IF (IPHASE.EQ.4) THEN
*       Switch to unperturbed three-body regularization.
      call cputim(tt1)
          ISUB = 0 
          CALL TRIPLE(ISUB)
      call cputim(tt2)
      if(rank.eq.0)ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.5) THEN
*       Switch to unperturbed four-body regularization.
      call cputim(tt1)
          ISUB = 0
          CALL QUAD(ISUB)
      call cputim(tt2)
      if(rank.eq.0)ttks = ttks + (tt2-tt1)*60.
*
*       Adopt c.m. approximation for inner binary in hierarchical triple.
      ELSE IF (IPHASE.EQ.6) THEN
      call cputim(tt1)
          CALL MERGE
      call cputim(tt2)
      if(rank.eq.0)ttks = ttks + (tt2-tt1)*60.
*
      ELSE IF (IPHASE.EQ.7) THEN
*       Restore old binary in hierarchical configuration.
      call cputim(tt1)
          CALL RESET
      call cputim(tt2)
      if(rank.eq.0)ttks = ttks + (tt2-tt1)*60.
*
*       Begin chain regularization.
      ELSE IF (IPHASE.EQ.8) THEN
      call cputim(tt1)
          ISUB = 0
          TIMENB = TTOT
          CALL CHAIN(ISUB,TIMENB)
      call cputim(tt2)
      if(rank.eq.0)ttks = ttks + (tt2-tt1)*60.
      END IF
*
*       Continue integration.
      GO TO 1
*
      END
