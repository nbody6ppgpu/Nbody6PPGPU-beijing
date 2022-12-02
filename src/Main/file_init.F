      SUBROUTINE FILE_INIT
*
*       Opening of Files with proper names.
*       ------------------------------------
*
      INCLUDE 'common6.h'
*
      CHARACTER*10 FILE(100)
      CHARACTER*7 FILESTAT
      INTEGER STAT
      LOGICAL LSKIP
*
#ifdef PARALLEL
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
*
      DO IFILE = 1,100
      FILE(IFILE) = '          '
      END DO
*
*      Warning following files not used, opened in other routines, some with time flags (R.Sp. Aug. 2021)
*      1 and 2 in mydump.F for full restarts
      FILE(1)  = 'comm.1    '
      FILE(2)  = 'comm.2_...'
*      3,9,19,40,50,82,83,87 opened with time flags for multiple files
      FILE(3)  = 'conf.3_...'
      FILE(9)  = 'bdat.9_...'
      FILE(19) = 'bwdat.19_.'
      FILE(40) = 'single.40 '  ! also binary.40 and merger.40
      FILE(50) = 'tail.50_..'
      FILE(82) = 'bev.82_...'
      FILE(83) = 'sev.83_...'
      FILE(87) = 'hidat.87_.'
*
      FILE(40) = 'single.40 '  ! also binary.40 and merger.40
      FILE(42) = 'ibcoll.42 '
      FILE(52) = 'itid3.52  '
*
      FILE(4)  = 'degen.4   '
      FILE(7)  = 'lagr.7    '
      FILE(8)  = 'bdat.8    '
      FILE(10) = 'dat.10    '
      FILE(11) = 'esc.11    '
      FILE(12) = 'hiarch.12 '
      FILE(13) = 'coll.13   '
      FILE(14) = 'shrink.14 '
      FILE(15) = 'mix.15    '
      FILE(16) = 'hirect.16 '
      FILE(17) = 'binev.17  '
      FILE(18) = 'pbin.18   '
      FILE(20) = 'symb.20   '
      FILE(21) = 'datsev.21 '
      FILE(22) = 'rocdeg.22 '
      FILE(24) = 'coal.24   '
      FILE(25) = 'sediag.25 '
      FILE(26) = 'cmbody.26 '
      FILE(29) = 'highv.29  '
      FILE(30) = 'global.30 '
      FILE(31) = 'escbin.31 '
      FILE(33) = 'ns.33     '
      FILE(34) = 'bh.34     '
      FILE(35) = 'event.35  '
      FILE(36) = 'status.36 '
      FILE(38) = 'sediag.38 '
      FILE(39) = 'hbin.39   '
      FILE(41) = 'nbflow.41 '
      FILE(43) = 'sediag.43 '
      FILE(44) = 'hinc.44   '
      FILE(45) = 'mbh.45    '
      FILE(46) = 'mbhnb.46  '
      FILE(54) = 'close.54  '
      FILE(55) = 'rpmax.55  '
      FILE(71) = 'cirdiag.71'
      FILE(73) = 'histab.73 '
      FILE(75) = 'cirdiag.75'
      FILE(77) = 'kscrit.77'
      FILE(81) = 'chstab.81 '
      FILE(84) = 'merger.84 '
      FILE(85) = 'roche.85  '
      FILE(86) = 'cmbody.86 '
      FILE(88) = 'dtmin.88  '
      FILE(89) = 'quastab.89'
      FILE(91) = 'bs.91     '
      FILE(95) = 'wdcirc.95 '
      FILE(96) = 'cirdiag.96'
#ifdef PARALLEL
      if(rank.eq.0) then
#endif
* At beginning of run make sure all files are deleted, except those to
* read from or stdin, stdout (5,6,10,21) .
      IF(TTOT.EQ.0.D0)THEN
         DO IFILE = 1,100
         LSKIP = IFILE.EQ.5.OR.IFILE.EQ.6.OR.IFILE.EQ.10.OR.IFILE.EQ.21
         IF(LSKIP) GOTO 123
      OPEN(UNIT=IFILE,IOSTAT=STAT,FILE=FILE(IFILE),STATUS='OLD',ERR=123)
         IF(STAT.EQ.0) CLOSE(IFILE,STATUS='delete')
 123     CONTINUE
         END DO
*
         FILESTAT = 'NEW    '
      ELSE
         FILESTAT = 'UNKNOWN'
      END IF
*
      
      IF (KZ(9).GE.3)
     &OPEN (UNIT=4,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(4),
     &             ACCESS='APPEND')
      IF (KZ(7).GE.2)
     &OPEN (UNIT=7,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(7),
     &             ACCESS='APPEND')
      IF (KZ(8).GT.0.OR.NBIN0.GT.0)
     &OPEN (UNIT=8,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(8),
     &             ACCESS='APPEND')
      IF (KZ(22).GT.0)
     &OPEN (UNIT=10,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(10))
      IF (KZ(23).EQ.2.OR.KZ(23).EQ.4.OR.KZ(23).EQ.3)
     &OPEN (UNIT=11,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(11),
     &             ACCESS='APPEND')
      IF (KZ(18).EQ.1.OR.KZ(18).EQ.3)
     &OPEN (UNIT=12,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(12),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=13,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(13),
     &             ACCESS='APPEND')

      OPEN (UNIT=14,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(14),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=15,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(15),
     &             ACCESS='APPEND')
      IF (KZ(27).EQ.2.or.KZ(34).GT.0.or.KZ(28).GT.0)
     &OPEN (UNIT=16,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(16),
     &             ACCESS='APPEND')
      IF (KZ(9).GE.2)
     &OPEN (UNIT=17,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(17),
     &             ACCESS='APPEND')
      IF (KZ(8).GT.0)
     &OPEN (UNIT=18,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(18),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=20,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(20),
     &             ACCESS='APPEND')
      IF(KZ(19).EQ.2.OR.KZ(19).EQ.4)
     &OPEN (UNIT=21,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FILE(21))
      IF (KZ(34).GT.0)
     &OPEN (UNIT=22,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(22),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3) 
     &OPEN (UNIT=24,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(24),
     &             ACCESS='APPEND')
      IF (KZ(27).GT.0) 
     &OPEN (UNIT=25,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(25),
     &             ACCESS='APPEND')

      OPEN (UNIT=26,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(26),
     &             ACCESS='APPEND')

      OPEN (UNIT=29,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(29),
     &             ACCESS='APPEND')

      OPEN (UNIT=30,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(30),
     &             ACCESS='APPEND')
      IF (KZ(23).EQ.2.OR.KZ(23).EQ.4)
     &OPEN (UNIT=31,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(31),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=33,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(33),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=34,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(34),
     &             ACCESS='APPEND')
      IF (KZ(19).GT.0.OR.KZ(27).GT.0) 
     &OPEN (UNIT=35,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(35),
     &             ACCESS='APPEND')

      OPEN (UNIT=36,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(36),
     &     ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=38,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(38),
     &             ACCESS='APPEND')
      IF (KZ(9).EQ.1.or.KZ(9).EQ.3)
     &OPEN (UNIT=39,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(39),
     &             ACCESS='APPEND')

      OPEN (UNIT=41,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(41),
     &             ACCESS='APPEND')
      IF (KZ(8).EQ.1.or.KZ(8).GE.3)
     &OPEN (UNIT=42,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(42),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=43,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(43),
     &             ACCESS='APPEND')
      IF (KZ(27).GT.0)
     &OPEN (UNIT=44,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(44),
     &             ACCESS='APPEND')
      IF (KZ(24).EQ.1)
     &OPEN (UNIT=45,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(45),
     &             ACCESS='APPEND')
      IF (KZ(24).EQ.1)
     &OPEN (UNIT=46,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(46),
     &             ACCESS='APPEND')
      IF (KZ(14).EQ.3)
     &OPEN (UNIT=52,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(52),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=54,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(54),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=55,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(55),
     &             ACCESS='APPEND')
      IF (KZ(27).GT.0)
     &OPEN (UNIT=71,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(71),
     &             ACCESS='APPEND')
      IF (KZ(15).GT.0)
     &OPEN (UNIT=73,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(73),
     &             ACCESS='APPEND')
      IF (KZ(27).EQ.2)
     &OPEN (UNIT=75,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(75),
     &             ACCESS='APPEND')
      IF (KZ(16).GT.2.AND.(KZ(30).GT.0.OR.KZ(30).EQ.-1))
     &OPEN (UNIT=77,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(77),
     &             ACCESS='APPEND')
      IF (KZ(15).GT.0.AND.(KZ(30).GT.0.or.KZ(30).EQ.-1))
     &OPEN (UNIT=81,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(81),
     &             ACCESS='APPEND')
      IF (KZ(8).GE.2)
     &OPEN (UNIT=84,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(84),
     &             ACCESS='APPEND')
      IF (KZ(34).GT.0)
     &OPEN (UNIT=85,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(85),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3)
     &OPEN (UNIT=86,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(86),
     &             ACCESS='APPEND')
      OPEN (UNIT=88,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(88),
     &             ACCESS='APPEND')
      IF (KZ(15).GE.3)
     &OPEN (UNIT=89,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(89),
     &             ACCESS='APPEND')
      IF (KZ(19).GE.3) 
     &OPEN (UNIT=91,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(91),
     &             ACCESS='APPEND')
      IF (KZ(27).GT.0) 
     &OPEN (UNIT=95,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(95),
     &             ACCESS='APPEND')
      IF (KZ(27).GT.0) 
     &OPEN (UNIT=96,STATUS=FILESTAT,FORM='FORMATTED',FILE=FILE(96),
     &             ACCESS='APPEND')
#ifdef PARALLEL
      end if
#endif
*
*     DO jfile=0,9
*     WRITE(6,556)(rank,10*jfile+ifile,file(10*jfile+ifile),ifile=1,10)
*556  FORMAT(1X,' rank,ifile,file=',10(2I4,A11))
*     END DO
*
      RETURN
      END
