*       pulsar_ctrl.h
*       -------------
*       Pulsar/neutron-star birth-distribution and CE-accretion control
*       parameters. Set once via the INPULSAR namelist in READPULSAR
*       (input.F) and read (never written) by the physics routines in
*       pulsar.F. This is the single declaration of COMMON/PULSAR/;
*       every routine that touches these parameters must INCLUDE this
*       file rather than re-declaring the COMMON block, so that layout
*       and types cannot drift out of sync between callers.
*
*       PSR_ACC_CE, PSR_PMODE, PSR_BMODE, PSRM_PMODE, PSRM_BMODE are
*       mode selectors and must be integer (see manual for the field
*       list and allowed values); everything else is REAL*8.
*
      REAL*8  PSR_SPINA,PSR_SPINB,PSR_BMAGA,PSR_BMAGB,
     &        PSRM_SPINA,PSRM_SPINB,PSRM_BMAGA,PSRM_BMAGB,
     &        PSR_BMIN,PSR_RAD,PSR_TAU,PSR_KAPPA,PSR_EP
      INTEGER PSR_PMODE,PSR_BMODE,PSRM_PMODE,PSRM_BMODE,PSR_ACC_CE
*
      COMMON/PULSAR/ PSR_SPINA,PSR_SPINB,PSR_BMAGA,PSR_BMAGB,
     &               PSRM_SPINA,PSRM_SPINB,PSRM_BMAGA,PSRM_BMAGB,
     &               PSR_BMIN,PSR_RAD,PSR_TAU,PSR_KAPPA,PSR_EP,
     &               PSR_PMODE,PSR_BMODE,PSRM_PMODE,PSRM_BMODE,
     &               PSR_ACC_CE
