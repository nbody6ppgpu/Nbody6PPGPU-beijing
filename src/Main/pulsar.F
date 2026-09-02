      SUBROUTINE PULSAREVO(mass, mass_comp, rad_comp, B_i,
     & period_i, pdot_i, psrmdot, dtm, nsincl, ce, psrn,
     & B_f, period_f, pdot_f)
* Calculates updated pulsar parameters (spin, spin down rate,
* and magnetic field) as described in Chattopadhyay+2020.

* Input:
* mass [Msun] - Mass of the pulsar.
* mass_comp [Msun] - Mass of the companion (non-degenerate).
* rad_comp [Rsun] - Radius of the companion (non-degenerate).
* B_i [Tesla] - Initial magnetic field strength for the time step
* period_i [second]- Initial spin period of the pulsar.
* pdot_i [second/second]- Initial spin down rate (period derivative) of the pulsar.
* psrmdot [Msun/second] - Mass transfer rate
* dtm [second] - Time-step.
* nsincl [radian] - The angle between pulsar's spin axis and magnetic field axis.
* ce -  Logical variable indicating whether ce is occuring.
* psrn - name of pulsar (NAMENS entry); identifies the object in the
*   warnings this routine can still emit (RLOF cap reached, bad RLOF
*   step, unhandled input).
*
* Outputs:
* B_f [Tesla] - Final magnetic field strength for the time step
* period_f [second]- Final spin period of the pulsar.
* pdot_f [second/second]- Initial spin down rate (period derivative) of the pulsar.

* NBODY6++GPU header files. Comment out if not using NBODY6++GPU.
      INCLUDE 'common6.h'
* MOCCA header files. Comment out if not using MOCCA.
* Please add appropriate header files here for MOCCA:


      integer :: psrn ! NAMENS entry; used only to identify objects in warnings.
      REAL(kind=8) :: mass,mass_comp,rad_comp,psrmdot,dtm,nsincl

      LOGICAL :: ce ! if common envelope happens or not.

      REAL(kind=8) :: PULSAR_TAU, PULSAR_KAPPA, PULSAR_BMIN, PULSAR_EP

      !counters for splitting timestep in spin-down calculations.
      integer :: counter, psrloop
      REAL(8) :: SDTSCALE !pulsar spin-down timescale, P/Pdot

      ! Lower limit of mass transfer rate to consider spin-up.
      REAL(kind=8), PARAMETER :: acc = 1e-20
      ! Constants
      REAL(kind=8), PARAMETER :: pi = 3.14159265358979323D0
      REAL(kind=8), PARAMETER :: mu_0 = 4.0D-7 * pi ! H/m Permeability of free space.
      REAL(kind=8), PARAMETER :: c = 2.99792458D8  ! m / s
      REAL(kind=8), PARAMETER :: G = 6.67430D-11  ! N * m^2 * kg^-2
      REAL(kind=8), PARAMETER :: G_TO_T = 1.0d-4   !
      REAL(kind=8), PARAMETER :: MSOL_TO_KG = 1.989D30 !
      REAL(kind=8), PARAMETER :: KM_TO_M = 1.0D3
      REAL(kind=8), PARAMETER :: myr_in_sec = 3.1536d13
      REAL(kind=8), PARAMETER :: RSOL_TO_M = 6.957d8

      ! variables used to calculate RLOF spin-up in smaller steps.
      INTEGER :: divide_timestep_by
      REAL(8) :: PERIOD_COPY
      REAL(8) :: PDOT_COPY
      REAL(8) :: BMAG_COPY
      REAL(8) :: MASS_COPY
      LOGICAL :: RLOF_LOOP_DONE
      INTEGER :: COUNT, K
      REAL(8) :: new_dtm, mass_f

      REAL(kind=8) :: mdot_Msun_sec, consAlf, p, q, alfRad, vKep,
     &                omegaKep, vDiff, J, Jacc, MoIf

      ! common envolope mass transfer parameters:
      real(kind=8) :: cea, ceb, dm, psrmdot_eff
      REAL(kind=8), PARAMETER :: MAX_CE_DM = 0.1D0

      INTEGER :: MAX_RLOF_DIVIDE
      LOGICAL :: BAD_RLOF_STEP,RLOF_CAP_REACHED

      ! Gather pulsar inputs (single declaration in pulsar_ctrl.h).
      INCLUDE 'pulsar_ctrl.h'
      ! PSR_BMIN (Gauss): LOG10 of minimum magnetic field of neutron star at the end of magnetic field decay
      ! PSR_RAD (km) : Radius of neutron star, in the future it will be transferred to an option of NS EOS to choose NS radius.
      ! PSR_TAU (Myr): NS Magnetic decay time scale.
      ! PSR_KAPPA (Msun): NS Magnetic decay mass scale.
      ! PSR_EP : Efficiency of transfer of angular momentum in RLOF.
      ! PSR_ACC_CE : Determining the pulsar accretion scenario in common envelope evolution.
      !               0: 'ZERO': No accretion onto NS during CE.
      !               1: 'DISK': Forming accretion disk around CE similar to RLOF
      !               2: 'SURFACE': Accreted materials fall directly onto surface of NS.

      ! Converting the inputs from input file into SI units for calculations.
      PULSAR_TAU = PSR_TAU * myr_in_sec
      PULSAR_BMIN = 10.0d0**PSR_BMIN * G_TO_T
      PULSAR_KAPPA = PSR_KAPPA * MSOL_TO_KG
      PULSAR_EP = PSR_EP

      MAX_RLOF_DIVIDE = 100000
      RLOF_CAP_REACHED = .FALSE.

      B_f = B_i
      period_f = period_i
      pdot_f = pdot_i
      IF (dtm.LE.0.0D0) THEN
*        A non-positive timestep would make ISOSPINDOWN's
*        EXP(-step/PULSAR_TAU) exceed 1 (B growing instead of
*        decaying) and can drive Pfsquared negative, so SQRT would
*        return NaN; the RLOF/CE branches divide by dtm and are no
*        better off. Should not happen in normal forward-time
*        integration; guard rather than silently propagate Inf/NaN.
         WRITE(6,*)'WARNING: PULSAREVO called with non-positive dtm,',
     &        ' pulsar unchanged this step:', psrn, dtm
         RETURN
      ENDIF
      ! Update isolated (no mass transfer case) pulsar spin, spin down rate, and magnetic field
      ! Calculate surface magnetic field due to magnetic braking, n=3
       IF (((psrmdot.LT.acc).or.(psrmdot.gt.acc
     &  .and.ce.and.PSR_ACC_CE.eq.0)).and.pdot_i.gt.0)THEN
            SDTSCALE = period_i/pdot_i ! neutron star spin-down time scale.
            psrloop = 1
            ! If the timestep is larger than the spin-down timescale,
            ! then the spin-down calculation needs to be broken into smaller steps
            ! to avoid potential numerical overshooting problems.
            do while ((dtm/psrloop).gt.SDTSCALE)
               psrloop = psrloop * 5
            end do
            B_f = B_i
            period_f = period_i
            pdot_f = pdot_i
            TSSIZE = dtm/psrloop
            do counter = 1, psrloop
               period_i = period_f
               B_i = B_f
               pdot_i = pdot_f
               CALL ISOSPINDOWN(psrn,B_i, period_i,pdot_i,
     &                B_f, period_f, pdot_f, TSSIZE)
            end do

      ! Update stable RLOF(mass transfer, not CE)  pulsar spin, spin down rate, and magnetic field
      ! Calculate surface magnetic field due to mass accretion and magnetic burying
        ELSEIF (((psrmdot.gt.acc.and.ce.eqv..false.) .or.
     &    (psrmdot.gt.acc.and.ce.and.PSR_ACC_CE.eq.1))
     &    .and.pdot_i.gt.0) THEN
         !As implemented in Chattopadhyay+2020
            !Save copies of initial condition for the second round of calculations.
            PERIOD_COPY = period_i
            PDOT_COPY = pdot_i
            BMAG_COPY = B_i
            MASS_COPY = mass
            divide_timestep_by = 100
            period_f = period_i
            pdot_f= pdot_i
            B_f = B_i
            mass_f = mass
            RLOF_LOOP_DONE = .false.
            psrmdot_eff = psrmdot
            IF (ce.AND.psrmdot_eff*dtm.GT.MAX_CE_DM) THEN
               psrmdot_eff = MAX_CE_DM/dtm
            ENDIF
            DO WHILE (RLOF_LOOP_DONE.eqv..false.)
               COUNT = 0
               divide_timestep_by = divide_timestep_by * 2
               IF (divide_timestep_by.GT.MAX_RLOF_DIVIDE) THEN
                  RLOF_CAP_REACHED = .TRUE.
                  WRITE(6,*)'WARNING: PULSAR RLOF cap reached ',
     &                 psrn, divide_timestep_by
                  B_f = BMAG_COPY
                  period_f = PERIOD_COPY
                  pdot_f = PDOT_COPY
                  mass_f = MASS_COPY
                  EXIT
               ENDIF

               new_dtm = dtm / divide_timestep_by
               BAD_RLOF_STEP = .FALSE.

               period_f = PERIOD_COPY
               pdot_f = PDOT_COPY
               B_f = BMAG_COPY
               mass_f = MASS_COPY

               DO K = 1, divide_timestep_by

                  period_i = period_f
                  pdot_i = pdot_f
                  B_i = B_f
                  mass = mass_f

                  CALL NSACCRETE(psrmdot_eff, new_dtm, mass,
     &                 PSR_RAD, period_i, B_i, pdot_i,
     &                 period_f, B_f, pdot_f, mass_f)

                  IF (period_f.LE.0.0D0 .OR. pdot_f.LE.0.0D0
     &               .OR. B_f.LE.0.0D0 .OR. mass_f.LE.0.0D0) THEN
                     BAD_RLOF_STEP = .TRUE.
                     EXIT
                  ENDIF

               END DO

               IF (.NOT.BAD_RLOF_STEP) THEN
                  RLOF_LOOP_DONE = .true.
               ENDIF

            END DO

            new_dtm = dtm / divide_timestep_by
            period_i = PERIOD_COPY
            pdot_i = PDOT_COPY
            B_i = BMAG_COPY
            mass = MASS_COPY

            IF (RLOF_CAP_REACHED.eqv..FALSE.) THEN

               DO K = 1, divide_timestep_by
                  CALL NSACCRETE(psrmdot_eff, new_dtm, mass, PSR_RAD,
     &            period_i,B_i,pdot_i,period_f,B_f,pdot_f,mass_f)
                  IF (period_f.LE.0.0D0 .OR. pdot_f.LE.0.0D0
     &               .OR. B_f.LE.0.0D0 .OR. mass_f.LE.0.0D0) THEN
                  WRITE(6,*)'WARNING: bad final RLOF pulsar step ',
     &               psrn, period_f, pdot_f, B_f, mass_f
                     B_f = B_i
                     period_f = period_i
                     pdot_f = pdot_i
                     mass_f = mass
                     EXIT
                  ENDIF
                  period_i=period_f
                  B_i=B_f
                  pdot_i=pdot_f
                  mass = mass_f
               END DO
            ELSE
               WRITE(6,*)'WARNING: PULSAR RLOF ',
     &          'cap reached, no spin-up applied'
            ENDIF

      ! Update CE pulsar spin, spin down rate, and magnetic field
      ! Calculate surface magnetic field due to mass accretion and magnetic burying
      ELSEIF ((psrmdot.gt.acc.and.ce.and.PSR_ACC_CE.eq.2)
     &   .and.pdot_i.gt.0) THEN
         !Fitting formuale from MacLeod&Ramirez-Ruiz 2015 as implemented in Chattopadhyay+2020
         cea=(-1.1*1e-5)*mass_comp+(1.2*1e-4)
         ceb=(1.5*1e-2)*mass_comp+(-1.5*1e-1)
         dm=cea*rad_comp+ceb
         dm=MIN(dm,MAX_CE_DM)
        !Update magnetic field
        B_f=PULSAR_BMIN+EXP(-dm*MSOL_TO_KG/PULSAR_KAPPA)
     &   *(B_i-PULSAR_BMIN)
        ! WRITE(*,*) dm, B_i, B_f
        MoI = NSMOI(mass, PSR_RAD)
        ! Calculate angular momentum increase due to accretion
        Jacc = MoI*dm*SQRT(G*mass*MSOL_TO_KG/(PSR_RAD*KM_TO_M))/mass
        J = (2.0d0 * pi / period_i) * MoI   ! GW
        J = J + Jacc
        period_f = 2.0d0 * pi * MoI / J
        mass = mass + dm
        MoIf = NSMOI(mass, PSR_RAD)
        pdot_f = NSSPINDOWN(period_f, B_f, PSR_RAD, MoIf)
      ELSE
*        None of the three evolution branches above matched. With
*        PSR_ACC_CE forced to 0 at input and all PSR_*MODE=0 birth
*        options removed, pdot_i should be strictly positive for the
*        life of the object (see manual/source-map.md); reaching here
*        means that invariant broke. Return the inputs unchanged
*        (already done by the B_f/period_f/pdot_f = *_i assignments
*        above) rather than continue on an unvalidated state.
         WRITE(6,*)'WARNING: PULSAREVO invariant violated (pdot',
     &        ' non-positive, or psrmdot exactly at the acc',
     &        ' threshold); pulsar unchanged this step:', psrn,
     &        psrmdot, pdot_i, ce
      ENDIF

      CONTAINS

      SUBROUTINE ISOSPINDOWN(psrn, B_i, period_i,pdot_i,
     &                     B_f, period_f, pdot_f,step)
      !This subroutine calculates the updated magnetic field,
      !spin, and spin-down rate for an ISOlated pulsar experiencing
      !normal spindown.

         REAL(kind=8) :: atrm, up, MoI, dwn, Pfsquared
         REAL(kind=8) :: omega, cons, trm1, trm2, trm3, pDotTop, pDotBottom

         REAL(kind=8) :: period_i,period_f,B_i,B_f,pdot_i,pdot_f,step
         INTEGER(4) :: psrn
         B_f = PULSAR_BMIN + EXP(-step / PULSAR_TAU)
     &      * (B_i - PULSAR_BMIN)
          MoI = NSMOI(mass, PSR_RAD)
          atrm = SIN(nsincl)**2.0d0
          up = 8.0d0 * pi * pi * ((PSR_RAD*KM_TO_M)**6.0d0)
     &       * atrm

          dwn = 3.0d0 * mu_0 * (c**3.0d0) * MoI

          cons = 2.0d0 * up / dwn
          trm1 = (PULSAR_BMIN**2.0d0) * step
          trm2 = PULSAR_TAU * PULSAR_BMIN * (B_f - B_i)
          trm3 = (PULSAR_TAU / 2.0d0) * (B_f**2.0d0 - B_i**2.0d0)
          Pfsquared = cons*(trm1 - trm2 - trm3) + period_i * period_i
          period_f = SQRT(Pfsquared)
          pdot_f = NSSPINDOWN(period_f, B_f, PSR_RAD, MoI)
      END SUBROUTINE ISOSPINDOWN

      SUBROUTINE NSACCRETE(psrmdot, dtm, mass, PSR_RAD,
     &    period_i, B_i, pdot_i, period_f, B_f, pdot_f, mass_f)
      ! This function updates the spin, spindown rate and
      ! magnetic field of the neutron star when the binary
      ! undergoes RLOF.
         REAL(8) :: psrmdot, dtm, mass, PSR_RAD, period_i,
     &     B_i, pdot_i, MoI, MoIf, psrmdot_kg_sec,
     &     consAlf, period_f, B_f, pdot_f, mass_f
         ! Converting the mass tranfser rate to SI units.
         psrmdot_kg_sec = psrmdot*MSOL_TO_KG
         ! Calculate Alfven radius.
         consAlf = (2.0d0*pi*pi/G/mu_0/mu_0)**(1.0d0/7.0d0)
         p = ((PSR_RAD*KM_TO_M)**6.0d0) /
     &     ((mass*MSOL_TO_KG)**0.5d0 * psrmdot_kg_sec)
         q = p**(2.0d0 / 7.0d0)
         alfRad = consAlf * q * B_i**(4.0d0 / 7.0d0)

         ! Update magnetic field
         B_f = PULSAR_BMIN + EXP(-(psrmdot_kg_sec * dtm)
     &   / PULSAR_KAPPA) * (B_i - PULSAR_BMIN)
         ! Calculate keplerian velocity
         vKep = SQRT(2.0d0*G*mass*MSOL_TO_KG / alfRad)
         omegaKep = 2.0d0 * vKep / alfRad
         ! Calculate difference of velocity
         vDiff = omegaKep - (2*pi/period_i)
         MoI = NSMOI(mass, PSR_RAD)
         ! Calculate increase of angular momentum due to accretion
         Jacc = 0.25d0 * PULSAR_EP * (alfRad**2.0d0)
     &      * (psrmdot_kg_sec) * vDiff * dtm
         J = (2*pi/period_i) * MoI
         ! New angular momentum
         J = J + Jacc
         ! Calculate updated spin period
         period_f = 2.0d0 * pi * MoI / J
         ! Mass after accretion.
         mass_f = mass + psrmdot * dtm
         MoIf = NSMOI(mass_f,PSR_RAD)
         ! Updated spin down rate
         pdot_f = NSSPINDOWN(period_f, B_f, PSR_RAD, MoIf)
      END SUBROUTINE NSACCRETE

      FUNCTION NSMOI(mass, PSR_RAD)
      ! Calculate the moment of inertia of the neutron star given:
      ! mass[Msun]: mass of NS
      ! PSR_RAD[km]: radius of NS.
         IMPLICIT NONE
         REAL(kind=8), INTENT(IN) :: mass, PSR_RAD
         REAL(kind=8) :: M_R, NSMOI
         M_R = mass / PSR_RAD

         NSMOI = 0.237d0 * (mass*MSOL_TO_KG) * ((PSR_RAD*KM_TO_M)**2)

         NSMOI = NSMOI * (1.0d0 + (4.2d0 * M_R)
     &    + 90.0d0 * (M_R**4))

      END FUNCTION NSMOI

      FUNCTION NSSPINDOWN(P, B, R, MOI)
      ! Calculate the spindown rate of the neutron star given:
      ! P [sec]: spin period
      ! B [Tesla]: magnetic field
      ! R [km]: radius
      ! MOI [kg m^2]: moment of inertia
         REAL(kind=8), INTENT(IN) :: P, B, R, MOI
         REAL(kind=8) :: NSSPINDOWN
         REAL(kind=8) :: radius_6, bfield_2
         REAL(kind=8) :: const1, const2, pDTop, pDBot

         radius_6 = (R*KM_TO_M) ** 6
         bfield_2 = B**2
         const1 = 8.0 * pi * pi
         const2  = 3.0d0 * mu_0 * c * c * c
         pDTop = const1 * radius_6 * bfield_2
         pDBot = const2 * MOI * P
         NSSPINDOWN = pDTop / pDBot
      END FUNCTION NSSPINDOWN

      END SUBROUTINE PULSAREVO

* ====================================================================================
      SUBROUTINE PULSARINIT(PULSARINITB, PULSARINITP,
     &     PULSARINITPD, PSRMASS, ISAIC)
* Initialize pulsar parameters based on the selections
* and distributions chosen in the setup.
* Input:
* PSRMASS [Msol] - mass of the pulsar
* ISAIC [bool] - if the pulsar is formed through AIC as a PSRM.
*
* Output:
* PULSARINITB [Tesla] - Magnetic field strength at pulsar birth
* PULSARINITP [second] - Spin period of the pulsar at birth
* PULSARINITPD [second/second] - Spin down rate of the pulsar at birth
*

* NBODY6++GPU header file(s), random number generator, and inputs.
* Comment out if not using NBODY6++GPU.
      INCLUDE 'common6.h'
      REAL*8  :: RAN2
      INCLUDE 'pulsar_ctrl.h'
      ! PSR_RAD : Radius of neutron star, in the future it will be transferred to an option of NS EOS to choose NS radius.
      ! PSR_TAU : NS Magnetic decay time scale.
      ! PSR_KAPPA : NS Magnetic decay mass scale.
      REAL(kind=8) :: MOI, radius_6, bfield_2, pDTop, pDBot
      REAL(kind=8) :: const1, const2
      REAL(8) :: PULSARINITB, PULSARINITP, PULSARINITPD, PSRMASS
      REAL(8) :: log10P, log10B, u, diff, bmaga10,bmagb10
      LOGICAL :: ISAIC
* constants
      REAL(kind=8), PARAMETER :: pi = 3.14159265358979323846264338D0
      REAL(kind=8), PARAMETER :: mu_0 = 4.0D-7 * pi ! H/m Permeability of free space.
      REAL(kind=8), PARAMETER :: c = 2.99792458D8  ! m / s
      REAL(kind=8), PARAMETER :: G = 6.67430D-11  !N m^2 kg-2
      REAL(kind=8), PARAMETER :: G_TO_T = 1.0d-4   !Gauss to Tesla conversion.
      REAL(kind=8), PARAMETER :: MSOL_TO_KG = 1.989D30 !
      REAL(kind=8), PARAMETER :: KM_TO_M = 1.0D3
      REAL(kind=8), PARAMETER :: RSOL_TO_KM = 6.957D5


      ! PSR_SPINA & PSR_SPINB / PSRM_SPINA & PSRM_SPINB:  first and second definitions for pulsar birth spin period. See below.
      ! PSR_PMODE / PSRM_PMODE: Distribution of pulsar spin period at birth.
      !        1 = FIXED distribution, all pulsar born with spin period at PSR_SPINA sec.
      !        2 = UNIFORM distribution, PSR_SPINA is the lower limit and PSR_SPINB is the upper limit of the distribution.
      !        3 = NORMAL distribution, PSR_SPINA is the mean and PSR_SPINB is the stdev of the distribution.
      !        4 = LOGUNIFORM / Flat-in-log distribution, PSR_SPINA is the lower limit and PSR_SPINB is the upper limit of the distribution.
      ! PSR_BMAGA & PSR_BMAGB / PSRM_BMAGA & PSRM_BMAGB:  first and second definitions for pulsar birth magfield. See below.
      ! PSR_BMODE & PSRM_BMODE: Distribution of pulsar magnetic field strength at birth.
      !        1 = FIXED distribution, all pulsar born with magfield at 10^PSR_BMAGA Gauss.
      !        2 = UNIFORM distribution, 10^PSR_BMAGA is the lower limit and 10^PSR_BMAGB is the upper limit of the distribution.
      !        4 = LOGUNIFORM / Flat-in-log distribution, 10^PSR_BMAGA is the lower limit and 10^PSR_BMAGB is the upper limit of the distribution.
      !        5 = LOGNORMAL distribution, PSR_BMAGA is the LOG10(mean) and PSR_BMAGB is the LOG10(stdev) of the distribution.
      !

* Please add MOCCA header files or anything I missed below
* Comment out if not using MOCCA.
      ! EXTERNAL ran3
      ! REAL(4) :: ran3

      ! converting MOCCA input variables to NBODY variables.
      ! REAL(8) PSR_SPINA,PSR_SPINB,PSR_BMAGA,PSR_BMAGB,
      !    PSR_BMIN,PSR_RAD,PSR_TAU,PSR_KAPPA,PSR_EP,PSR_PMODE,PSR_BMODE

      ! INTEGER :: PULSAR_INITIAL_SPIN_MODEL
      ! COMMON /C_PULSAR_INITIAL_SPIN_MODEL/ PULSAR_INITIAL_SPIN_MODEL
      ! REAL(8) :: PULSAR_INITIAL_SPIN_FIXED, PULSAR_INITIAL_SPIN_MIN, PULSAR_INITIAL_SPIN_MAX
      ! REAL(8) :: PULSAR_INITIAL_SPIN_MEAN, PULSAR_INITIAL_SPIN_STD
      ! COMMON /C_PULSAR_INITIAL_SPIN/ &
      !    PULSAR_INITIAL_SPIN_FIXED, PULSAR_INITIAL_SPIN_MIN, PULSAR_INITIAL_SPIN_MAX, &
      !    PULSAR_INITIAL_SPIN_MEAN, PULSAR_INITIAL_SPIN_STD

      ! INTEGER :: PULSAR_INITIAL_B_MODEL
      ! COMMON /C_PULSAR_INITIAL_B_MODEL/ PULSAR_INITIAL_B_MODEL
      ! REAL(8) :: PULSAR_INITIAL_B_FIXED, PULSAR_INITIAL_B_MIN, PULSAR_INITIAL_B_MAX
      ! REAL(8) :: PULSAR_INITIAL_B_MEAN, PULSAR_INITIAL_B_STD
      ! COMMON /C_PULSAR_INITIAL_B/ &
      !    PULSAR_INITIAL_B_FIXED, PULSAR_INITIAL_B_MIN, PULSAR_INITIAL_B_MAX, &
      !    PULSAR_INITIAL_B_MEAN, PULSAR_INITIAL_B_STD

      ! PSR_PMODE = PULSAR_INITIAL_SPIN_MODEL
      ! IF (PSR_PMODE.EQ.1) THEN
      !   PSR_SPINA = PULSAR_INITIAL_SPIN_FIXED
      ! ELSEIF (PSR_PMODE.EQ.2) THEN
      !   PSR_SPINA = PULSAR_INITIAL_SPIN_MIN
      !   PSR_SPINB = PULSAR_INITIAL_SPIN_MAX
      ! ELSEIF (PSR_PMODE.EQ.3) THEN
      !   PSR_SPINA = PULSAR_INITIAL_SPIN_MEAN
      !   PSR_SPINB = PULSAR_INITIAL_SPIN_STD
      ! ELSEIF (PSR_PMODE.EQ.4) THEN
      !   PSR_SPINA = PULSAR_INITIAL_SPIN_MIN
      !   PSR_SPINB = PULSAR_INITIAL_SPIN_MAX

      ! PSR_BMODE = PULSAR_INITIAL_B_MODEL
      ! IF (PSR_BMODE.EQ.1) THEN
      !   PSR_BMAGA = PULSAR_INITIAL_B_FIXED
      ! ELSEIF (PSR_BMODE.EQ.2) THEN
      !   PSR_BMAGA = PULSAR_INITIAL_B_MIN
      !   PSR_BMAGB = PULSAR_INITIAL_B_MAX
      ! ELSEIF (PSR_BMODE.EQ.3) THEN
      !   PSR_BMAGA = PULSAR_INITIAL_B_MEAN
      !   PSR_BMAGB = PULSAR_INITIAL_B_STD
      ! ELSEIF (PSR_BMODE.EQ.4) THEN
      !   PSR_BMAGA = PULSAR_INITIAL_B_MIN
      !   PSR_BMAGB = PULSAR_INITIAL_B_MAX
      ! ELSEIF (PSR_BMODE.EQ.5) THEN
      !   PSR_BMAGA = PULSAR_INITIAL_B_MEAN
      !   PSR_BMAGB = PULSAR_INITIAL_B_STD


* Setting up random number seed in NBODY6++GPU
* Comment out if not using NBODY6++GPU
      KDUM = IDUM1
      if (ISAIC.eqv..false.) then

      ! Check distribution type and compute pulsar spin at birth accordingly
      IF(PSR_PMODE.EQ.1) THEN
         ! Set NS initial spin to a
         ! fixed constant value as indicated by PSR_SPINA
         PULSARINITP = PSR_SPINA

      ELSEIF(PSR_PMODE.EQ.2) THEN
         ! Set NS initial spin to Uniform distribution,
         ! minimum as PSR_SPINA, and maximum as PSR_SPINB
         diff = PSR_SPINB - PSR_SPINA
         ! using random number generator in NBODY6++GPU
         PULSARINITP = PSR_SPINA + RAN2(KDUM) * diff
         ! using random number generator in MOCCA
         ! PULSARINITP = PSR_SPINA + ran3() * diff

      ELSEIF (PSR_PMODE.EQ.3) THEN
         ! Set to Normal distribution from  Faucher-Giuere&Kaspi 2006 (https://iopscience.iop.org/article/10.1086/501516)
         ! mean as PSR_SPINA, stdev as PSR_SPINB
         PULSARINITP = -1.0d0
         do while (PULSARINITP < 0.0d0)
            PULSARINITP = GAUS_RAND(PSR_SPINA, PSR_SPINB)
         end do
         !write(6,*)"PulsarINITP ",PULSARINITP,PSR_SPINA,u,PSR_RAD

      ELSEIF (PSR_PMODE.EQ.4) THEN
         ! Set to Flat-in-log distribution
         ! minimum as PSR_SPINA, maximum as PSR_SPINB
         diff = LOG10(PSR_SPINB) - LOG10(PSR_SPINA)
         ! Generating random number in NBODY6++GPU
         log10P = LOG10(PSR_SPINA)+RAN2(KDUM)*diff
         ! Generating random number in MOCCA
         ! log10P = LOG10(PSR_SPINA)+ran3()*diff
         PULSARINITP = 10.d0 ** log10P
         ! write(6,*)"PulsarINITP ",PULSARINITP,PSR_SPINA,PSR_RAD

      ELSE
        ! Warning for NBODY6++GPU
        WRITE(*,*) 'PSR Invalid PSR_PMODE provided'
        ! Warning for MOCCA
        ! WRITE(*,*) 'PSR Invalid PULSAR_INITIAL_SPIN_MODEL provided'
        ! should we end program at this point and raise an error message?
      ENDIF

      ! check and calculate pulsar magnetic field at birth accordingly in Gauss.
      if (PSR_BMODE.EQ.1) then
         ! Set birth magnetic field to a fixed constant value
         ! that's set by PSR_BMAGA
         log10B = PSR_BMAGA
         PULSARINITB = 10.d0 ** log10B

      else if (PSR_BMODE.EQ.2) then
         ! Set birth magnetic field to a Uniform distribution
         ! mean as PSR_BMAGA, stdev as PSR_BMAGB
         bmaga10 = 10.0d0 ** PSR_BMAGA
         bmagb10 = 10.0d0 ** PSR_BMAGB
         diff = bmagb10 - bmaga10
         ! Generating random number in NBODY6++GPU
         PULSARINITB = bmaga10+ RAN2(KDUM)*diff
         ! Generating random number in MOCCA
         ! PULSARINITB = bmaga10+ ran3()*diff

      else if (PSR_BMODE.EQ.4) then
         ! Set to Flat-in-log distribution
         ! PSR_BMAGB set as maximum, and PSR_BMAGA set as minimum
         diff=PSR_BMAGB - PSR_BMAGA
         ! Generating random number in NBODY6++GPU
         log10B = PSR_BMAGA+RAN2(KDUM)*diff
         ! Generating random number in MOCCA
         ! log10B = PSR_BMAGA+ran3()*diff
         PULSARINITB = 10.d0 ** log10B

      else if (PSR_BMODE.EQ.5) then
         ! Log-normal distribution, from Faucher-Giuere&Kaspi 2006 (https://iopscience.iop.org/article/10.1086/501516)
         ! PSR_BMAGA as the mean and PSR_BMAGB as stdev
         log10B = -1.0d0
         do while (log10B.LT.PSR_BMIN)
            log10B = PSR_BMAGA + PSR_BMAGB * GAUS_RAND(0.0d0, 1.0d0)
         end do
         PULSARINITB = 10.d0 ** log10B
         ! write(6,*)"PulsarINITB ",log10B,PULSARINITB,PSR_RAD

      else
         ! Warning for NBODY6++GPU
         WRITE(*,*) 'PSR Invalid PSR_BMODE provided'
         ! Warning for MOCCA
         ! WRITE(*,*) 'PSR Invalid PULSAR_INITIAL_B_MODEL provided'
      endif
      endif

      if(ISAIC.eqv..true.)then
*     Birth properties of PSRMs formed through AICs.

      IF(PSRM_PMODE.EQ.1) THEN
         ! Set NS initial spin to a
         ! fixed constant value as indicated by PSRM_SPINA
         PULSARINITP = PSRM_SPINA

      ELSEIF(PSRM_PMODE.EQ.2) THEN
         ! Set NS initial spin to Uniform distribution,
         ! minimum as PSRM_SPINA, and maximum as PSRM_SPINB
         diff = PSRM_SPINB - PSRM_SPINA
         ! using random number generator in NBODY6++GPU
         PULSARINITP = PSRM_SPINA + RAN2(KDUM) * diff
         ! using random number generator in MOCCA
         ! PULSARINITP = PSRM_SPINA + ran3() * diff

      ELSEIF (PSRM_PMODE.EQ.3) THEN
         ! Set to Normal distribution
         ! mean is PSRM_SPINA, stdev is PSRM_SPINB
         PULSARINITP = -1.0d0
         do while (PULSARINITP < 0.0d0)
            PULSARINITP = GAUS_RAND(PSRM_SPINA, PSRM_SPINB)
         end do

      ELSEIF (PSRM_PMODE.EQ.4) THEN
         ! Set to Flat-in-log distribution
         ! minimum as PSRM_SPINA, maximum as PSRM_SPINB
         diff = LOG10(PSRM_SPINB) - LOG10(PSRM_SPINA)
         ! Generating random number in NBODY6++GPU
         log10P = LOG10(PSRM_SPINA)+RAN2(KDUM)*diff
         ! Generating random number in MOCCA
         ! log10P = LOG10(PSRM_SPINA)+ran3()*diff
         PULSARINITP = 10.d0 ** log10P

      ELSE
        ! Warning for NBODY6++GPU
        WRITE(*,*) 'PSR Invalid PSRM_PMODE provided'
        ! Warning for MOCCA
        ! WRITE(*,*) 'PSR Invalid PULSAR_INITIAL_SPIN_MODEL provided'
        ! should we end program at this point and raise an error message?
      ENDIF

      ! check and calculate pulsar magnetic field at birth accordingly in Gauss.
      if (PSRM_BMODE.EQ.1) then
         ! Set birth magnetic field to a fixed constant value
         ! that's set by PSRM_BMAGA
         log10B = PSRM_BMAGA
         PULSARINITB = 10.d0 ** log10B

      else if (PSRM_BMODE.EQ.2) then
         ! Set birth magnetic field to a Uniform distribution
         ! mean as PSRM_BMAGA, stdev as PSRM_BMAGB
         bmaga10 = 10.0d0 ** PSRM_BMAGA
         bmagb10 = 10.0d0 ** PSRM_BMAGB
         diff = bmagb10 - bmaga10
         ! Generating random number in NBODY6++GPU
         PULSARINITB = bmaga10+ RAN2(KDUM)*diff
         ! Generating random number in MOCCA
         ! PULSARINITB = bmaga10+ ran3()*diff

      else if (PSRM_BMODE.EQ.4) then
         ! Set to Flat-in-log distribution
         ! PSRM_BMAGB set as maximum, and PSRM_BMAGA set as minimum
         diff=PSRM_BMAGB - PSRM_BMAGA
         ! Generating random number in NBODY6++GPU
         log10B = PSRM_BMAGA+RAN2(KDUM)*diff
         ! Generating random number in MOCCA
         ! log10B = PSRM_BMAGA+ran3()*diff
         PULSARINITB = 10.d0 ** log10B

      else if (PSRM_BMODE.EQ.5) then
         ! Log-normal distribution
         ! PSRM_BMAGA as the mean and PSRM_BMAGB as stdev
         log10B = -1.0d0
         do while (log10B.LT.PSR_BMIN)
            log10B = PSRM_BMAGA + PSRM_BMAGB * GAUS_RAND(0.0d0, 1.0d0)
         end do
         PULSARINITB = 10.d0 ** log10B
         ! write(6,*)"PulsarINITB ",log10B,PULSARINITB,PSR_RAD

      else
         ! Warning for NBODY6++GPU
         WRITE(*,*) 'PSR Invalid PSRM_BMODE provided'
         ! Warning for MOCCA
         ! WRITE(*,*) 'PSR Invalid PULSAR_INITIAL_B_MODEL provided'
      endif
      endif


* Calculate initial spin-down

      MOI = NSMOI(PSRMASS, PSR_RAD)

      ! Calculate spin down rate, need to check units.

      radius_6 = (PSR_RAD*KM_TO_M) ** 6

      bfield_2 = (PULSARINITB * G_TO_T)**2
      const1 = 8.0 * pi * pi
      const2  = 3.0d0 * mu_0 * c * c * c
      pDTop = const1 * radius_6 * bfield_2
      pDBot = const2 * MOI * PULSARINITP
      PULSARINITPD = pDTop / pDBot
      PULSARINITB = PULSARINITB * G_TO_T

      IDUM1 = KDUM
      CONTAINS

      FUNCTION NSMOI(PSRMASS, PSR_RAD)
         ! Calculate moment of inertia of a neutron star
         IMPLICIT NONE
         REAL(kind=8), INTENT(IN) :: PSRMASS, PSR_RAD
         REAL(kind=8) :: M_R, NSMOI
         M_R = PSRMASS / PSR_RAD

         NSMOI = 0.237d0 * (PSRMASS*MSOL_TO_KG) * ((PSR_RAD*KM_TO_M)**2)

         NSMOI = NSMOI * (1.0d0 + (4.2d0 * M_R)
     &    + 90.0d0 * (M_R**4))

      END FUNCTION NSMOI

      FUNCTION GAUS_RAND(mu, stdev)
         IMPLICIT NONE
         ! Use the Box-Muller transform to create random Gaussian numbers.
         REAL(KIND=8), INTENT(IN) :: mu, stdev
         REAL(KIND=8) :: GAUS_RAND
         REAL(KIND=8) :: u1, u2

         ! Generating random number in NBODY6++GPU
         u1 = RAN2(KDUM)
         u2 = RAN2(KDUM)
         IF (u2.LE.0.0d0) u2 = 1.0d-12
         ! Generating random number in MOCCA
         ! u1 = ran3()
         ! u2 = ran3()

         GAUS_RAND = mu + stdev * SQRT(-2.0 *
     &    LOG(u2)) * COS(2.0 * pi * u1)

      END FUNCTION GAUS_RAND
      END SUBROUTINE PULSARINIT

* ====================================================================================
      SUBROUTINE PSRREG(NAMEI, MASSMS, KWI, TNOW, ISAICR, REGISTD)
* Register a newly-formed neutron star as a pulsar. Idempotent: a
* no-op (REGISTD=.FALSE.) if NAMEI is already registered in NAMENS.
* Single choke point for every lifecycle hook that can create a PSR
* object (kick.F AIC path, mdot.F ordinary core-collapse/merger
* paths), so the NSCOUNT<NMAX capacity check and dedup search live in
* exactly one place.
*
* Input:
* NAMEI  - NAME() of the star becoming a neutron star.
* MASSMS [Msun] - its current mass.
* KWI    - its current stellar type (KW/KSTAR), stored as NSTYPE.
* TNOW [Myr] - current physical time (TPHYS or equivalent), stored as
*   AGENSX; must be on the same time base later differenced against
*   TPHYS by the spin-down hook.
* ISAICR - true if formed via accretion-induced collapse (selects the
*   PSRM_* birth distribution in PULSARINIT; see KZ(29)=2).
*
* Output:
* REGISTD - .TRUE. if a new NSCOUNT slot was actually used.
*
      INCLUDE 'common6.h'
#ifdef PARALLEL
#define MPIINIT 1
#else
#define MPIINIT 0
#endif
      INTEGER NAMEI, KWI, PIND
      REAL*8 MASSMS, TNOW, BINIT, PINIT, PDINIT
      LOGICAL ISAICR, REGISTD, FOUND

      REGISTD = .FALSE.
      FOUND = .FALSE.
      DO PIND = 1, NSCOUNT
         IF (NAMENS(PIND).EQ.NAMEI) THEN
            FOUND = .TRUE.
            EXIT
         ENDIF
      END DO
      IF (FOUND) RETURN

      IF (NSCOUNT.GE.NMAX) THEN
         WRITE(6,*)'FATAL ERROR: NSCOUNT reached NMAX registering',
     &        ' NS, NAME=', NAMEI
#if MPIINIT
         CALL MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
         STOP
#endif
      ENDIF

      NSCOUNT = NSCOUNT + 1
      CALL PULSARINIT(BINIT, PINIT, PDINIT, MASSMS, ISAICR)
      NAMENS(NSCOUNT) = NAMEI
      XMNS(NSCOUNT) = MASSMS
      NSTYPE(NSCOUNT) = KWI
      AGENS0(NSCOUNT) = 0.0D0
      AGENSX(NSCOUNT) = TNOW
      BMAGNS0(NSCOUNT) = 0.0D0
      XMNS0(NSCOUNT) = 0.0D0
      NSSTAT(NSCOUNT) = 1
      BMAGNS(NSCOUNT) = BINIT
      PERIODNS(NSCOUNT) = PINIT
      PDOTNS(NSCOUNT) = PDINIT
      REGISTD = .TRUE.
      CALL PSREVENT(1, NSCOUNT)

      RETURN
      END

      SUBROUTINE PSREVENT(ICODE, NSINDEX)
* Write one pulsar.204 event-bank record for the neutron star at
* NAMENS/XMNS/BMAGNS/.../PDOTNS slot NSINDEX. Single choke point for
* every lifecycle hook (PSRREG registration, isolated spin-down, RLOF
* spin-up, CE evolution) so the field list, FORMAT and the
* KZ(50)/rank guard cannot drift out of sync between call sites.
*
* ICODE:
*   1 - new NS registered (PSRREG)
*   2 - isolated spin-down step (mdot.F)
*   3 - stable RLOF spin-up step (roche.f)
*   4 - common-envelope step (roche.f)
*
      INCLUDE 'common6.h'
      INCLUDE 'pulsar_ctrl.h'
      INTEGER ICODE, NSINDEX

      IF (KZ(50).NE.1.OR.rank.NE.0) RETURN

      WRITE (204,61) ICODE, TTOT*TSTAR, NAMENS(NSINDEX),
     &     XMNS(NSINDEX), BMAGNS(NSINDEX), PERIODNS(NSINDEX),
     &     PDOTNS(NSINDEX), AGENSX(NSINDEX), AGENS0(NSINDEX),
     &     BMAGNS0(NSINDEX), XMNS0(NSINDEX), PSR_TAU
   61 FORMAT('PULSAR ',I2,1P,E17.9,I10,9E17.9)

      RETURN
      END
