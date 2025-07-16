! from mocca Lx.f90
FUNCTION get_Lx(Kacc, Racc, Macc, dMdonor) RESULT(Lx)
    ! Wrapper function around different models of X-ray luminosity calculations
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: Kacc, Racc, Macc, dMdonor
    REAL(KIND=8) :: Lx

    INTEGER LX_MODEL
    LX_MODEL = 0  ! later make a switch or input; temporarily hard-coded to 0

    IF (LX_MODEL == 0) THEN
        Lx = get_Lx_startrack(Kacc, Racc, Macc, dMdonor)
    ELSEIF (LX_MODEL == 1) THEN
        Lx = get_Lx_SS(Macc, dMdonor)
    ELSE
        write(*,*) "Wrong value in get_Lx()", LX_MODEL
    ENDIF

    CONTAINS

        FUNCTION get_Lx_startrack(Kacc, Racc, Macc, dMdonor) RESULT(Lx)
          ! Calculates X-ray luminosity from accretion onto NS/BH as in startrack code
          ! Args:
          !   Kacc 
          !   Racc [Rsun]
          !   Macc [Msun]
          !   dMdon [Msun/yr]
          ! Return:
          !   [erg/s] - X-ray luminosity
          !
          !  based on startrack code
          IMPLICIT NONE
          REAL(KIND=8), INTENT(IN) :: Kacc, Racc, Macc, dMdonor
          REAL(KIND=8) :: Lx
          REAL(KIND=8) :: Rfactor, etaA
          !REAL(KIND=8) :: Ys, Lxdonor, GG, Rsun, Msun, etaA, Rfactor

          ! Constants
          REAL(KIND=8), PARAMETER :: GG = 3.93d+08  !2938.408721d0 * (365.25d0**2.0) 
                                                    ! [Rsun^(3) Msun^(-1) yr^(-2)]
          REAL(KIND=8), PARAMETER :: Rsun = 6.957d+10                      ! [cm], sun radius in cgs units
          REAL(KIND=8), PARAMETER :: Msun = 1.988d+33                      ! [g], sun mass in cgs units
          REAL(KIND=8), PARAMETER :: Ys = 31557600.0d0  ! Year in seconds

          ! Check Kb
          IF (Kacc == 14) THEN
            etaA = 0.5d0                        ! Disk accretion onto BH
            Rfactor = 3.0d0                   ! Disk ends at 3 Schwarzschild radii of BH
          ELSE
            etaA = 1.0d0                        ! Surface accretion onto NS
            Rfactor = 1.0d0
          END IF

          ! Calculate Lxdonor (not Eddington limited X-ray luminosity)
          ! write (*, *) "Lx0", Kacc, Racc, Macc, dMdonor
          ! write (*, *) "Lx00", etaA, Rfactor
          Lx = etaA * GG * Macc * dMdonor / (Rfactor * Racc)
          ! write (*, *) "Lx1", Lx
          ! write (*, *) "Lx2", (Rsun * Rsun * Msun / (Ys**3.0d0))
          Lx = Lx * (Rsun * Rsun * Msun / (Ys**3.0d0))  ! [erg/s]
          ! write (*, *) "Lx3", Lx
          
          !get_Lx = Lx

        END FUNCTION get_Lx_startrack

        FUNCTION get_Lx_SS(M_acc, Mdot) RESULT(Lx)
            ! Calculate the X-ray luminosity from the accretion disk following Shakura & Sunyaev (1973)
            !    see also Lasota & King (2023) for general introduction
            !   
            !   Note: It's the accretion luminosity and may not be equal to the apparent
            !   luminosity if the emission is not isotropic. This can happen e.g. due to
            !   beaming and need to be accounted for in post-processing
            !
            !   FIXME: the calculation of eddington limit follows the one from Lasota &
            !   King (2023) but it's different to general BSE formula in
            !   get_Eddington_Luminosity()
            !
            ! Args:
            !   Macc [Msun] : accretor mass
            !   Mdot [Msun/yr] : mass transfer rate into the disk (not necessirly the
            !       accretion onto the compact object which may be lower by the disk winds>)
            ! Returns:
            !   [erg/s] : the X-ray luminosity of the accretion disk 

            IMPLICIT NONE
            REAL(kind=8), INTENT(IN) :: M_acc, Mdot
            REAL(kind=8) :: Lx
            REAL(kind=8) :: eta01, Mdot_Edd, L_Edd, Mdot_by_Mdot_Edd

            eta01 = 1  ! FIXME: find the exact meaning of eta01 - it's not defined in Lasota & King (2023) 
            Mdot_Edd = 2.2e-8 * M_acc / eta01  ! Eddinton limit [Msun/yr]  eq. 3 Lasota & King (2023)
            L_Edd = 5.66295811d45 * Mdot_Edd  ! Eddinton luminosity [erg/s]
                                              ! the constant was calculated as (0.1 * c.c**2 * u.Msun / u.yr).to('erg/s') using astropy
            Mdot_by_Mdot_Edd = Mdot / Mdot_Edd  ! descriptive name because FORTRAN is case insensitive and mdot cannot be used here
            IF (Mdot_by_Mdot_Edd > 1) THEN
                Lx = L_Edd * ( 1 + log(Mdot_by_Mdot_Edd))
            ELSE
                Lx = L_Edd * Mdot_by_Mdot_Edd
            ENDIF

        END FUNCTION get_Lx_SS

END FUNCTION get_Lx
