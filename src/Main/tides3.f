      SUBROUTINE TIDES3(SEMI,M1,M2,ECC,CLIGHT,DTIME,DE)
*
*
*       GR tidal energy loss for interacting stars.
*       OUTPUT:
*              DE(1) energy loss per time unit
*              DE(3) separation change per time unit
*              DE(4) ecc. change pet time unit
*              DE(5) ang. momentum change per time unit	
*       -------------------------------------------
*  
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8  DE(5),CLIGHT
	  REAL, PARAMETER :: INVPI = (4.D0*DATAN(1.D0))**(-1)
*

      ECC2 = ECC*ECC
      ECC4 = ECC2*ECC2
      MTOT = M1 + M2
      COST = CLIGHT**(-5) * M1 * M2 * MTOT

* Elliptic case

	  IF (ECC.LT.1.0D+0) THEN

		FE = (1.0D+0 + 7.3D+1/2.4D+1*ECC2 + 3.7D+1/9.6D+1*ECC4)
		FE1 = (1 - ECC2)**(-3.5)*FE
		FE2 = ECC*(1 - ECC2)**(-2.5)*(1.0D+0 + 1.21D+2/3.04D+2*ECC2)
		FE3 = (1 - ECC2)**(-2.0)*(1.0D+0 + 7.0D+0/8.0D+0*ECC2)

*       DENERGY due to GW
		DE(1) = 3.2D+1/5.0D+0*COST*M1*M2*SEMI**(-5.0)*FE1*DTIME
*
		DE(2) = 0.0
*       DSEMI due to GW
		DE(3) = 6.4D+1/5.0D+0*COST*SEMI**(-3)*FE1*DTIME 
*       DECC due to GW
		DE(4) = 3.04D+2/1.5D+1*COST*SEMI**(-4)*FE2*DTIME
*       DJ due to GW
		DE(5) = 3.2D+1/5.0D+0*COST*M1*M2*(M1 + M2)**(-0.5)*SEMI**(-3.5)*FE3
		DE(5) = DE(5)*DTIME
		
* Hyperbolic case JFNS

	  ELSE IF (ECC.GT.1.0D+0) THEN
	
		FE1 = 3.0D+0*(9.6D+1 + 2.92D+2*ECC2 + 3.7D+1*ECC4)*DACOS(-1/ECC)
	    FE1 = FE1 + (6.73D+2*ECC2+6.02D+2)*DSQRT(ECC2-1.0D0)
		FE1 = INVPI*(ECC2-1)**(-3.5)*FE1
		
		FE2 = 3.0D+0*ECC2*(3.04D+2 + 1.21D+2*ECC2)*DACOS(-1/ECC)
		FE2 = FE2 + (1.34D+2+1.069D+3*ECC2+7.2D+1*ECC4)*DSQRT(ECC2-1.0D0)
		FE2 = INVPI*ECC**(-1.0)*(ECC2-1)**(-2.5)*FE2
		
		
		FE3 = (7.0D+0*ECC2+8.0D+0)*DACOS(-1/ECC)
		FE3 = FE3 + (2.0D+0*ECC2+1.3D+1)*DSQRT(ECC2-1.0D0)
		FE3 = INVPI*(ECC2-1)**(-2.0)*FE3
		
*       DENERGY due to GW
		DE(1) = 1.0D+0/4.5D+1*COST*M1*M2*DABS(SEMI)**(-5.0)*FE1*DTIME
*
		DE(2) = 0.0
*       DSEMI due to GW
		DE(3) = 2.0D+0/4.5D+1*COST*DABS(SEMI)**(-3.0)*FE1*DTIME 
*       DECC due to GW
		DE(4) = 1.0D+0/4.5D+1*COST*DABS(SEMI)**(-4)*FE2*DTIME
*       DJ due to GW
		DE(5) = 4.0D+0/5.0D+0*COST*M1*M2*MTOT**(-0.5)*DABS(SEMI)**(-3.5)
		DE(5) = DE(5)*FE3*DTIME
	  
	  END IF
      RETURN
*
      END
