      SUBROUTINE TIDES3(SEP,M1,M2,ECC,VSTAR,DTIME,DE)
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
      REAL*8  DE(5)
*
*
    
      C = 3.0D+05/VSTAR
      ECC2 = ECC*ECC
      ECC4 = ECC2*ECC2
      MTOT = M1 + M2
      COST = C**(-5) * M1 * M2 * MTOT


      FE = (1.0D+0 + 7.3D+1/2.4D+1*ECC2 + 3.7D+1/9.6D+1*ECC4)
      FE1 = (1 - ECC2)**(-3.5)*FE
      FE2 = ECC*(1 - ECC2)**(-2.5)*(1.0D+0 + 1.21D+2/3.04D+2*ECC2)
      FE3 = (1 - ECC2)**(-2.0)*(1.0D+0 + 7.0D+0/8.0D+0*ECC2)

*       DENERGY due to GW
      DE(1) = 3.2D+1/5.0D+0*COST*M1*M2*SEP**(-5.0)*FE1*DTIME
*
      DE(2) = 0.0
*       DSEP due to GW
      DE(3) = 6.4D+1/5.0D+0*COST*SEP**(-3)*FE1*DTIME 
*       DECC due to GW
      DE(4) = 3.04D+2/1.5D+1*COST*SEP**(-4)*FE2*DTIME
*       DJ due to GW
      DE(5) = 3.2D+1/5.0D+0*COST*M1*M2*(M1 + M2)**(-0.5)*SEP**(-3.5)*FE3
      DE(5) = DE(5)*DTIME
       
*      PRINT *, 'TIDES3 41 ',ECC 
*

      RETURN
*
      END
