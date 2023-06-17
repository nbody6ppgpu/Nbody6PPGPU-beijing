	subroutine bisection(RR, I)
*
*
*	
*
*
	INCLUDE 'common6.h'
	real*8:: error,f,r,temp
	integer:: a,b,c,I
*	real*8:: radii(NMAX), massArray(NMAX)
*
	error= 1.0
	a=1
*	b=count(RMR/=0.)
	b=NParticles
*
*
*	XMR cumulated mass profile
*	RMR are the radii of the previously cumulated mass profile
*	
*	write(*,*) "Radius"
*
* 10 	read(*,*) RR
*
*
 15	temp=((RMR(a)-RR)*(RMR(b)-RR))
*
	if (temp .eq. 0) then
		if (RMR(a)-RR .eq. 0) then
		c=a
		goto 20
		else
		c=b
		goto 20
		end if
*
  	else if (temp .le. 0) then
		c=ceiling((a+b)/2.0)
	else
*		write(*,*)"a, b = ",a,b
*		write(*,*)"RMR = ", RMR(a),RMR(b)
*		write(*,*)"RR = ",RR
*		write(*,*)"Try with another values of a and b"
		goto 30
	end if
*
	if (((RMR(a)-RR)*(RMR(c)-RR)) .le. 0) then
		b=c
	else
		a=c
*
	end if
	if (abs(b-a) .gt. error) then
		goto 15
	endif
*
*
 20	if ((RR-RMR(c)) .le. 0) then
		c=c-1
	endif
*	
	if (RR .lt. RMR(100)) then
		call interpolate_powerlaw(RR, I)
*
	else if (RR .gt. RMR(NMAX-6)) then
        	call interpolate_linearly(RR,
     &			   RMR(c-5:c+5), XMR(c-5:c+5), I, 11)
*
	else	
 		call interpolate_linearly(RR,
     &                     RMR(c-5:c+5), XMR(c-5:c+5), I, 6)
	endif
*
 30	end subroutine bisection
*	
*
*
	subroutine interpolate_linearly(interRadius, RMR_sel,
     &					XMR_sel, I, J)
*
*
*	
*
*
	INCLUDE 'common6.h'
*	!interpolation between the two adjacent array values
	real*8:: interRadius, XMR_sel(11), RMR_sel(11) 
        real*8:: dR, dM
	integer:: I, J
*	
	dR = RMR_sel(11)-RMR_sel(1)
	dM = XMR_sel(11)-XMR_sel(1)
*	MASSSLOPE(I) = dM/dR
	XSLOPE(I) = dM/dR
	XINTERMASS(I)=XMR_sel(J)+XSLOPE(I)
     &   *(interRadius-RMR_sel(J))
*	DIAGNOSTICS - REMOVE LATER!
*	write(*,*)"interpol. Mass, XMR1, XMR2, slope = ",
*     &   XINTERMASS(I),XMR_sel(5:7), XSLOPE(I)
*	write(*,*)"I,J, radius, RMR1, RMR2 = ",I,J,
*     &   interRadius,RMR_sel(5:7)

	end subroutine interpolate_linearly
*
*
        subroutine interpolate_powerlaw(interRadius, I)
*
*
*
*
*
        INCLUDE 'common6.h'
	REAL*8:: delta,interRadius
	INTEGER:: I
*
	delta = log(XMR(100))/log(RMR(100))
	XINTERMASS(I)=interRadius**delta
	end subroutine interpolate_powerlaw
