*===============================================================================
*	INITIALIZE STARTING CONFIGURATION
*===============================================================================
* file: init-cube1.f
* date: 09-nov-2001
*-------------------------------------------------------------------------------
* This routine is only partially funtional and has some limitations
* it sets a^3 atoms in a cubic primitive lattice, sets the total momentum
* to zero and scales the kinetic temperature to the desired value.
* Only system sizes with correct cubic values are permissible!! If illegal
* values are set in the input file, the nearest correct value is set and
* used.
* 04/10/2007(JW):
*   - changed to accept arbitrary N and maximum spacing in the box 
*   - also randomly swaps particles at the start if there is more than one sort
*   - performing some warm-up steps before the actual MD main run
*-------------------------------------------------------------------------------
* notes on units of measure inside the program:
* length (l)	 : in  A (10^-10 m)
* time   (t)  	 : in fs (10^-15 s)
* mass   (m)	 : in kg (10^-27 kg)
* 		  Input in in amu!!  (1.6605655 10^-27 kg)
* temperature (T): in K
* Energy      (J): in m*l^2*t^-2
*		   = 10^-17 J
*-------------------------------------------------------------------------------
* Parameters for argon:
* Diameter (sigma)	: 0.3405 nm
* Energy-well (epsilon)	: 1.6567944 10^-21 J (120K)
* Parameters Kari	: 3.400 A, 120.7728 K
*-------------------------------------------------------------------------------

	subroutine initconf()
	
	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'boxpp.inc'
	include 'energ.inc'

	integer ix,iy,iz,i,a,b
	integer index,nat
	double precision urx,ury,urz,temp
	double precision kpsRND
	
*-------------------------------------------------------------------------------
* ix,iy,iz	: loop variables for cube coordinates
* index		: index for atom in the cube, calculated from ix,iy,iz
* nat		: side length of the argon-cube
* urx,ury,urz	: translation vector to the box center
*-------------------------------------------------------------------------------

	if (rndctrl.eq.'X') then
	  call kpsinit(0)
	  write (*,*) 'Using random velocities'
	else
	  call kpsinit(1)
	  write (*,*) 'Using static velocities'
	endif

	urx=boxx/2.0d0
	ury=boxy/2.0d0
	urz=boxz/2.0d0

	nat=nint(real(natom)**(1.0d0/3.0d0)+0.5d0) 
*	if (natom.ne.nat*nat*nat) then
*	  write (*,*) 'WARNING --- Number of Atoms changed to',nat*nat*nat
*	endif

	if (natom.gt.mxa) then
	  write (*,*) 'ERROR --- Too many atoms requested'
	  write (*,*) 'Requested:',natom*natom*natom
	  write (*,*) 'Available:',mxa
	  stop
	endif
	if ((boxx-2.0*atr)/(DBLE(nat)-1.0d0).lt.atr) then 
	   write (*,*)
	   temp = ((boxx-2.0*atr)/(DBLE(nat)-1.0d0))/atr
	write (*,'(A,F3.1)') '-> CAREFUL! Very dense! Dist/sig less than ',temp
	   write (*,*)
	endif

	do ix=0,nat-1
	 do iy=0,nat-1
	  do iz=0,nat-1
	   index=ix*nat*nat+iy*nat+iz+1
	   if (index.eq.natom+1) goto 10
	   x(index)=atr+real(ix)*((boxx-2.0*atr)/(DBLE(nat)-1.0))
	   y(index)=atr+real(iy)*((boxy-2.0*atr)/(DBLE(nat)-1.0))
	   z(index)=atr+real(iz)*((boxz-2.0*atr)/(DBLE(nat)-1.0))
	   vx(index)=(kpsRND()-0.5d0)
	   vy(index)=(kpsRND()-0.5d0)
	   vz(index)=(kpsRND()-0.5d0)
	  enddo
	 enddo
	enddo
	
10	continue


*** Swap some particles	if more than one type
	if (atsorts.ge.2) then
	  do i=1,10*natom
	    a = INT(kpsRND()*natom)+1
	    b = INT(kpsRND()*natom)+1
	    if (a.ne.b) call swap(a,b)
	  enddo
	endif
	
*** Initalize everything	
	call stopdead(-1)
	call neighbor(0)
	call AccAtom()
	call AccWall()
	
*** A few warmup steps...
	write (*,*) 'doing 100 warmup steps...'
	do i=1,100
	   call integrator(i)
	enddo

*** Resetting center-of-mass movement again, just to be sure...
	call stopdead(-1)
	
	return
	end
	
*-------------------------------------------------------------------------------	
*** SUBROUTINES
*** SWAP 

	subroutine swap(a,b)
	
	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'boxpp.inc'
	include 'energ.inc'
	
	integer a,b
	real tx,ty,tz
	
	tx = x(a)
	ty = y(a)
	tz = z(a)
	
	x(a) = x(b)
	y(a) = y(b)
	z(a) = z(b)
	
	x(b) = tx
	y(b) = ty
	z(b) = tz

	return
	end	
	
*-------------------------------------------------------------------------------
