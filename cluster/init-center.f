*===============================================================================
*	INITIALIZE STARTING CONFIGURATION
*===============================================================================
* file: init-center.f
* date: 09-nov-2001
*-------------------------------------------------------------------------------
* Sets the first particle right in the center and the rest on an cube of optimal
* length
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

	integer ix,iy,iz
	integer index,nat
	double precision urx,ury,urz
	double precision distancex,distancey,distancez,halfcube
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

	nat=nint(real(natom)**(1.0d0/3.0d0)) 

	distancex = REAL((boxx-1.0)/nat)	
	if (distancex.lt.atr) write (*,*) 'Carful! Very dense!'
	distancey = REAL((boxx-1.0)/nat)	
	if (distancey.lt.atr) write (*,*) 'Carful! Very dense!'
	distancez = REAL((boxx-1.0)/nat)	
	if (distancez.lt.atr) write (*,*) 'Carful! Very dense!'

	if (nat*nat*nat.gt.mxa) then
	  write (*,*) 'ERROR --- Too many atoms requested'
	  write (*,*) 'Requested:',nat*nat*nat
	  write (*,*) 'Available:',mxa
	  stop
	endif

	halfcube = (boxx-1.2*atrcut)/2

	do ix=0,nat-1
	  do iy=0,nat-1
	    do iz=0,nat-1
	      index=ix*nat*nat+iy*nat+iz+1
	      if (index.eq.natom+1) goto 10
	      x(index)=real(ix)*distancex+0.5
	      y(index)=real(iy)*distancey+0.5
	      z(index)=real(iz)*distancez+0.5
	      vx(index)=(kpsRND()-0.5d0)
	      vy(index)=(kpsRND()-0.5d0)
	      vz(index)=(kpsRND()-0.5d0)
	    enddo
	  enddo
	enddo
	
10	continue

	call stopdead(-1)
	call stopdead(-1)
	call neighbor(0)
	call AccAtom()
	call AccWall()

	return
	end
*-------------------------------------------------------------------------------
