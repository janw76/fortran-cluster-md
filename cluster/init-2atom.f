*===============================================================================
*	INITIALIZE STARTING CONFIGURATION
*===============================================================================
* file: init-2atom.f
*-------------------------------------------------------------------------------
* This specialized routine sets two atoms next to each other on a (more or
* less) collision course. To alter the positions, the source code has to
* be edited... INCONVENIENT, ONLY FOR TESTING PURPOSES!
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

	double precision urx,ury,urz
	
*-------------------------------------------------------------------------------
* urx,ury,urz	: translation vector to the box center
*-------------------------------------------------------------------------------

	urx=boxx/2.0d0
	ury=boxy/2.0d0
	urz=boxz/2.0d0

	x(1)=urx+10.0
	y(1)=ury+0.0
	z(1)=urz
	vx(1)=1.0
	vy(1)=0.0d0
	vz(1)=0.0d0

	x(2)=urx-10.0
	y(2)=ury-0.0
	z(2)=urz-0.0
	vx(2)=-1.0
	vy(2)=0.0d0
	vz(2)=0.0d0

	x(3)=urx-2.0
	y(3)=ury-0.0
	z(3)=urz-0.0
	vx(3)=0.0d0
	vy(3)=0.0d0
	vz(3)=0.0d0

	x(4)=urx+2.0
	y(4)=ury+0.0
	z(4)=urz+0.0
	vx(4)=0.0d0
	vy(4)=0.0d0
	vz(4)=0.0d0
	natom=4

	call stopdead(-1)
	call neighbor(0)	
	call AccAtom()
	call AccWall()

	return
	end
*-------------------------------------------------------------------------------
