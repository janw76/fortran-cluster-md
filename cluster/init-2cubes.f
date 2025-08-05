*===============================================================================
*	INITIALIZE STARTING CONFIGURATION
*===============================================================================
* file: init-2cubes.f
*-------------------------------------------------------------------------------
* This routine sets two latices in tho box, one in the left half, one in
* the right. One lattice has a size of 10^3, the other 5^3
* FOR TESTING PURPOSES ONLY
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
	integer index,nat1,nat2
	double precision urx,ury,urz
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

	urx=boxx/4.0d0
	ury=boxy/2.0d0
	urz=boxz/2.0d0

	nat1=10

	do ix=0,nat1-1
	  do iy=0,nat1-1
	    do iz=0,nat1-1
	      index=ix*nat1*nat1+iy*nat1+iz+1
	      x(index)=real(ix-int(nat1/2))*atr*1.5d0+urx
	      y(index)=real(iy-int(nat1/2))*atr*1.5d0+ury
	      z(index)=real(iz-int(nat1/2))*atr*1.5d0+urz
	      vx(index)=(kpsRND()-0.5d0)
	      vy(index)=(kpsRND()-0.5d0)
	      vz(index)=(kpsRND()-0.5d0)
	    enddo
	  enddo
	enddo

	urx=3.0d0*boxx/4.0d0
	ury=boxy/2.0d0
	urz=boxz/2.0d0

	nat2=0

	do ix=0,nat2-1
	  do iy=0,nat2-1
	    do iz=0,nat2-1
	      index=ix*nat2*nat2+iy*nat2+iz+1+nat1*nat1*nat1
	      x(index)=real(ix-int(nat2/2))*atr*1.5d0+urx
	      y(index)=real(iy-int(nat2/2))*atr*1.5d0+ury
	      z(index)=real(iz-int(nat2/2))*atr*1.5d0+urz
	      vx(index)=(kpsRND()-0.5d0)
	      vy(index)=(kpsRND()-0.5d0)
	      vz(index)=(kpsRND()-0.5d0)
	    enddo
	  enddo
	enddo

	natom=nat1*nat1*nat1+nat2*nat2*nat2

	call stopdead(-1)
	call neighbor(0)	
	call AccAtom()
	call AccWall()

	return
	end
*-------------------------------------------------------------------------------
