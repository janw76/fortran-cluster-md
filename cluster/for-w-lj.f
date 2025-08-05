*===============================================================================
*	ACCELERATION ATOM-WALL
*===============================================================================
* file: for-w-lj.f
*-------------------------------------------------------------------------------
* In this routine the influence of the walls on the simulated particles is
* calculated. The origin of the coordinate system in placed in the front
* lower left corner of the box.
* There are two LJ-Walls at x and x+boxx, the other walls are missing
*-------------------------------------------------------------------------------

	subroutine AccWall()
	
	implicit none

	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'boxpp.inc'
	include 'energ.inc'
	
	integer i
	double precision c1,c2,c3
	double precision dxl,dxr,dxr2,dxl2,dxl6,dxr6
	double precision fwxl,fwxr,ewxl,ewxr
	
*-------------------------------------------------------------------------------
* i		: Allgemeine Schleifenvariable
* c1		: particle diameter to the ninth power
* c2		: force constant for wall potential
* dxl,dyl,dzl	: distance of the atom to the left wall in X,Y,Z-direction
* dxr,dyr,dzr	: distance of the atom to the right wall in X,Y,Z-direction
* dxl9,dyl9,dzl9: distance to the left wall, ninth power
* dxr9,dyr9,dzr9: distance to the right wall, ninth power
* fwxl,fwyl,fwzl: forces on the atoms from the left walls
* fwxr,fwyr,fwzr: forces on the atoms from the right walls
*-------------------------------------------------------------------------------

	c1=atr*atr*atr*atr*atr*atr
	c2=c1*c1
	c3=0.5d0*c1

	ewall=0.0d0
	fwall=0.0d0

	do i=1,natom
	  dxl=x(i)
	  dxr=(boxx-x(i))
	  if (dxl.le.atrcut) then
	    dxl2=1.0d0/(dxl*dxl)
	    dxl6=dxl2*dxl2*dxl2
	    fwxl=dxl*dxl2*dxl6*(c2*dxl6-c3)
	    ewxl=dxl6*(c2*dxl6-c1)
	  else
	    fwxl=0.0d0
	    ewxl=0.0d0
	  endif
	  if (dxr.le.atrcut) then
	    dxr2=1.0d0/(dxr*dxr)
	    dxr6=dxr2*dxr2*dxr2
	    fwxr=dxr*dxr2*dxr6*(c2*dxr6-c3)
	    ewxr=dxr6*(c2*dxr6-c1)
	  else
	    fwxr=0.0d0
	    ewxr=0.0d0
	  endif

	  ax(i)=ax(i)+beng*48.0d0*(fwxl-fwxr)/atm

	  fwall=fwall+fwxl+fwxr
	  ewall=ewall+ewxl+ewxr
	enddo

	ewall=ewall*beng*4.0d0
	epot=epot+ewall
	fwallsum=fwallsum+fwall
	ewallsum=ewallsum+ewall
	return
	end
*-------------------------------------------------------------------------------
