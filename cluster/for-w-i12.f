*===============================================================================
*	ACCELERATION ATOM-WALL
*===============================================================================
* file: for-w-i12.f
*-------------------------------------------------------------------------------
* In this routine the influence of the walls on the simulated particles is
* calculated. The origin of the coordinate system in placed in the front
* lower left corner of the box.
* This routine is written for a fixed 1/r^12 potential. (some optimization
* potential existing)
*-------------------------------------------------------------------------------

	subroutine AccWall()
	
	implicit none

	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'boxpp.inc'
	include 'energ.inc'
	
	integer i
	double precision dxl,dxr,dyl,dyr,dzl,dzr,c
	double precision dxl9,dxr9,dyl9,dyr9,dzl9,dzr9
	
*-------------------------------------------------------------------------------
* i		: Allgemeine Schleifenvariable
* c		: particle diameter to the twelfth power
* dxl,dyl,dzl	: distance of the atom to the left wall in X,Y,Z-direction
* dxr,dyr,dzr	: distance of the atom to the right wall in X,Y,Z-direction
* dxl9,dyl9,dzl9: distance to the left wall, twelfth power
* dxr9,dyr9,dzr9: distance to the right wall, twelfth power
*-------------------------------------------------------------------------------

	c=atr*atr*atr*atr*atr*atr*atr*atr*atr*atr*atr*atr
	do i=1,natom
	  dxl=1.0d0/x(i)
	  dxr=1.0d0/(boxx-x(i))
	  dyl=1.0d0/y(i)
	  dyr=1.0d0/(boxy-y(i))
	  dzl=1.0d0/z(i)
	  dzr=1.0d0/(boxz-z(i))
	  dxl9=dxl*dxl*dxl*dxl*dxl*dxl*dxl*dxl*dxl*dxl*dxl*dxl
	  dxr9=dxr*dxr*dxr*dxr*dxr*dxr*dxr*dxr*dxr*dxr*dxr*dxr
	  dyl9=dyl*dyl*dyl*dyl*dyl*dyl*dyl*dyl*dyl*dyl*dyl*dyl
	  dyr9=dyr*dyr*dyr*dyr*dyr*dyr*dyr*dyr*dyr*dyr*dyr*dyr
	  dzl9=dzl*dzl*dzl*dzl*dzl*dzl*dzl*dzl*dzl*dzl*dzl*dzl
	  dzr9=dzr*dzr*dzr*dzr*dzr*dzr*dzr*dzr*dzr*dzr*dzr*dzr
	  
	  epot=epot+beng*c*(dxl9+dxr9+dyl9+dyr9+dzl9+dzr9)
    
	  ax(i)=ax(i)+9d0*beng*c*((dxl9*dxl)-(dxr9*dxr))/atm
	  ay(i)=ay(i)+9d0*beng*c*((dyl9*dyl)-(dyr9*dyr))/atm
	  az(i)=az(i)+9d0*beng*c*((dzl9*dzl)-(dzr9*dzr))/atm
	enddo

	return
	end
*-------------------------------------------------------------------------------
