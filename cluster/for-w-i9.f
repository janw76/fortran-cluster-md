*===============================================================================
*	ACCELERATION ATOM-WALL
*===============================================================================
* file: for-w-i9.f
*-------------------------------------------------------------------------------
* In this routine the influence of the walls on the simulated particles is
* calculated. The origin of the coordinate system in placed in the front
* lower left corner of the box.
* This routine is written for a fixed 1/r^9 potential. (some optimization
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
	double precision dxl,dxr,dyl,dyr,dzl,dzr,c1,c2
	double precision dxl9,dxr9,dyl9,dyr9,dzl9,dzr9
	double precision fwxl,fwxr,fwyl,fwyr,fwzl,fwzr
	
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

	c1=atr*atr*atr*atr*atr*atr*atr*atr*atr
	c2=c1*beng*9.0d0
	ewall=0.0d0
	fwall=0.0d0

	do i=1,natom
	  dxl=1.0d0/x(i)
	  dxr=1.0d0/(boxx-x(i))
	  dyl=1.0d0/y(i)
	  dyr=1.0d0/(boxy-y(i))
	  dzl=1.0d0/z(i)
	  dzr=1.0d0/(boxz-z(i))

	  dxl9=dxl*dxl*dxl*dxl*dxl*dxl*dxl*dxl*dxl
	  dxr9=dxr*dxr*dxr*dxr*dxr*dxr*dxr*dxr*dxr
	  dyl9=dyl*dyl*dyl*dyl*dyl*dyl*dyl*dyl*dyl
	  dyr9=dyr*dyr*dyr*dyr*dyr*dyr*dyr*dyr*dyr
	  dzl9=dzl*dzl*dzl*dzl*dzl*dzl*dzl*dzl*dzl
	  dzr9=dzr*dzr*dzr*dzr*dzr*dzr*dzr*dzr*dzr
	  
	  fwxl=dxl9*dxl
          fwxr=dxr9*dxr
          fwyl=dyl9*dyl
          fwyr=dyr9*dyr
          fwzl=dzl9*dzl
          fwzr=dzr9*dzr

	  ax(i)=ax(i)+c2*(fwxl-fwxr)/atm
          ay(i)=ay(i)+c2*(fwyl-fwyr)/atm
          az(i)=az(i)+c2*(fwzl-fwzr)/atm

	  fwall=fwall+c2*(fwxl+fwxr+fwyl+fwyr+fwzl+fwzr)
	  ewall=ewall+dxl9+dxr9+dyl9+dyr9+dzl9+dzr9
	enddo
	ewall=ewall*beng*c1
	epot=epot+ewall
	fwallsum=fwallsum+fwall
	ewallsum=ewallsum+ewall
	return
	end
*-------------------------------------------------------------------------------
