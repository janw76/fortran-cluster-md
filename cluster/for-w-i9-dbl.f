*===============================================================================
*	ACCELERATION ATOM-WALL
*===============================================================================
* file: for-w-i9-dbl.f
*-------------------------------------------------------------------------------
* In this routine the influence of the walls on the simulated particles is
* calculated. The origin of the coordinate system in placed in the front
* lower left corner of the box.
* This routine is written for a fixed 1/r^9 potential. (some optimization
* potential existing)
* This routine additionally has a seperating wall at (x/2) with an
* ajustable hole in it. This is for checking an adiabatic expansion!
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
	double precision hx,hy,hz,mx,my,mz,r,rr,dx,uw,atri2
	double precision ex0,ex1,ex11,ex2,ex21,fwr,fwx

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

* f(x)=(1/(1+exp((x+x0)/dx)))+(1/(1+exp((x0-x)/dx)))


	atri2=1.0d0/(atr*atr)
	c1=atr*atr*atr*atr*atr*atr*atr*atr*atr
	c2=c1*beng*9.0d0
	ewall=0.0d0
	fwall=0.0d0

	r=60.0d0
	dx=2.0d0
	hx=boxx/2.0d0
	hy=boxy/2.0d0
	hz=boxz/2.0d0

	do i=1,natom
	  dxl=1.0d0/x(i)
	  dxr=1.0d0/(boxx-x(i))
	  dyl=1.0d0/y(i)
	  dyr=1.0d0/(boxy-y(i))
	  dzl=1.0d0/z(i)
	  dzr=1.0d0/(boxz-z(i))

	  mx=x(i)-hx
	  my=y(i)-hy
	  mz=z(i)-hz
	  rr=dsqrt(my*my+mz*mz)

*	  ex0=exp(-abs(mx)/atr)
	  ex0=exp(-mx*mx*atri2)
	  ex1=exp((rr+r)/dx)
	  ex11=1.0d0+ex1
	  ex2=exp((r-rr)/dx)
	  ex21=1.0d0+ex2

	  uw =100.0d0*beng*ex0*(1.0d0/ex11 + 1.0d0/ex21)
	  fwr=100.0d0*beng*ex0*(ex2/(ex21*ex21) - ex1/(ex11*ex11))/dx
	  fwx=2.0d0*mx*atri2*uw
*	  fwx=2.0d0*mx*uw/atr

	  if (rr.lt.1.0d-2) then
	    rr=1.0d-2
	    my=0.0d0
	    mz=0.0d0
	  endif

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

	  ax(i)=ax(i)+c2*(fwxl-fwxr)/atm+fwx/atm
          ay(i)=ay(i)+c2*(fwyl-fwyr)/atm-fwr*my/(atm*rr)
          az(i)=az(i)+c2*(fwzl-fwzr)/atm-fwr*mz/(atm*rr)

	  fwall=fwall+c2*(fwxl+fwxr+fwyl+fwyr+fwzl+fwzr)
	  ewall=ewall+dxl9+dxr9+dyl9+dyr9+dzl9+dzr9
	enddo
	ewall=ewall*beng*c1
	epot=epot+ewall
	fwallsum=fwallsum+fwall
	ewallsum=ewallsum+ewall

*	if(mod(step,5).eq.0) then
*	write (42,10) step,x(2),vx(2),ax(2),fwx/atm,
*     $                     y(2),vy(2),ay(2),fwr*my/(atm*rr),
*     $                     z(2),vz(2),az(2),fwr*mz/(atm*rr)
*	endif
*10	format (I10,3(2X,F15.5,3(2X,F15.12)))

	return
	end
*-------------------------------------------------------------------------------
