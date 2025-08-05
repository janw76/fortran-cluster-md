*===============================================================================
* 	RESIZE BOX
*===============================================================================
* file: resize.f
*-------------------------------------------------------------------------------
* This routine expands or shrinks the box according to parameters
*-------------------------------------------------------------------------------

	subroutine resize(step)

	implicit none

	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'boxpp.inc'

	integer step
	integer i
	double precision xmin,xmax,ymin,ymax,zmin,zmax
	double precision dx,dy,dz
	double precision vxs,vys,vzs

*-------------------------------------------------------------------------------
* step		: Time step
*...............................................................................
* i		: standard loop variable
* xmin		: variable for the smallest X-coordinate
* ymin		: variable for the smallest Y-coordinate
* zmin		: variable for the smallest Z-coordinate
* xmax		: variable for the largest X-coordinate
* ymax		: variable for the largest Y-coordinate
* zmax		: variable for the largest Z-coordinate
* dx,dy,dz	: shifting parameters for coordinates
*-------------------------------------------------------------------------------

	xmin=boxx
	ymin=boxy
	zmin=boxz
	xmax=0.0d0
	ymax=0.0d0
	zmax=0.0d0
	dx=0.0d0
	dy=0.0d0
	dz=0.0d0

	do i=1,natom
          x(i)=dmod(dmod(x(i),boxx)+boxx,boxx)
          y(i)=dmod(dmod(y(i),boxy)+boxy,boxy)
          z(i)=dmod(dmod(z(i),boxz)+boxz,boxz)
	  xmin=min(x(i),xmin)
	  ymin=min(y(i),ymin)
	  zmin=min(z(i),zmin)
	  xmax=max(x(i),xmax)
	  ymax=max(y(i),ymax)
	  zmax=max(z(i),zmax)
	enddo

	xmax=min(boxx,xmax+0.001d0*(boxx-xmax))
	ymax=min(boxy,ymax+0.001d0*(boxy-ymax))
	zmax=min(boxz,zmax+0.001d0*(boxz-zmax))
	
	if (boxx.lt.boxxn) then
	  dx=min(boxxs,(boxxn-boxx))/2.0d0
	  boxx=min(boxx+boxxs,boxxn)
	endif
	if (boxy.lt.boxyn) then
	  dy=min(boxys,(boxyn-boxy))/2.0d0
	  boxy=min(boxy+boxys,boxyn)
	endif
	if (boxz.lt.boxzn) then
	  dz=min(boxzs,(boxzn-boxz))/2.0d0
	  boxz=min(boxz+boxzs,boxzn)
	endif
	
	if (boxx.gt.boxxn) then
	  xmin=min(xmin,boxxs/2.0d0)
	  dx=-xmin
	  xmax=boxx-max(xmax,boxx-boxxs/2.0d0)
	  boxx=max(boxx-xmin-xmax,boxxn)
	endif
	if (boxy.gt.boxyn) then
	  ymin=min(ymin,boxys/2.0d0)
	  dy=-ymin
	  ymax=boxy-max(ymax,boxy-boxys/2.0d0)
	  boxy=max(boxy-ymin-ymax,boxyn)
	endif
	if (boxz.gt.boxzn) then
	  zmin=min(zmin,boxzs/2.0d0)
	  dz=-zmin
	  zmax=boxz-max(zmax,boxz-boxzs/2.0d0)
	  boxz=max(boxz-zmin-zmax,boxzn)
	endif

	do i=1,natom
	  x(i)=x(i)+dx
	  y(i)=y(i)+dy
	  z(i)=z(i)+dz
	enddo

	if (boxx.eq.boxxn .and.
     $	    boxy.eq.boxyn .and.
     $	    boxz.eq.boxzn) then
	  bshr='-'
	endif
	
	write (*,'(A,3(2X,F10.5))') 'Neue Box:',boxx,boxy,boxz

	call stopdead()

        vxs=0.0d0
        vys=0.0d0
        vzs=0.0d0

        do i=1,natom
          vxs=vxs+vx(i)*1d12
          vys=vys+vy(i)*1d12
          vzs=vzs+vz(i)*1d12
        enddo
	
	write (*,'(A,3(2X,F12.10))') 'Systemgeschw.:',vxs,vys,vzs	

	call neighbor(-step)

	return
	end

*-------------------------------------------------------------------------------
