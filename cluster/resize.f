*===============================================================================
* 	RESIZE BOX
*===============================================================================
* file: resize.f
*-------------------------------------------------------------------------------
* This routine expands or shrinks the box according to parameters
* still not very optimized! (but functional, at least)
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
*-------------------------------------------------------------------------------

	xmin=boxx
	ymin=boxy
	zmin=boxz
	xmax=0
	ymax=0
	zmax=0

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

	if (boxx.lt.boxxn) then
	  do i=1,natom
	    x(i)=x(i)+boxxs/2.0d0
	  enddo
	  boxx=min(boxx+boxxs,boxxn)
	endif

	if (boxy.lt.boxyn) then
	  do i=1,natom
	    y(i)=y(i)+boxys/2.0d0
	  enddo
	  boxy=min(boxy+boxys,boxyn)
	endif

	if (boxz.lt.boxzn) then
	  do i=1,natom
	    z(i)=z(i)+boxzs/2.0d0
	  enddo
	  boxz=min(boxz+boxzs,boxzn)
	endif

	if (boxx.gt.boxxn) then
	  xmin=min(xmin,boxxs/2.0d0)
	  xmax=boxx-max(xmax,boxx-boxxs/2.0d0)
	  do i=1,natom
	    x(i)=x(i)-xmin
	  enddo
	  boxx=max(boxx-xmin-xmax,boxxn)
	endif

	if (boxy.gt.boxyn) then
	  ymin=min(ymin,boxys/2.0d0)
	  ymax=boxy-max(ymax,boxy-boxys/2.0d0)
	  do i=1,natom
	    y(i)=y(i)-ymin
	  enddo
	  boxy=max(boxy-ymin-ymax,boxyn)
	endif

	if (boxz.gt.boxzn) then
	  zmin=min(zmin,boxzs/2.0d0)
	  zmax=boxz-max(zmax,boxz-boxzs/2.0d0)
	  do i=1,natom
	    z(i)=z(i)-zmin
	  enddo
	  boxz=max(boxz-zmin-zmax,boxzn)
	endif

	if (boxx.eq.boxxn .and.
     $	    boxy.eq.boxyn .and.
     $	    boxz.eq.boxzn) then
	  bshr='-'
	endif
	
	write (*,'(A,3(2X,F10.5))') 'Neue Box:',boxx,boxy,boxz
	call neighbor(-step)

	return
	end

*-------------------------------------------------------------------------------
