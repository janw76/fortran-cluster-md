*===============================================================================
* 	CHANGE NATOM
*===============================================================================
* file: chg-natom.f
* date: 18-jan-2001
*-------------------------------------------------------------------------------
* This routine adjusts the total particle number by +/- 1 by adding a
* single atom at a random position (checking for overlap!), or removing
* a single monomer (!).
* Removal takes place only for monomers, while the minimum distance for an
* inserted atom from everthing else is the Stillinger-Radius.
* Neighborlists, Temperature, total momentum etc are adjusted accordingly.
*-------------------------------------------------------------------------------

	subroutine chgNatom(step)

	implicit none

	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'cvtpp.inc'
	include 'boxpp.inc'

	integer step
	integer i,ntry
	integer nm,nc,xa
	double precision xx,yy,zz,dx,dy,dz,rr
	integer kpsMAX
	double precision kpsRND

*-------------------------------------------------------------------------------
* step		: Time step
*...............................................................................
* i		: standard loop variable
* ntry		: number of tries for atom insertion
* nm		: number of monemer to be destroyed
* nc		: cluster number of monomer nm
* xa		: atom number of monomer nm
* xx,yy,zz	: test coordinates for inserted atom
* dx,dy,dz	: components of distance between two atoms
* rr		: minimum distance to next neigbor of inserted atom
*-------------------------------------------------------------------------------

	if (nnatom.lt.natom) then

	  if (step.ne.cstep) call stoddard(step)
	  if (csi(1).eq.0) then
	    write(*,*) '*** INFO *** No monomers to remove! Skipping.'
	    return
	  endif
	  nm=kpsMAX(csi(1))+1
	  do i=1,natom
	    if (cs(i).eq.1) then
	      nm=nm-1
	      if (nm.eq.0) nc=i
	    endif
	  enddo
	  do i=1,natom
	    if (cl(i).eq.nc) xa=i
	  enddo
	  write (*,*) 'Removing atom #',xa,'... Natom now:',natom-1
	  x0(xa)=x0(natom)
	  y0(xa)=y0(natom)
	  z0(xa)=z0(natom)
	  x1(xa)=x1(natom)
	  y1(xa)=y1(natom)
	  z1(xa)=z1(natom)
	  x2(xa)=x2(natom)
	  y2(xa)=y2(natom)
	  z2(xa)=z2(natom)
	  x3(xa)=x3(natom)
	  y3(xa)=y3(natom)
	  z3(xa)=z3(natom)
	  x4(xa)=x4(natom)
	  y4(xa)=y4(natom)
	  z4(xa)=z4(natom)
	  x5(xa)=x5(natom)
	  y5(xa)=y5(natom)
	  z5(xa)=z5(natom)
	  natom=natom-1

	else

	  ntry=0
10	  ntry=ntry+1
	  if (ntry.gt.10) then
	    write(*,*) '*** INFO *** No free space found after 10 tries'
	    return
	  endif
	  xx=boxx*kpsRND()
	  yy=boxy*kpsRND()
	  zz=boxz*kpsRND()
	  rr=boxx*boxx+boxy*boxy+boxz*boxz
	  do i=1,natom
	    dx=x0(i)-xx
	    dy=y0(i)-yy
	    dz=z0(i)-zz
	    rr=min(rr,dx*dx+dy*dy+dz*dz)
	  enddo
	  if (dsqrt(rr).le.atskin) goto 10
	  i=kpsMAX(natom)+1
	  natom=natom+1
	write(*,*) 'Inserting after',ntry,' tries ... Natom now: ',natom
	  x0(natom)=xx
	  y0(natom)=yy
	  z0(natom)=zz
	  x1(natom)=x1(i)
	  y1(natom)=y1(i)
	  z1(natom)=z1(i)
	  x2(natom)=x2(i)
	  y2(natom)=y2(i)
	  z2(natom)=z2(i)
	  x3(natom)=x3(i)
	  y3(natom)=y3(i)
	  z3(natom)=z3(i)
	  x4(natom)=x4(i)
	  y4(natom)=y4(i)
	  z4(natom)=z4(i)
	  x5(natom)=x5(i)
	  y5(natom)=y5(i)
	  z5(natom)=z5(i)
	endif

	call stopdead(step)
	call neighbor(-step)

	if (nnatom.eq.natom) then
	  chgNctl='-'
	endif

	return
	end

*-------------------------------------------------------------------------------
