*===============================================================================
* 	NEIGHBOR LIST
*===============================================================================
* file: neighbor.f
*-------------------------------------------------------------------------------
* This routine determines the neighbors for each atom
* this saves a lot of (nearly) useless force calculations
*-------------------------------------------------------------------------------

	subroutine neighbor(step)

	implicit none

	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'timep.inc'
	include 'dmmpp.inc'

	integer step
	integer i,j,ind,istep
	double precision tx,ty,tz,dx,dy,dz,r2,rs2

*-------------------------------------------------------------------------------
* step		: timestep number
*...............................................................................
* istep		: local |step|; never write to the dummy arg -- callers pass
*		  literals (init-cube1.f: call neighbor(0)), which live in
*		  read-only memory (SIGBUS on write)
* i		: standard loop variable
* j		: standard loop variable
* ind		: index for neighborlist
* tx,ty,tz	: temporary variables for position of atom i
* dx,dy,dz	: components for distance of two atoms
* r2		: Square of distance between two atoms
* rs2		: square of skin radius
*-------------------------------------------------------------------------------

	ind=1
	rs2=atrskin*atrskin

	istep=abs(step)

*** no elapsed steps yet -> no average (guard against /0, cf. neighbor-pbc-mav.f)
	if (istep.gt.btim) then
	  dmmav=dmmsu*dble(nbtim)/dble(istep-btim)
	else
	  dmmav=0.0d0
	endif

	if (dmms.gt.(atrskin-atrcut)) then
	  write (*,*) '*** FATAL ERROR! ***  Skin radius exceeded!!'
	  write (*,*) 'Maximum:',atrskin-atrcut
	  write (*,*) 'Sum    :',dmms
*	  call wrtXYZ(step,42)
*	  call wrtVEL(step,43)
*	  call wrtACC(step,44)
	  stop 1
	endif

	if (dmmst.eq.nbtim.and.
     $	dmms.lt.dmmav.and.dmmav.lt.0.95d0*(atrskin-atrcut)) then
	  write (*,*) '*** INFO *** Reducing Skin Radius by 1%'
	  atrskin=atrskin-0.01d0*(atrskin-atrcut)
	endif

	do i=1,natom
	  atnlist(ind)=-i
	  atnidx(i)=ind
	  ind=ind+1
	  tx=x(i)
	  ty=y(i)
	  tz=z(i)
	  do j=i+1,natom
	    dx=tx-x(j)
	    dy=ty-y(j)
            dz=tz-z(j)
            r2=dx*dx+dy*dy+dz*dz
	    if (r2.le.rs2) then
	      atnlist(ind)=j
	      ind=ind+1
	    endif
	  enddo
	  if (ind+natom.gt.mxnlist) then
	    write (*,*) '*** FATAL ERROR! ***  '//
     $		        'Neighborlist is nearly full... stopping'
	    stop 1
	  endif
	enddo

	dmmax=max(dmms,dmmax)
	write (*,10) istep,ind-1-natom,dmms,dmmav,dmmax,
     $	             atrskin-atrcut, atrskin/atr
	dmms=0.0d0
	dmmst=0
10	format (2(I10,2X),5(F8.4,2X))

	return
	end

*-------------------------------------------------------------------------------
