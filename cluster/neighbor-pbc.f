*===============================================================================
* 	NEIGHBOR LIST - PBC
*===============================================================================
* file: neighbor-pbc.f
*-------------------------------------------------------------------------------
* This routine determines the neighbors for each atom
* this saves a lot of (nearly) useless force calculations
* This routine additionally uses the periodic boundary conditions
*-------------------------------------------------------------------------------

	subroutine neighbor(step)

	implicit none

	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'boxpp.inc'
	include 'timep.inc'
	include 'dmmpp.inc'

	integer step
	integer i,j,ind
	double precision tx,ty,tz,dx,dy,dz,r2,rs2
	double precision ibx,iby,ibz

*-------------------------------------------------------------------------------
* step		: timestep number
*...............................................................................
* i		: standard loop variable
* j		: standard loop variable
* ind		: index for neighborlist
* tx,ty,tz	: temporary variables for position of atom i
* dx,dy,dz	: components for distance of two atoms
* r2		: Square of distance between two atoms
* rs2		: square of skin radius
* ibx,iby,ibz	: inverse of boxlengths (needed for PBC)
*-------------------------------------------------------------------------------

	ind=1
*	step=abs(step)
	if (step.lt.0) step=-step
	
	ibx=1.0d0/boxx
	iby=1.0d0/boxy
	ibz=1.0d0/boxz

	dmmav=dmmsu*dble(nbtim)/dble(step-btim)

	if (dmms.gt.(atrskin-atrcut)) then
	  write (*,*) '*** FATAL ERROR! ***  Skin radius exceeded!!'
	  stop
	endif

	if (dmmst.eq.nbtim.and.
     $	dmms.lt.dmmav.and.dmmav.lt.0.95d0*(atrskin-atrcut)) then
	  write (*,*) '*** INFO *** Reducing Skin Radius by 1%'
	  atrskin=atrskin-0.01d0*(atrskin-atrcut)
	endif
	rs2=atrskin*atrskin

	do i=1,natom
	  x(i)=dmod(dmod(x(i),boxx)+boxx,boxx)
	  y(i)=dmod(dmod(y(i),boxy)+boxy,boxy)
	  z(i)=dmod(dmod(z(i),boxz)+boxz,boxz)
	enddo

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
	    dx=dx+boxx-boxx*int(1.5d0+dx*ibx)
	    dy=dy+boxy-boxy*int(1.5d0+dy*iby)
	    dz=dz+boxz-boxz*int(1.5d0+dz*ibz)

	    r2=dx*dx+dy*dy+dz*dz
	    if (r2.le.rs2) then
	      atnlist(ind)=j
	      ind=ind+1
	    endif
	  enddo
	  if (ind+natom.gt.mxnlist) then
	    write (*,*) '*** FATAL ERROR *** '//
     $                  'Neighborlist is nearly full!'
	    stop
	  endif
	enddo

	dmmax=max(dmms,dmmax)
	write (*,10) step,ind-1-natom,dmms,dmmav,dmmax,
     $	             atrskin-atrcut, atrskin/atr
	dmms=0.0d0
	dmmst=0
10	format (2(I10,2X),5(F8.4,2X))

	return
	end

*-------------------------------------------------------------------------------
