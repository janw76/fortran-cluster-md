*===============================================================================
*	ACCELERATION ATOM-ATOM
*===============================================================================
* file: for-a-lj1.f
*-------------------------------------------------------------------------------
* In this routine the forces between the individual atoms are summed and
* converted to accelerations. The Lennard-Jones 12-6-potential is used.
* This routine uses the *full* LJ-potential and is only usable for
* finite systems with box walls.
*-------------------------------------------------------------------------------

	subroutine AccAtom()
	
	implicit none

	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'energ.inc'

	integer i,j
	double precision atr2,rij2,rij6
	double precision eij,fij
	double precision dx,dy,dz
	double precision c1,c2,c3

*-------------------------------------------------------------------------------
* i,j		: standard loop variables
* dx,dy,dz	: components of distance between two atoms i,j
* c1		: constant for the 1st LJ-term
* c2		: constant for the 2nd LJ-term
* c3		: constant for the 1st LJ-term * 1/2
* atr2		; inverse of squared diameter of the atoms
* rij2		: inverse of the distance between two atoms, squared 
* rij6		: inverse of the distance between two atoms, 6th power
* eij		: potential energy between two atoms i,j
* fij		: force between two atoms i,j
*-------------------------------------------------------------------------------

	atr2=1/(0.2*0.2*atr*atr)
	c1=atr*atr*atr*atr*atr*atr
	c2=c1*c1
	c3=0.5d0*c1
	    	
	do i=1,natom
	  do j=i+1,natom
	    dx=x(i)-x(j)
	    dy=y(i)-y(j)
	    dz=z(i)-z(j)
	    rij2=1.0d0/(dx*dx+dy*dy+dz*dz)
	    if (rij2.gt.atr2) then
	      write (*,*) 'Eeeeeeeeeek.....'
	      write (*,*) i,j,1/sqrt(rij2),atr
*	      write (*,'(4F8.4)') x(1),x(2),y(1),y(2)
	      write (*,'(3F8.4)') dx,dy,dz
	    endif
	    rij6=rij2*rij2*rij2
	    eij=ate*4d0*rij6*(c2*rij6-c1)
	    fij=ate*rij2*48d0*rij6*(c2*rij6-c3)
	    ax(i)=ax(i)+fij*dx/atm
	    ax(j)=ax(j)-fij*dx/atm
	    ay(i)=ay(i)+fij*dy/atm
	    ay(j)=ay(j)-fij*dy/atm
	    az(i)=az(i)+fij*dz/atm
	    az(j)=az(j)-fij*dz/atm
	    epot=epot+eij
	  enddo
	enddo
	return
	end
*-------------------------------------------------------------------------------
