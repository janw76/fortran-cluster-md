*===============================================================================
*	ACCELERATION ATOM-ATOM
*===============================================================================
* file: for-a-lj3.f
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
	double precision eij,fij,esum
	double precision dx,dy,dz,tx,ty,tz
	double precision c1,c2,c3
	double precision fsumx(mxa),fsumy(mxa),fsumz(mxa)

*-------------------------------------------------------------------------------
* i,j		: standard loop variable
* dx,dy,dz	: components of distance between two atoms i,j
* tx,ty,tz	: temporary variables for coordinates of atom i
* c1		: constant for the 1st LJ-term
* c2		: constant for the 2nd LJ-term
* c3		: constant for the 1st LJ-term * 1/2
* atr2		; inverse of squared diameter of the atoms
* rij2		: inverse of the distance between two atoms, squared 
* rij6		: inverse of the distance between two atoms, 6th power
* eij		: potential energy between two atoms i,j
* fij		: force between two atoms i,j
* esum		: sum of potential, converted into epot afterwards
* fsumx()	: field for sum of forces, converted into ax() afterwards
* fsumy()	: field for sum of forces, converted into ax() afterwards
* fsumz()	: field for sum of forces, converted into ax() afterwards
*-------------------------------------------------------------------------------

	atr2=1.0d0/(25.0d0*atr*atr)
	c1=atr*atr*atr*atr*atr*atr
	c2=c1*c1
	c3=0.5d0*c1
	esum=0.0d0
	do i=1,natom
	  fsumx(i)=0.0d0
	  fsumy(i)=0.0d0
	  fsumz(i)=0.0d0
	enddo

	do i=1,natom
	  tx=x(i)
	  ty=y(i)
  	  tz=z(i)
	  do j=i+1,natom
	    dx=tx-x(j)
	    dy=ty-y(j)
	    dz=tz-z(j)
	    rij2=1.0d0/(dx*dx+dy*dy+dz*dz)
*	    if (rij2.gt.atr2) then
*	      write (*,*) 'Eeeeeeeeeek.....'
*	      write (*,*) i,j,1/sqrt(rij2),atr
*	      write (*,'(4F8.4)') x(1),x(2),y(1),y(2)
*	      write (*,'(3F8.4)') dx,dy,dz
*	    endif
	    rij6=rij2*rij2*rij2
	    eij=rij6*(c2*rij6-c1)
	    fij=rij2*rij6*(c2*rij6-c3)
	    fsumx(i)=fsumx(i)+fij*dx
	    fsumx(j)=fsumx(j)-fij*dx
	    fsumy(i)=fsumy(i)+fij*dy
	    fsumy(j)=fsumy(j)-fij*dy
	    fsumz(i)=fsumz(i)+fij*dz
	    fsumz(j)=fsumz(j)-fij*dz
	    esum=esum+eij
	  enddo
	enddo

	epot=esum*ate*4.0d0
	do i=1,natom
	  ax(i)=ax(i)+ate*48.0d0*fsumx(i)/atm
	  ay(i)=ay(i)+ate*48.0d0*fsumy(i)/atm
	  az(i)=az(i)+ate*48.0d0*fsumz(i)/atm
	enddo
	return
	end
*-------------------------------------------------------------------------------
