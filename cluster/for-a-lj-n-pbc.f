*===============================================================================
*	ACCELERATION ATOM-ATOM
*===============================================================================
* file: for-a-lj-n-pbc.f
*-------------------------------------------------------------------------------
* In this routine the forces between the individual atoms are summed and
* converted to accelerations. The Lennard-Jones 12-6-potential is used.
* This routine uses the LJ-potential with a cutoff distance.
* Additionally, the routine uses a neighborlist to cut down the numeber of
* time consuming force calculations.
* Furthermore, periodic boundary conditions are in use.
*-------------------------------------------------------------------------------

	subroutine AccAtom()
	
	implicit none

	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'energ.inc'
	include 'boxpp.inc'

	integer i,j,k
	double precision atr2,rij2,rij6
	double precision fij,esum
	double precision dx,dy,dz,tx,ty,tz
	double precision c1,c2,c3,bf1,bf2
	double precision fsumx(mxa),fsumy(mxa),fsumz(mxa)
	double precision ibx,iby,ibz,iatm
	double precision virial

*-------------------------------------------------------------------------------
* i,j		: standard loop variables
* k		: index for neighborlist
* dx,dy,dz	: components of distance between two atoms i,j
* tx,ty,tz	: temporary variables for coordinate of atom i
* c1		: constant for the 1st LJ-term
* c2		: constant for the 2nd LJ-term
* c3		: constant for the 1st LJ-term * 1/2
* atr2		; inverse of squared diameter of the atoms
* rij2		: inverse of the distance between two atoms, squared 
* rij6		: inverse of the distance between two atoms, 6th power
* fij		: force between two atoms i,j
* esum		: sum of potential, converted to epot afterwards
* fsumx()	: field for sum of forces, converted into ax() afterwards
* fsumy()	: field for sum of forces, converted into ax() afterwards
* fsumz()	: field for sum of forces, converted into ax() afterwards
* bf1,bf2	: box-factors for minimum image convention
* ibx,iby,ibz	: inverse of Boxlengths (needed for PBC)
* iatm		: inverse of the atomic mass (needed for force-->vel)
*-------------------------------------------------------------------------------

	atr2=atrcut*atrcut
	c1=atr*atr*atr*atr*atr*atr
	c2=c1*c1
	c3=0.5d0*c1
	bf1=2.0d0+int(atrskin/atr)
	bf2=bf1+0.5d0
	esum=0.0d0

	ibx=1.0d0/boxx
	iby=1.0d0/boxy
	ibz=1.0d0/boxz

	do i=1,natom
	  fsumx(i)=0.0d0
	  fsumy(i)=0.0d0
	  fsumz(i)=0.0d0
	enddo

	k=1

10	i=-atnlist(k)
	if (i.eq.natom) goto 30

	  k=k+1
	  tx=x(i)
	  ty=y(i)
  	  tz=z(i)

20	  j=atnlist(k)
	  if (j.gt.0) then
	    k=k+1

	    dx=tx-x(j)
	    dy=ty-y(j)
	    dz=tz-z(j)
	    dx=dx+boxx*(bf1-int(bf2+dx*ibx))
	    dy=dy+boxy*(bf1-int(bf2+dy*iby))
	    dz=dz+boxz*(bf1-int(bf2+dz*ibz))

	    rij2=dx*dx+dy*dy+dz*dz

	    if (rij2.le.atr2) then
	      rij2=1.0d0/rij2
	      rij6=rij2*rij2*rij2
	      fij=rij2*rij6*(c2*rij6-c3)
	      fsumx(i)=fsumx(i)+fij*dx
	      fsumx(j)=fsumx(j)-fij*dx
	      fsumy(i)=fsumy(i)+fij*dy
	      fsumy(j)=fsumy(j)-fij*dy
	      fsumz(i)=fsumz(i)+fij*dz
	      fsumz(j)=fsumz(j)-fij*dz
	      esum=esum+rij6*(c2*rij6-c1)
*	      virial=virial+rij6*(c2*rij6-c3)
	    endif
	    goto 20
	  endif
	goto 10

30	continue

*	virial=virial*48.0d0*ate
*	write (44,*) virial
*	virial= (1.0d13/(boxx*boxy*boxz))*
*     $          (dble(natom)*cikb*att-virial/3.0d0)

*	write (44,*) virial

	iatm=1.0d0/atm
	epot=esum*ate*4.0d0
	do i=1,natom
	  ax(i)=ax(i)+ate*48.0d0*fsumx(i)*iatm
	  ay(i)=ay(i)+ate*48.0d0*fsumy(i)*iatm
	  az(i)=az(i)+ate*48.0d0*fsumz(i)*iatm
	enddo
	return
	end
*-------------------------------------------------------------------------------
