*===============================================================================
*	 INTIALIZE STARTING CONFIGURATION
*===============================================================================
* file: init-ball.f
*-------------------------------------------------------------------------------
* This routine set an any desired number of atoms in a layered sphere.
*-------------------------------------------------------------------------------

	subroutine initconf()
	
	implicit none

	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'boxpp.inc'
	include 'energ.inc'

	integer i,j
	integer atomnr,nring,layer,nlayer,nphi
	double precision urx,ury,urz
	double precision phi,phi0,phiD
	double precision theta,theta0,thetaD
	double precision kf,rad
	double precision kpsRND

*-------------------------------------------------------------------------------
* i,j		: standard loop variables
* atomnr	: index of the target atom
* nring		: number of atoms in target ring
* layer		: index of spherelayer
* nlayer	: number of atoms in target layer
* nphi 		: number of rings in a layer
* urx,ury,urz	: translation vector to box center
* phi		: elevation angle for ringlevel
* phi0		: elevation angle start
* phiD		: elevation angle step
* theta		: azimuth
* theta0	: azimuth starting angle (random)
* thetaD	: azimuth step
* rad		: radius of layer/level
* kf		: correction factor for atom numbers
*-------------------------------------------------------------------------------

	if (rndctrl.eq.'X') then
	  call kpsinit(0)
	  write(*,*) 'Using random velocities'
	else
	  call kpsinit(1)
	  write(*,*) 'Using static velocities'
	endif

	urx=boxx/2.0d0
	ury=boxy/2.0d0
	urz=boxz/2.0d0

	x(1)=0.0d0
	y(1)=0.0d0
	z(1)=0.0d0
	vx(1)=(kpsRND()-0.5d0)
	vy(1)=(kpsRND()-0.5d0)
	vz(1)=(kpsRND()-0.5d0)

	layer=1
	phi=0.0d0
	atomnr=2

10	continue
	phiD=acos(1.0d0-(1.0d0/(2.0d0*layer*layer)))
	nphi=int(pi/phiD)
	phi0=(pi-real(nphi)*phiD)/2.0d0
	nlayer=0

        do i=0,nphi
          phi=phi0+real(i)*phiD
          rad=abs(real(layer)*atskin*sin(phi))
          if (phi.lt.1d-10.or.pi-phi.lt.1d-10) then
            nring=1
          else
            nring=int(real(rad)*2.0d0*pi)/atskin
          endif
          nlayer=nlayer+nring
        enddo
	
	if (natom-atomnr.gt.nlayer) then
	  kf=1.0d0
	else
	  kf=1.05d0*real(natom-atomnr+1)/real(nlayer)
	endif
	
	do i=0,nphi
	  phi=phi0+real(i)*phiD
	  rad=abs(real(layer)*atskin*sin(phi))
	  if (phi.lt.1d-10.or.pi-phi.lt.1d-10) then
	    nring=1
	  else
	    nring=int(real(rad)*2.0d0*pi)/atskin
	  endif
	  nring=int(kf*real(min0(natom-atomnr+1,nring))+0.5d0)
	  if (nring.eq.0) nring=1

	  thetaD=2.0d0*pi/real(nring)
	  theta0=kpsRND()*thetaD
	  
	  rad=real(layer)*atskin
          do j=0,nring-1
	    theta=theta0+real(j)*thetaD
            x(atomnr)=rad*sin(phi)*cos(theta)+urx
            y(atomnr)=rad*sin(phi)*sin(theta)+ury
            z(atomnr)=rad*cos(phi)           +urz
	    vx(atomnr)=(kpsRND()-0.5d0)
	    vy(atomnr)=(kpsRND()-0.5d0)
	    vz(atomnr)=(kpsRND()-0.5d0)
            atomnr=atomnr+1
          enddo

	  if (atomnr.gt.natom) goto 30
	enddo
	layer=layer+1
        if (atomnr.le.natom) goto 10
		
30	continue
	
	call stopdead(-1)
	call neighbor(0)
	call AccAtom()
	call AccWall()

	return
	end
*-------------------------------------------------------------------------------
