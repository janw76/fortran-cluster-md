*===============================================================================
* 	INTEGRATOR VERLET
*===============================================================================
* file: int-verlet-dmm.f
*-------------------------------------------------------------------------------
* This routine integrates the equations of motion with the Verlet-method,
* modified by Swope, Anderson, Berens and Wilson.
* 	x(t+   dt) = x(t) + (v(t) + 0.5 a(t)dt)dt
*	v(t+0.5dt) = v(t) + 0.5 a(t)dt
*	calculate a(t+dt)
*	v(t+   dt) = v(t+0.5dt) + 0.5 a(t+dt)
* This routine is commonly used in simulations and is universally regarded
* as being fairly accurate.
* This Routine additionally calculates the sum over the maximum
* displacements for use in the neighborlist
*-------------------------------------------------------------------------------

	subroutine integrator(step)
	
	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'energ.inc'
	include 'dmmpp.inc'

	integer i,step
	double precision dx,dy,dz

*-------------------------------------------------------------------------------
* step		: time step
*...............................................................................
* i		: standard loop variable
* dx,dy,dz	: components of the displacement
*-------------------------------------------------------------------------------

	dmmx=0.0d0
	dmmy=0.0d0
	dmmz=0.0d0
	do i=1,natom
	  dx=(vx(i)+0.5d0*ax(i)*dt)*dt
	  dy=(vy(i)+0.5d0*ay(i)*dt)*dt
	  dz=(vz(i)+0.5d0*az(i)*dt)*dt
	  x(i)=x(i)+dx
	  y(i)=y(i)+dy
	  z(i)=z(i)+dz
	  vx(i)=vx(i)+0.5*ax(i)*dt
	  vy(i)=vy(i)+0.5*ay(i)*dt
	  vz(i)=vz(i)+0.5*az(i)*dt
	  dmmx=max(dmmx,dx*dx)
	  dmmy=max(dmmy,dy*dy)
	  dmmz=max(dmmz,dz*dz)
	enddo

	dmms=dmms+2.0d0*dsqrt(dmmx+dmmy+dmmz) 
	dmmsu=dmmsu+2.0d0*dsqrt(dmmx+dmmy+dmmz)
	dmmst=dmmst+1

	if (mod(step,nbtim).eq.0) call neighbor(step)

	if (dmms.gt.0.95d0*(atrskin-atrcut)) then
	  write (*,*) '*** WARNING! *** Skin Radius nearly exceeded!'//
     $                ' Increasing by 1%'
	  atrskin=atrskin+(atrskin-atrcut)*0.01d0
	  call neighbor(step)
	endif

	epot=0.0d0
	do i=1,natom
	  ax(i)=0.0d0
	  ay(i)=0.0d0
	  az(i)=0.0d0
	enddo
	
	call AccAtom()
	call AccWall()
	
	ekin=0.0d0
	do i=1,natom
	  vx(i)=vx(i)+0.5*ax(i)*dt
	  vy(i)=vy(i)+0.5*ay(i)*dt
	  vz(i)=vz(i)+0.5*az(i)*dt
	  ekin=ekin+0.5d0*atm*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
	enddo

	return
	end
*-------------------------------------------------------------------------------
