*===============================================================================
* 	INTEGRATOR EULER 2
*===============================================================================
* file: int-euler2.f
*-------------------------------------------------------------------------------
* This routine integrates the trajectories by using the improved Euler
* algorithm:
* 	x(t+dt) = x(t) + (v(t) + a(t)dt)dt
*	v(t+dt) = v(t) + a(t)dt
* This routine is for testing purposes only, it is still way too inaccurate 
* for real simulations.
*-------------------------------------------------------------------------------

	subroutine integrator(step)
	
	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'energ.inc'

	integer i,step
	double precision kf

*-------------------------------------------------------------------------------
* step		: Time step
*...............................................................................
* i		: standard loop variable
* kf		: constant for acceleration
*-------------------------------------------------------------------------------

	do i=1,natom
	  ax(i)=0
	  ay(i)=0
	  az(i)=0
	enddo

	epot=0
	ekin=0
	kf=0.5
		
	call AccAtom()
	call AccWall()
	
	do i=1,natom
	  ekin=ekin+0.5d0*atm*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
	  x(i)=x(i)+(vx(i)+kf*ax(i)*dt)*dt
	  vx(i)=vx(i)+ax(i)*dt
	  y(i)=y(i)+(vy(i)+kf*ay(i)*dt)*dt
	  vy(i)=vy(i)+ay(i)*dt
	  z(i)=z(i)+(vz(i)+kf*az(i)*dt)*dt
	  vz(i)=vz(i)+az(i)*dt
*	  ekin=ekin+0.5d0*atm*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
	enddo

	return 1
	end
*-------------------------------------------------------------------------------
