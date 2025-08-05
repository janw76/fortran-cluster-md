*===============================================================================
* 	INTEGRATOR EULER 1
*===============================================================================
* file: int-euler1.f
*-------------------------------------------------------------------------------
* This routine integrates the trajectories by using the simple Euler
* algorithm:
* 	x(t+dt) = x(t) + v(t)dt
*	v(t+dt) = v(t) + a(t)dt
* This routine is for testing purposes only, it is way too inaccurate for
* real simulations.
*-------------------------------------------------------------------------------

	subroutine integrator(step)
	
	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'energ.inc'

	integer step,i

*-------------------------------------------------------------------------------
* step		: time step
*...............................................................................
* i		: standard loop variable
*-------------------------------------------------------------------------------

	do i=1,natom
	  ax(i)=0
	  ay(i)=0
	  az(i)=0
	enddo

	epot=0
	ekin=0
		
	call AccAtom()
	call AccWall()
	
	do i=1,natom
	  ekin=ekin+.5*atm*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
	  x(i)=x(i)+vx(i)*dt
	  vx(i)=vx(i)+ax(i)*dt
	  y(i)=y(i)+vy(i)*dt
	  vy(i)=vy(i)+ay(i)*dt
	  z(i)=z(i)+vz(i)*dt
	  vz(i)=vz(i)+az(i)*dt
	enddo

	return
	end
*-------------------------------------------------------------------------------
