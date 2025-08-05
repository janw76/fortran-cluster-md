*===============================================================================
* 	INTEGRATOR VERLET
*===============================================================================
* file: int-verlet.f
*-------------------------------------------------------------------------------
* This routine integrates the equations of motion with the Verlet-method,
* modified by Swope, Anderson, Berens and Wilson.
* 	x(t+   dt) = x(t) + (v(t) + 0.5 a(t)dt)dt
*	v(t+0.5dt) = v(t) + 0.5 a(t)dt
*	calculate a(t+dt)
*	v(t+   dt) = v(t+0.5dt) + 0.5 a(t+dt)
* This routine is commonly used in simulations and is universally regarded
* as being fairly accurate.
*-------------------------------------------------------------------------------

	subroutine integrator(step)
	
	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'energ.inc'

	integer i,step

*-------------------------------------------------------------------------------
* step		: time step
*...............................................................................
* i		: standard loop variable
*-------------------------------------------------------------------------------

	do i=1,natom
	  x(i)=x(i)+(vx(i)+0.5d0*ax(i)*dt)*dt
	  y(i)=y(i)+(vy(i)+0.5d0*ay(i)*dt)*dt
	  z(i)=z(i)+(vz(i)+0.5d0*az(i)*dt)*dt
	  vx(i)=vx(i)+0.5*ax(i)*dt
	  vy(i)=vy(i)+0.5*ay(i)*dt
	  vz(i)=vz(i)+0.5*az(i)*dt
	enddo

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
