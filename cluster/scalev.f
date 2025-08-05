*===============================================================================
* 	THEROMSTAT BY VELOCITY SCALING
*===============================================================================
* file: scalev.f
*-------------------------------------------------------------------------------
* This routine scales the velocities of the particles to the correct
* kinetic energy, and thus sets the correct temperature. A better
* alternative would be using the Nose-Hoover thermostat, but this
* essentially does the same, just has a more fancy way of arriving at the
* scaling factor.
*-------------------------------------------------------------------------------

	subroutine thermostat(step)

	implicit none
	
	include 'const.inc'
	include 'energ.inc'
	include 'atomp.inc'
	include 'atomc.inc'

	integer step
	integer i
	double precision kf,ekinT

*-------------------------------------------------------------------------------
* step		: time step number
*...............................................................................
* i		: standard loop variable
* kf		: correction factor for velocities
* ekinT		: kinetic energy required for temperature T
*-------------------------------------------------------------------------------

	ekinT=1.5d0*cikb*att*natom
	kf=sqrt(ekinT/ekin)

	do i=1,natom
	  vx(i)=vx(i)*kf
	  vy(i)=vy(i)*kf
	  vz(i)=vz(i)*kf
	enddo
	ekin=ekinT
	
	return
	end
*-------------------------------------------------------------------------------
