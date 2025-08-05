*===============================================================================
*	STOP OVERALL SYSTEM MOVEMENT
*===============================================================================
* file: stopdead.f
*-------------------------------------------------------------------------------
* This routine sets the system movement vector to zero. 
*-------------------------------------------------------------------------------

	subroutine stopdead(step)
	
	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'energ.inc'

	integer i,step
	double precision vxs,vys,vzs
	
*-------------------------------------------------------------------------------
* i		: standard loop variable
* vxs,vys,vzs	: components of total system velocity
*-------------------------------------------------------------------------------

	vxs=0.0d0
	vys=0.0d0
	vzs=0.0d0

	do i=1,natom
	  vxs=vxs+vx(i)
	  vys=vys+vy(i)
	  vzs=vzs+vz(i)
	enddo

	vxs=vxs/dble(natom)
	vys=vys/dble(natom)
	vzs=vzs/dble(natom)

	ekin=0.0d0
	do i=1,natom
	  vx(i)=vx(i)-vxs
	  vy(i)=vy(i)-vys
	  vz(i)=vz(i)-vzs
	  ekin=ekin+vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
	enddo
	ekin=0.5d0*atm*ekin

	call thermostat(step)

*10	format (I15,3(2X,F17.14))

	return
	end

*-------------------------------------------------------------------------------
