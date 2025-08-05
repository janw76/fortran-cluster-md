*===============================================================================
*	STOP OVERALL SYSTEM MOVEMENT
*===============================================================================
* file: stopdead-m.f
*-------------------------------------------------------------------------------
* This routine sets the system movement vector to zero. 
* Changed it to run several times to get a better "reset" (JW, 4/10/07)
*-------------------------------------------------------------------------------

	subroutine stopdead(step)
	
	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'energ.inc'

	integer i,step,loop
	double precision vxs,vys,vzs,mas,imas
	
*-------------------------------------------------------------------------------
* i		: standard loop variable
* vxs,vys,vzs	: components of total system velocity
*-------------------------------------------------------------------------------

*** Doing the same stuff 10 times (JW)
	DO LOOP=1,10

	vxs=0.0d0
	vys=0.0d0
	vzs=0.0d0

	do i=1,natom
	  mas=apmas(api(i))
	  vxs=vxs+vx(i)*mas
	  vys=vys+vy(i)*mas
	  vzs=vzs+vz(i)*mas
	enddo

	vxs=vxs/dble(natom)
	vys=vys/dble(natom)
	vzs=vzs/dble(natom)

*	if (step.eq.-1.or.mod(step,1000).eq.-42) then
	ekin=0.0d0
	do i=1,natom
	  mas=apmas(api(i))
	  imas=1.0d0/mas
	  vx(i)=vx(i)-vxs*imas
	  vy(i)=vy(i)-vys*imas
	  vz(i)=vz(i)-vzs*imas
	  ekin=ekin+mas*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
	enddo
	ekin=0.5d0*ekin

	call thermostat(step)

	vxs=0.0d0
	vys=0.0d0
	vzs=0.0d0
	do i=1,natom
	  mas=apmas(api(i))
	  vxs=vxs+vx(i)*mas
	  vys=vys+vy(i)*mas
	  vzs=vzs+vz(i)*mas
	enddo
	
	ENDDO
	
	write (*,*) 'called stopdead !'
	write (*,10) step,vxs/dble(natom),vys/dble(natom),vzs/dble(natom)

10	format (I15,3(2X,F20.18))
	
	return

	end

*-------------------------------------------------------------------------------
