*===============================================================================
* 	THEROMSTAT BY VELOCITY SCALING
*===============================================================================
* file: scalev-clu.f
*-------------------------------------------------------------------------------
* This routine scales the velocities of the particles to the correct
* kinetic energy, and thus sets the correct temperature. A better
* alternative would be using the Nose-Hoover thermostat, but this
* essentially does the same, just has a more fancy way of arriving at the
* scaling factor.
* This routine additionally does a seperate thermostating of the biggest
* cluster and the rest of the system in order to faster arrive in an
* equilibrium
*-------------------------------------------------------------------------------

	subroutine thermostat(step)

	implicit none
	
	include 'const.inc'
	include 'energ.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'cvtpp.inc'

	integer step
	integer i,ci,cb
	double precision kfC,kfG,ekinG,ekinT

*-------------------------------------------------------------------------------
* step		: Time step
*...............................................................................
* i		: standard loop variable
* ci		: index of biggest cluster
* kfG		: correction factor for GAS-phase
* kfC		: correction factor for CLUSTER
* ekinG		: Kinetic energy of gas phase
* ekinT		: required kinetic energy for temperature T
*-------------------------------------------------------------------------------

	ekinT=1.5d0*cikb*att*natom

	if (eclth.eq.'X'.and.mod(step,ecltim).eq.0) then

	  if (cstep.ne.step) call stoddard(step)
	  cb=0
	  do i=1,natom
	    if (cb.lt.cs(i)) then
	      cb=cs(i)
	      ci=i
	    endif
	  enddo

	  ekinG=0.0d0
 	  do i=1,natom
	    if (i.ne.ci) ekinG=ekinG+cek(i)
	  enddo

	  kfG=sqrt(1.5d0*cikb*att*(natom-cs(ci))/ekinG)
	  kfC=sqrt(1.5d0*cikb*att*cs(ci)/cek(ci))
	  do i=1,natom
	    if (cl(i).eq.ci) then
	      vx(i)=vx(i)*kfC
	      vy(i)=vy(i)*kfC
	      vz(i)=vz(i)*kfC
	    else
	      vx(i)=vx(i)*kfG
	      vy(i)=vy(i)*kfG
	      vz(i)=vz(i)*kfG
	    endif
	  enddo

	else

	  kfG=sqrt(ekinT/ekin)
	  do i=1,natom
	    vx(i)=vx(i)*kfG
	    vy(i)=vy(i)*kfG
	    vz(i)=vz(i)*kfG
	  enddo
	endif

	ekin=ekinT

	return
	end
*-------------------------------------------------------------------------------
