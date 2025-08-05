*===============================================================================
* 	THEROMSTAT BY VELOCITY SCALING
*===============================================================================
* file: scalev-m.f
*-------------------------------------------------------------------------------
* This routine scales the velocities of the particles to the correct
* kinetic energy, and thus sets the correct temperature. A better
* alternative would be using the Nose-Hoover thermostat, but this
* essentially does the same, just has a more fancy way of arriving at the
* scaling factor.
* last changed:	14-jul-2005 Jan Wedekind
* change:	corrected overcounting of atoms in lines: do j=n,n+apnum(i)-1
*-------------------------------------------------------------------------------

	subroutine thermostat(step)

	implicit none
	
	include 'const.inc'
	include 'energ.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'

	integer step
	integer i,j,n
	double precision kf,ek,ekinT

*-------------------------------------------------------------------------------
* step		: time step number
* i		: standard loop variable
* kf		: correction factor for velocities
* ek		: kinetic energy of given atsort(i)
*-------------------------------------------------------------------------------
	n=1
	ekin=0.0d0

	do i=1,atsorts
	  ekinT=1.5d0*cikb*att*apnum(i)
	  ek=0.0d0
	  do j=n,n+apnum(i)-1
	    ek=ek+0.5d0*apmas(i)*
     $	       (vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j))
	  enddo
	  ekin=ekin+ek

	  kf=dsqrt(ekinT/ek)

*	if (mod(step,20000).eq.0) then
*	  write (42,'(2X,F15.10,2X,F8.4,$)') ek,ek/(apnum(i)*1.5d0*cikb)
*	write (*,'(3(2X,F7.4))') 2d0*ekin/(3*natom*cikb),
*    $ 2d0*ekinT/(3*natom*cikb),kf
*	
*	endif
*
*	write (88,*) (ekin-ekinT)/(cikb*att)*(ekin-ekinT)/(cikb*att)

	  if (apth(i).eq.'X'.or.step.eq.-1) then
	  do j=n,n+apnum(i)-1
	      vx(j)=vx(j)*kf
	      vy(j)=vy(j)*kf
	      vz(j)=vz(j)*kf
	    enddo
	    ekin=ekin-ek+ekinT
	  endif
	  n=n+apnum(i)
	enddo

	return
	end
*-------------------------------------------------------------------------------
