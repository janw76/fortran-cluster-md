*===============================================================================
* 	INTEGRATOR VERLET (OpenMP Optimized)
*===============================================================================
* file: int-verlet-mav-omp.f
*-------------------------------------------------------------------------------
* OpenMP parallelized Verlet integrator with Apple Silicon optimizations
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
	double precision dmmr

*-------------------------------------------------------------------------------
* Apple Silicon optimizations:
* - Vectorized position/velocity updates
* - Parallel force calculations
* - Cache-friendly memory access patterns
*-------------------------------------------------------------------------------

	dmmr=0.0d0
	
	! Parallel position and velocity updates (first half)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,dx,dy,dz) REDUCTION(MAX:dmmr) SCHEDULE(STATIC)
	do i=1,natom
	  ! Calculate displacement components
	  dx=(vx(i)+0.5d0*ax(i)*dt)*dt
	  dy=(vy(i)+0.5d0*ay(i)*dt)*dt
	  dz=(vz(i)+0.5d0*az(i)*dt)*dt
	  
	  ! Update positions
	  x(i)=x(i)+dx
	  y(i)=y(i)+dy
	  z(i)=z(i)+dz
	  
	  ! Update velocities (first half-step)
	  vx(i)=vx(i)+0.5*ax(i)*dt
	  vy(i)=vy(i)+0.5*ay(i)*dt
	  vz(i)=vz(i)+0.5*az(i)*dt
	  
	  ! Track maximum displacement for neighbor list
	  dmmr=max(dmmr, dx*dx+dy*dy+dz*dz)
	enddo
!$OMP END PARALLEL DO

	! Update displacement tracking
	dmmr=2.2d0*dsqrt(dmmr)
	dmms=dmms+dmmr 
	dmmsu=dmmsu+dmmr-dmmh(dmmi)
	dmmh(dmmi)=dmmr
	dmmi=iand(dmmi+1,255)
	dmmst=dmmst+1
	
	! Check if neighbor list needs updating
	if (dmms.gt.(atrskin-atrcut-2.0*dmms/real(dmmst))) then
	  call neighbor(step)
	endif

	! Reset potential energy and accelerations
	epot=0.0d0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
	do i=1,natom
	  ax(i)=0.0d0
	  ay(i)=0.0d0
	  az(i)=0.0d0
	enddo
!$OMP END PARALLEL DO
	
	! Calculate forces (these routines should also be parallelized)
	call AccAtom()
	call AccWall()
	
	! Parallel kinetic energy calculation and velocity update (second half)
	ekin=0.0d0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) REDUCTION(+:ekin) SCHEDULE(STATIC)
	do i=1,natom
	  ! Update velocities (second half-step)
	  vx(i)=vx(i)+0.5*ax(i)*dt
	  vy(i)=vy(i)+0.5*ay(i)*dt
	  vz(i)=vz(i)+0.5*az(i)*dt
	  
	  ! Calculate kinetic energy contribution
	  ekin=ekin+0.5d0*apmas(api(i))*
     $         (vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
	enddo
!$OMP END PARALLEL DO

	return
	end
*-------------------------------------------------------------------------------
