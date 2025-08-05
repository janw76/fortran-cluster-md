*===============================================================================
* 	INTEGRATOR VERLET Nose-Hoover
*===============================================================================
* file: int-verlet-nosehoover-mav.f
*-------------------------------------------------------------------------------
* This routine integrates the equations of motion with the Verlet-method
* 26-apr-2006: added Nose-Hoover
* if Nose-Hoover is specified as Thermostat with the 'N' flag in the .3d file
* Nose-Hoover factor is used. Otherwise S=PSI=0
* NOSE HOOVER ONLY WORKS WITH ONE SORT OF ATOMS!!!
*-------------------------------------------------------------------------------

	subroutine integrator(step)
	
	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'energ.inc'
	include 'dmmpp.inc'

	integer i,step,n,iter
	double precision dx,dy,dz,ek,dt2,dth,err
	double precision vxs,vys,vzs,mas,imas
	double precision dmmr, pso, psn,delps,di
	double precision vxn(mxa),vyn(mxa),vzn(mxa)
	double precision vxo(mxa),vyo(mxa),vzo(mxa)
	double precision bx(mxa),by(mxa),bz(mxa)
	logical ready

*-------------------------------------------------------------------------------
* step		: time step
* i,n,iter	: standard loop variable
* dx,dy,dz	: components of the displacement
* dmmr 		: Maximum Displacement Vector
* vxn(i),vxo(i)	: temporary vel.
* bx(i),delps,di: temporary variables for iteration
* psn, pso	: temporary PSI values
*-------------------------------------------------------------------------------

	dt2 = dt*dt/2
	dth = dt/2
	sumv2 = 0
	dmmr=0.0d0

*** calculate new pos and intermediate vel
	do i=1,natom
	   dx = dt*VX(i) + dt2 * ( ax(i) - PSI*VX(i) )
	   dy = dt*VY(i) + dt2 * ( ay(i) - PSI*VY(i) )
	   dz = dt*VZ(i) + dt2 * ( az(i) - PSI*VZ(i) )
	   x(i) = x(i) + dx
	   y(i) = y(i) + dy
	   z(i) = z(i) + dz
	   sumv2 = sumv2 + VX(i)*VX(i)+VY(i)*VY(i)+VZ(i)*VZ(i)
	   VX(i) = VX(i) + dth * ( ax(i) - PSI*VX(i) )
	   VY(i) = VY(i) + dth * ( ay(i) - PSI*VY(i) )
	   VZ(i) = VZ(i) + dth * ( az(i) - PSI*VZ(i) )
	   dmmr = max(dmmr,dx*dx+dy*dy+dz*dz)
	enddo
	
*** calculating Nose-Hoover coeff for t+dt
	if (ethst.eq.'N') then
	   s = s + PSI*dt + (apmas(1)*sumv2 - 3.0d0*natom*cikb*att) * dt2/Q
	   PSI = PSI + (apmas(1)*sumv2 - 3.0d0*natom*cikb*att) * dth/Q
	endif

*** maximum displacement vector from all single atoms, safety factor of 10%
	dmmr=2.2d0*dsqrt(dmmr)
	dmms=dmms+dmmr 
	dmmsu=dmmsu+dmmr-dmmh(dmmi)
	dmmh(dmmi)=dmmr
	dmmi=iand(dmmi+1,255)
	dmmst=dmmst+1
	if (dmms.gt.(atrskin-atrcut-2.0*dmms/real(dmmst))) then
	  call neighbor(step)
	endif

*** reset epot and acc to zero, call forces
	epot=0.0d0
	do i=1,natom
	   ax(i)=0.0d0
	   ay(i)=0.0d0
	   az(i)=0.0d0
	enddo
	call AccAtom()
	call AccWall()
	
*** final update of velocities (iterative)
	ekin=0.0d0
	n=1
	err=1.0D-10	
	sumv2 = 0
	do i=1,natom
	   vxn(i) = VX(i)
	   vyn(i) = VY(i)
	   vzn(i) = VZ(i)
	   sumv2 = sumv2 + vxn(i)*vxn(i) + vyn(i)*vyn(i) + vzn(i)*vzn(i)
	enddo

	psn = PSI
	ready = .FALSE.

*** begin iteration ***
	do iter=1,1000
	   pso = psn
	   delps = 0
	   
	   do i=1,natom
	      vxo(i) = vxn(i)
	      vyo(i) = vyn(i)
	      vzo(i) = vzn(i)
	      bx(i) = -dth*(ax(i)-pso*vxo(i)) - (VX(i)-vxo(i))
	      delps = delps + bx(i)*apmas(1)*vxo(i)*dt/Q
	      by(i) = -dth*(ay(i)-pso*vyo(i)) - (VY(i)-vyo(i))
	      delps = delps + by(i)*apmas(1)*vyo(i)*dt/Q
	      bz(i) = -dth*(az(i)-pso*vzo(i)) - (VZ(i)-vzo(i))
	      delps = delps + bz(i)*apmas(1)*vzo(i)*dt/Q
	   enddo

	   di = -(pso*dth + 1)

	   delps =delps-di*((-apmas(1)*sumv2+3*natom*cikb*att)*dth/Q-(PSI-pso))
	   delps = delps/(-dt2*apmas(1)*sumv2/Q + di)
	 
	   sumv2 = 0
	   do i=1,natom
	      vxn(i) = vxn(i) + (bx(i) + dth*vxo(i)*delps)/di
	      vyn(i) = vyn(i) + (by(i) + dth*vyo(i)*delps)/di
	      vzn(i) = vzn(i) + (bz(i) + dth*vzo(i)*delps)/di
	      sumv2 = sumv2 + vxn(i)*vxn(i) + vyn(i)*vyn(i) + vzn(i)*vzn(i)
	   enddo

	   psn = pso + delps

*** Test for convergence
	   ready = .true.
	   do i=1,natom
	      if (ABS((vxn(i)-vxo(i))/vxn(i)).GT.err) ready = .false.
	      if (ABS((vyn(i)-vyo(i))/vyn(i)).GT.err) ready = .false.
	      if (ABS((vzn(i)-vzo(i))/vzn(i)).GT.err) ready = .false.
	   enddo
*	   if (ABS((psn-pso)/psn).GT.err*100000) ready = .false.
	   if (ready) goto 10
	
	enddo
*** end iteration ***
	
	if (.not.ready) then
	   write (*,*) 'Numerical integration did not converge!'
	   write (*,*) '---> stopping simulation!'
	   write (*,*) ABS((psn-pso)/psn)
	   do i=1,natom
	      if (ABS((vxn(i)-vxo(i))/vxn(i)).GT.err) then
	       write (*,*) (ABS((vxn(i)-vxo(i))/vxn(i)))
	      endif
	      if (ABS((vyn(i)-vyo(i))/vyn(i)).GT.err) then
	       write (*,*) (ABS((vyn(i)-vyo(i))/vyn(i)))
	      endif
	      if (ABS((vzn(i)-vzo(i))/vzn(i)).GT.err) then
	       write (*,*) (ABS((vzn(i)-vzo(i))/vzn(i)))
	      endif
	   enddo
	   stop
	endif

10	continue

*** writing back final velocities and PSI
	sumv2 = 0
	do i=1,natom
	   VX(i) = vxn(i)
	   VY(i) = vyn(i)
	   VZ(i) = vzn(i)
	   sumv2 = sumv2 + VX(i)*VX(i)+VY(i)*VY(i)+VZ(i)*VZ(i)
	enddo
	
	PSI = psn
	ekin = 0.5d0*apmas(1)*sumv2

*** stopdead	
c	vxs=0.0d0
c	vys=0.0d0
c	vzs=0.0d0
c	do i=1,natom
c	  vxs=vxs+vx(i)
c	  vys=vys+vy(i)
c	  vzs=vzs+vz(i)
c	enddo
c	vxs=vxs/dble(natom)
c	vys=vys/dble(natom)
c	vzs=vzs/dble(natom)
c	
c	do i=1,natom
c	  vx(i)=vx(i)-vxs
c	  vy(i)=vy(i)-vys
c	  vz(i)=vz(i)-vzs
c	enddo
*** end stopdead

	return

	end
*-------------------------------------------------------------------------------
