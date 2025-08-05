*===============================================================================
* 	3rd ORDER GEAR-INTEGRATOR
*===============================================================================
* file: int-gear3.f
*-------------------------------------------------------------------------------
* 3th order Gear algorithm. The factors Fx2 are calculated for minimum
* loss of precision and maximum speed, with R2 (=xERR) calculated without
* dt^2/2. See J.M. Haile, 'Molecular Dynamics Simulation', John Wiley &Sons,
* 1992, p. 161/162
*-------------------------------------------------------------------------------

	subroutine integrator(step)
	
	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'energ.inc'

	integer i,step
	double precision dt1,dt2,dt3
	double precision f02,f12,f22,f32
	double precision axp(mxa),ayp(mxa),azp(mxa)
	double precision xerr,yerr,zerr

*-------------------------------------------------------------------------------
* step		: time step
*...............................................................................
* i		: standard loop variable
* dt1		: dt^1
* dt2		: dt^2
* dt3		: dt^3
* f02		: Factor for Gear-Algorithm
* f12		: Factor for Gear-Algorithm
* f22		: Factor for Gear-Algorithm
* f32		: Factor for Gear-Algorithm
* axp()		: predicted acceleration in X-direction
* ayp()		: predicted acceleration in Y-direction
* azp()		: predicted acceleration in Z-direction
* xerr		: Error of predicted position in X-direction
* yerr		: Error of predicted position in Y-direction
* zerr		: Error of predicted position in Z-direction
*-------------------------------------------------------------------------------

	dt1=dt
	dt2=dt*dt/2.0d0
	dt3=dt*dt*dt/6.0d0

	do i=1,natom
	  x0(i)=x0(i)+x1(i)*dt1+x2(i)*dt2+x3(i)*dt3
	  y0(i)=y0(i)+y1(i)*dt1+y2(i)*dt2+y3(i)*dt3
	  z0(i)=z0(i)+z1(i)*dt1+z2(i)*dt2+z3(i)*dt3
	  x1(i)=x1(i)+x2(i)*dt1+x3(i)*dt2
	  y1(i)=y1(i)+y2(i)*dt1+y3(i)*dt2
	  z1(i)=z1(i)+z2(i)*dt1+z3(i)*dt2
	  x2(i)=x2(i)+x3(i)*dt1
	  y2(i)=y2(i)+y3(i)*dt1
	  z2(i)=z2(i)+z3(i)*dt1
	  axp(i)=x2(i)
	  ayp(i)=y2(i)
	  azp(i)=z2(i)
	enddo

	epot=0.0d0
	do i=1,natom
	  x2(i)=0.0d0
	  y2(i)=0.0d0
	  z2(i)=0.0d0
	enddo

	call AccAtom()
	call AccWall()

	f02=dt*dt*1.0d0/12.0d0
	f12=dt*5.0d0/12.0d0
	f32=1.0d0/dt

	ekin=0.0d0
	do i=1,natom
	  xerr=x2(i)-axp(i)
	  yerr=y2(i)-ayp(i)
	  zerr=z2(i)-azp(i)
	  x0(i)=x0(i)+xerr*f02
          x1(i)=x1(i)+xerr*f12
          x3(i)=x3(i)+xerr*f32
          y0(i)=y0(i)+yerr*f02
          y1(i)=y1(i)+yerr*f12
          y3(i)=y3(i)+yerr*f32
          z0(i)=z0(i)+zerr*f02
          z1(i)=z1(i)+zerr*f12
          z3(i)=z3(i)+zerr*f32
	  ekin=ekin+0.5d0*atm*(x1(i)*x1(i)+y1(i)*y1(i)+z1(i)*z1(i))
	enddo

	return
	end
*-------------------------------------------------------------------------------
