*===============================================================================
*	ATOMIC COORDINATES AND DERIVATIVES
*===============================================================================
* file: atomc.inc
*-------------------------------------------------------------------------------

	double precision  x(mxa), y(mxa), z(mxa)
	double precision vx(mxa),vy(mxa),vz(mxa)
	double precision ax(mxa),ay(mxa),az(mxa)
	double precision x0(mxa),y0(mxa),z0(mxa)
	double precision x1(mxa),y1(mxa),z1(mxa)
	double precision x2(mxa),y2(mxa),z2(mxa)
	double precision x3(mxa),y3(mxa),z3(mxa)
	double precision x4(mxa),y4(mxa),z4(mxa)
	double precision x5(mxa),y5(mxa),z5(mxa)
	equivalence (x ,x0),(y ,y0),(z ,z0)
	equivalence (vx,x1),(vy,y1),(vz,z1)
	equivalence (ax,x2),(ay,y2),(az,z2)
	common /ATOMC/ x0,y0,z0,x1,y1,z1,x2,y2,z2,
     $		       x3,y3,z3,x4,y4,z4,x5,y5,z5

*-------------------------------------------------------------------------------
* 	Common Block ATOMC (atomc.inc)
* x0(),y0(),z0(): Positions of the particles
* x1(),y1(),z1(): 1st derivative in time of positions ('velocity')
* x2(),y2(),z2(): 2nd derivative in time of positions ('acceleration')
* x3(),y3(),z3(): 3rd derivative in time of positions ('jerk')
* x4(),y4(),z4(): 4th derivative in time of positions
* x5(),y5(),z5(): 5th derivative in time of positions
*  x(), y(), z(): alternate names for x0(),y0(),z0()
* vx(),vy(),vz(): alternate names for x1(),y1(),z1()
* ax(),ay(),az(): alternate names for x2(),y2(),z2()
*-------------------------------------------------------------------------------
