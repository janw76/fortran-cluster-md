*===============================================================================
*	ENERGY
*===============================================================================
* file: energ.inc
*-------------------------------------------------------------------------------

	integer ecltim
	double precision ekin,epot,H
	double precision ekinA(mxsorts)
	double precision ewall,ewallsum,fwall,fwallsum,sumv2
	character*1 ethst,eclth
	common /ENERG/ ekin,sumv2,epot,ekinA,H,
     $		       ewall,ewallsum,fwall,fwallsum,
     $		       ecltim,
     $		       ethst,eclth

*-------------------------------------------------------------------------------
* 	Common Block ENERG (energ.inc)
* ekin		: total kinetic energy at time t
* sumv2		: Sum of all v(i)^2
* epot		: total potential energy at time t
* H		: Hamiltonian of the system at time t
* ekinA()	: kinetic energy per atom sort
* ewall		: energy of interaction with the wall
* ewallsum	: summed wall-energy for pressure calculation
* fwall		: force with the wall for pressure calculation
* ethst		: flag for Thermostat 
*		  X: switched on
*		  -: switched off
* eclth		: flag for Cluster Thermostat 
*		  X: switched on
*		  -: switched off
* ecltim	: Timer for cluster thermostat
*-------------------------------------------------------------------------------
