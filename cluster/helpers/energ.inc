*===============================================================================
*	ENERGY
*===============================================================================
* file: energ.inc
*-------------------------------------------------------------------------------

	integer ecltim
	double precision ekin,epot
	double precision ekinA(mxsorts)
	double precision ewall,ewallsum,fwall,fwallsum
	character*1 ethst,eclth
	common /ENERG/ ekin,epot,ekinA,
     $		       ewall,ewallsum,fwall,fwallsum,
     $		       ecltim,
     $		       ethst,eclth

*-------------------------------------------------------------------------------
* 	Common Block ENERG (energ.inc)
* ekin		: total kinetic energy at time t
* epot		: total potential energy at time t
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
