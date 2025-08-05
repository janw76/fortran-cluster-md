*===============================================================================
*	CONSTANTS
*===============================================================================
* file: const.inc
*-------------------------------------------------------------------------------

	integer mxa,mxnlist,mxsorts
	double precision cokb,cona,comu,copl,cimu,cikb
	double precision pi
	parameter (mxa=4096)
	parameter (mxnlist=mxa*mxa/4)
	parameter (mxsorts=10)
	parameter (cokb=1.38066D-23)
	parameter (cona=6.022045D23)
	parameter (comu=1.66056D-27)
	parameter (copl=6.62618D-34)
	parameter (cikb=1.38066D-6)
	parameter (cimu=1.66056D00)
*	parameter (pi=4.0d0*datan(1.0d0))
	parameter (pi=3.14159265)

*-------------------------------------------------------------------------------
* mxa		: maximum number of atoms handled by the code
* mxnlist	: maximum number of entries in the neighbor list
* mxsorts	: maximum number of different species
*...............................................................................
* Various physical constants (Source: P.W. Atkins, Physikalische Chemie, 
* VCH Weinheim, 2. Auflage 1990
*...............................................................................
* cobk		: Boltzmann constant [J/K]
* cona		: Avogadro's number  [1/mol]
* comu		: Atomic mass units  [kg]
* copl		: Planck's constant  [Js]
* cikb		: Boltzmann constant, internal [10^-17 J/K]
* cimu		: Atomic mass units, internal  [10^-27 kg]
* pi		: constant pi
*-------------------------------------------------------------------------------
