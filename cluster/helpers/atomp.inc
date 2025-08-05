*===============================================================================
*	ATOMIC PARAMETERS
*===============================================================================
* file: atomp.inc
*-------------------------------------------------------------------------------

	integer natom,nnatom,atsorts
	integer atnlist(mxnlist),atnidx(mxa)
	double precision atm,atr,ate,att
	double precision atpr,atpa,atskin
	double precision atrcut,atrskin
	character*2 ats
	character*1 rndctrl,chgNctl

	integer api(mxa),apnum(mxsorts),apnnum(mxsorts)
	double precision apsig(mxsorts),apeps(mxsorts)
	double precision apmas(mxsorts)
	double precision appr(mxsorts),appa(mxsorts)
	character*2 apsym(mxsorts)
	character*1 apth(mxsorts),apNch(mxsorts)

	equivalence (ats,apsym(1)),(atm,apmas(1))
	equivalence (atr,apsig(1)),(ate,apeps(1))
	equivalence (atpr,appr(1)),(atpa,appa(1))
	equivalence (nnatom,apnnum(1)),(chgNctl,apNch(1))

	common /ATOMP/ natom,atsorts,apnum,apmas,apsig,apeps,
     $		       appr,appa,atskin,atrskin,atrcut,att,
     $		       atnlist,atnidx,api,apnnum,
     $		       apsym,apth,apNch,rndctrl

*-------------------------------------------------------------------------------
*	Common Block ATOMP (atomp.inc)
* natom		: total number of particles in the system
* atsorts	: number of particle types 
* api()		: sort of particle i
* apnum()	: number of particles of sort i
* apnnum()	: new number of particles of sort i
* apmas()	: mass of particle of sort i
* apsig()	: diameter of particle of sort i
* apeps()	: potential depth for particle of sort i
* apsym()	: atomic symbol of particle of sort i
* apth()	: flag for thermostat connection of sort i
* apNch()	: flag for atom number change of sort i
* appr()	: exponent for repulsive part of LJ-potential
* appa()	: exponent for attractive part of LJ-potential
* atm		: synonym for mass of sort 1
* atr		: synonym for diameter of sort 1
* ate		: synonym for potential depth of sort 1
* ats		: synonym for symbol of sort 1
* atpr		: synonym for repulsive exponent for sort 1
* atpa		: synonym for attractive exponent for sort 1
* att		: system temperature
* atskin	: radius for determination of the cluster size
* atrcut	: radius for cutoff of the potential 
* atrskin	: radius for neighborlist
* atnlist()	: neighborlist
* atnidx()	: atom-index for neighborlist (for stoddard-n)
* rndctlr	: Control for initial pure random velocities
*-------------------------------------------------------------------------------
