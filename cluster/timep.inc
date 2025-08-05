*===============================================================================
*	TIME PARAMETERS
*===============================================================================
* file: timep.inc
*-------------------------------------------------------------------------------

	double precision dt
	integer ntim,btim,sctim,nbtim,chgNtim
	common /TIMEP/ dt,ntim,btim,sctim,nbtim,chgNtim

*-------------------------------------------------------------------------------
*	Common Block TIMEP (timep.inc)
* dt		: simulation time step (fs, 10^-15 s)
* ntim		: number of time steps
* btim		: starting value of time step number
* sctim		: number of steps between screen output
* nbtim		: number of steps between neighbor-list updates
* chgNtim	: number of steps between changes of natom
*-------------------------------------------------------------------------------
