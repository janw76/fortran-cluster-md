*===============================================================================
*	OPTIMIZATION STAISTICS
*===============================================================================
* file: dmmpp.inc
*-------------------------------------------------------------------------------

	double precision dmmx,dmmy,dmmz,dmmav
	double precision dmmax,dmms,dmmsu
	double precision dmmh(0:255)
	integer dmmst,dmmi

	common /DMMPP/ dmmh,dmmx,dmmy,dmmz,dmms,dmmsu,
     $		       dmmax,dmmav,dmmst,dmmi

*-------------------------------------------------------------------------------
*	Common Block DMMPP (dmmpp.inc)
* dmmx,dmmy,dmmz: used for calculation of the maximum displacement 
* dmms		: sum of maximum displacement over nbtim steps
* dmmax		: Maximum of dmms
* dmmsu		: total sum of dmm for average/moving average
* dmmh()	: history field for moving average of dmms
* dmmi		: Index for dmmh()
* dmmst		: Number of steps in dmms-run
* dmmav		: (moving) average for dmm
*-------------------------------------------------------------------------------

