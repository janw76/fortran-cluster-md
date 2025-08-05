*===============================================================================
*	FILE PARAMETERS
*===============================================================================
* file: filep.inc
*-------------------------------------------------------------------------------

	character*40 fxyz,fvel,facc,frst,fcsi,fave,finp,frat
        character*80 fxxyz,fxvel,fxacc,fxrst,fxcsi,fxave,fxinp,fxrat
	character*40 fpath
	character*1  frctrl
	integer ftxyz,ftvel,ftacc,ftrst,ftcsi,ftave,ftrat
	common /FILEP/ fpath,fxyz,fvel,facc,frst,fcsi,fave,finp,frat,
     $		       ftxyz,ftvel,ftacc,ftrst,ftcsi,ftave,fxinp,fxrat,
     $		       fxxyz,fxvel,fxacc,fxrst,fxcsi,fxave,ftrat,
     $		       frctrl

*-------------------------------------------------------------------------------
* 	Common Block FILEP (filep.inc)
* fpath		: path for all In/Outputfiles, relative to the binary
* fxyz		: xyz-File for SciAn or other visualisation programs
* fvel		: file with velocities of particles (corresponding to XYZ)
* facc		: file with acceleratons of particles (corresponding to XYZ)
* frst		: file with all values for a restart of the program
* fcsi		: file with cluster distribution vs. time
* fave		: file with temperature, pressure, Ekin, Epot... vs. time
* frat		: file with evaporation/condensation rates of a cluster
* finp		: filename for input-restart-file
* ftxyz		: Time Interval vor output of xyz-file
* ftvel		: Time Interval vor output of vel-file
* ftacc		: Time Interval vor output of acc-file
* ftrst		: Time Interval vor output of rst-file
* ftcsi		: Time Interval vor output of csi-file
* ftave		: Time Interval vor output of ave-file
* ftrat		: Time Interval vor output of rat-file
* fxxyz		: xyz-filename expanded with pathname
* fxvel 	: vel-filename expanded with pathname
* fxacc 	: acc-filename expanded with pathname
* fxrst 	: rst-filename expanded with pathname
* fxcsi 	: csi-filename expanded with pathname
* fxave 	: ave-filename expanded with pathname
* fxinp		: inp-filename expanded with pathname
* fxrat		: rat-filename expanded with pathname
* frctrl	: use a restart-file as starting point
*		  X: use restart-file
*		  -: don't use restart-file
*-------------------------------------------------------------------------------
