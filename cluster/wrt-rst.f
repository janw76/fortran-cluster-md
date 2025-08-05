*===============================================================================
*	WRITE THE RESTART FILE
*===============================================================================
* file: wrt-rst.f
*-------------------------------------------------------------------------------
* This routine writes the relevant data for a restart-run of the program:
* positions, velocities, accelerations, energies, etc.
* The rest of the data is still supplied by the input-3d file!
*-------------------------------------------------------------------------------

	subroutine wrtRST(step,fnum)
	
	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'energ.inc'
	include 'filep.inc'

	integer step,fnum
	integer i

*-------------------------------------------------------------------------------
* step		: time step number
* fnum		: file number for output file
*...............................................................................
* i		: temporary variable for mxa
*-------------------------------------------------------------------------------

	open (fnum,file=fxrst,form='unformatted')
	
	i=step
	i=mxa
	write (fnum) i
	write (fnum) ekin,epot,ekin,epot,ekin
	write (fnum) ewall,ewallsum,fwall,fwallsum

	write (fnum) x0,y0,z0
	write (fnum) x1,y1,z1
	write (fnum) x2,y2,z2
	write (fnum) x3,y3,z3
	write (fnum) x4,y4,z4
	write (fnum) x5,y5,z5

	close (fnum)

	return
	end
*-------------------------------------------------------------------------------
