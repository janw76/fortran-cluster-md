*===============================================================================
*	WRITE VELOCITIES
*===============================================================================
* file: wrt-vel.f
*-------------------------------------------------------------------------------
* This routine writes the velocities of the atos to a file
*-------------------------------------------------------------------------------

	subroutine wrtVEL(step,fnum)

	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'

	integer step,fnum
	integer i
		
*-------------------------------------------------------------------------------
* step		: time step number
* fnum		: file number for output file
*...............................................................................
* i		: standard loop variable
*-------------------------------------------------------------------------------

	write (fnum,'(I15)') natom
	write (fnum,'(I15)') step
	do i=1,natom
	  write (fnum,10) ats,vx(i),vy(i),vz(i)
	enddo
10	format (A2,3(3X,F15.10))
	
	return
	end
*-------------------------------------------------------------------------------
