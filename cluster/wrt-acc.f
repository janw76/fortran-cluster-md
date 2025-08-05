*===============================================================================
*	WRITE ACCELERATIONS
*===============================================================================
* file: wrt-acc.f
*-------------------------------------------------------------------------------
* This routine writes the accelerations of the particles.
*-------------------------------------------------------------------------------

	subroutine wrtACC(step,fnum)

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'

	integer step,fnum
	integer i
		
*-------------------------------------------------------------------------------
* step		: time step number
* fnum		: file number for output
*...............................................................................
* i		: standard loop variable
* dt2		: squar of time step
*-------------------------------------------------------------------------------

	write (fnum,'(I15)') natom
	write (fnum,'(I15)') step
	do i=1,natom
	  write (fnum,10) ats,ax(i),ay(i),az(i)
	enddo
10	format (A2,3(3X,F15.10))
	
	return
	end
*-------------------------------------------------------------------------------
