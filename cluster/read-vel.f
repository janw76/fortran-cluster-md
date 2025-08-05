*===============================================================================
*	WRITE VELOCITIES
*===============================================================================
* file: read-vel.f
* date: 21-sep-2001
*-------------------------------------------------------------------------------
* This routine writes the velocities of the atoms to a file
*-------------------------------------------------------------------------------

	subroutine readVEL(step,fnum)

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

	read (fnum,'(I15)') natom
	read (fnum,'(I15)') step
	do i=1,natom
	  read (fnum,10) ats,vx(i),vy(i),vz(i)
	enddo
10	format (A2,3(3X,F15.10))
	
	return
	end
*-------------------------------------------------------------------------------
