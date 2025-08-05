*===============================================================================
*	READ RESTART-FILE
*===============================================================================
* file: read-rst.f
*-------------------------------------------------------------------------------
* This routine reads the important simulation data from the program
* restart file: positions, velocities, (and further derivatives) etc (see 
* below)
* The rest of the simulation data (like atomic parameters, temperature,
* etc.) is still given by the usual input file, which only controls wether
* a restart file is read or not.
*-------------------------------------------------------------------------------

	subroutine readRST(step,fnum)
	
	implicit none
	
	include 'const.inc'
	include 'filep.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'energ.inc'
	include 'timep.inc'
	
	integer step,fnum
	integer atoms
	character*1 in

*-------------------------------------------------------------------------------
* step		: time step number
* fnum		: file number for input
*...............................................................................
* atoms		: Number of atoms in restart-file at compile time
* in		: input variable for safety question
*-------------------------------------------------------------------------------

	open (fnum,file=fxinp,form='unformatted')
	
	read (fnum) atoms
	read (fnum) ekin,epot,ekin,epot,ekin
	read (fnum) ewall,ewallsum,fwall,fwallsum

	read (fnum) x0,y0,z0
	read (fnum) x1,y1,z1
	read (fnum) x2,y2,z2
	read (fnum) x3,y3,z3
	read (fnum) x4,y4,z4
	read (fnum) x5,y5,z5


	if (atoms.ne.mxa) then
	  write (*,*) 'This Restart-File is written by a version'
	  write (*,*) 'compiled for',step,' atoms. This version is'
	  write (*,*) 'compiled for',mxa,' atoms.'
	  write (*,'(A,$)') 'Read anyway? (yY, nN): '
	  read(*,'(A)') in
	  if (in.ne.'y'.and.in.ne.'Y') stop
	endif
	step=btim

	close (fnum)

	call neighbor(step)

	return
	end
*-------------------------------------------------------------------------------
