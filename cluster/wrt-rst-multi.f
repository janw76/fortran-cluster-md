*===============================================================================
*	WRITE MULTIPLE RESTART FILES
*===============================================================================
* file: wrt-rst-multi.f
*-------------------------------------------------------------------------------
* This routine writes the relevant data for a restart-run of the program:
* positions, velocities, accelerations, energies, etc.
* The rest of the data is still supplied by the input-3d file!
* There will be a different filename every writing step
*-------------------------------------------------------------------------------

	subroutine wrtRST(step,fnum)
	
	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'energ.inc'
	include 'filep.inc'

	integer step,fnum
	integer st
	integer i, count
	
	character*4 num,num0

*-------------------------------------------------------------------------------
* step		: time step number
* fnum		: file number for output file
* st		: Startposition fuer inserting
*...............................................................................
* i		: temporary variable for mxa
* count		: Counter for the Filename
*...............................................................................
* num		: insert of the filename
* num0		: searchpattern
*-------------------------------------------------------------------------------
	
	st=0

	count=step/ftrst
	num0='####'

	do i=1,77
	  if (fxrst(i:i+3).eq.num0) st=i
	enddo

	if (st.eq.0) then
	  write (*,*) '--- ERROR --- No pattern #### for'//
     $		      ' sequence number found!'
	  stop
	endif	

	write(num,'(I4.4)') count
	fxrst(st:st+3)=num

*	write(*,*) fxrst

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

	fxrst(st:st+3)=num0

	return
	end
*-------------------------------------------------------------------------------
