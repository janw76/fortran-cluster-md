*===============================================================================
*	WRITE FIRST PASSAGE TIMES
*===============================================================================
* file: wrt-fpt.f
*-------------------------------------------------------------------------------
* 05/10/07 (JW)
* This routine writes the first passage times of the Stillinger and the 
* tenWolde/Frenkel cluster 
* writes out in nanoseconds!
*-------------------------------------------------------------------------------

	subroutine wrtFPT()

	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'cvtpp.inc'	
	include 'timep.inc'
	
	integer i
		
*-------------------------------------------------------------------------------
* step		: time step number
* fnum		: file number for output
* i		: standard loop variable
*-------------------------------------------------------------------------------
	
	do i = 1,natom
	   if (cfpt(i).ne.-1) then
	      if (cfptf(i).eq.-1) cfptf(i) = REAL(cfptf(i))*REAL(dt*1.0d6)
	write (69,10) i,REAL(cfpt(i))/(dt*1.0d6),REAL(cfptf(i))/REAL(dt*1.0d6)
	   endif	   
	enddo

10	format (I7,3X,F12.4,3X,F12.4)
	
	return
	end
*-------------------------------------------------------------------------------
