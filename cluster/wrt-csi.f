*===============================================================================
*	WRITE CLUSTER SIZES
*===============================================================================
* file: wrt-csi.f
*-------------------------------------------------------------------------------
* This routine writes the cluster distribution in the system
*-------------------------------------------------------------------------------

	subroutine wrtCSI(step,fnum)

	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'cvtpp.inc'

	integer step,fnum,i
		
*-------------------------------------------------------------------------------
* step		: time step number
* fnum		: file number for output
*...............................................................................
* i		: standard loop variable
*-------------------------------------------------------------------------------

	if(cstep.ne.step) call stoddard(step)
	write (fnum,'(I15)') step
	do i=1,natom
	  if (csi(i).ne.0) then
	    write (fnum,10) i,csi(i),
     $	ceks(i)/(i*csi(i)),
     $  ceks(i)/(i*csi(i)*1.5d0*cikb)
	  endif
	enddo
10	format (2(3X,I7),3X,F17.15,3X,F10.5)

	return
	end
*-------------------------------------------------------------------------------
