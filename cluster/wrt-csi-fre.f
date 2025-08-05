*===============================================================================
*	WRITE CLUSTER SIZES FOR FRENKEL
*===============================================================================
* file: wrt-csi-fre.f
*-------------------------------------------------------------------------------
* This routine writes the Stillinger and the tenWolde/Frenkel cluster 
* distribution in the system into separate files
*
* 05/10/07 (JW):
*   - Changed the first-line output for tW/F to give zeroes if no monomers
*-------------------------------------------------------------------------------

	subroutine wrtCSI(step,fnum)

	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'cvtpp.inc'

	integer step,fnum,i,fnum2
		
*-------------------------------------------------------------------------------
* step		: time step number
* fnum		: file number for output
* i		: standard loop variable
*-------------------------------------------------------------------------------
	fnum2 = fnum + 3
	if(cstep.ne.step) call stoddard(step)

	write (fnum,'(I15)') step
	write (fnum2,'(I15)') step

*** First line output for tW/F
	if (csif(0).eq.0) then
	  write (fnum2,10) 0,csif(0),0.0d0,0.0d0,0.0d0
	else
	  write (fnum2,10) 0,csif(0),0.0d0,
     $	DBLE(ceksf(0))/DBLE(csif(0)),
     $  DBLE(ceksf(0))/DBLE(csif(0)*1.5d0*cikb)
	endif

** Here's the rest and Stillinger	
	do i=1,natom
	  if (csi(i).ne.0) then
	    write (fnum,10) i,csi(i),REAL(cgsa(i))/REAL(csi(i)),
     $	ceks(i)/(i*csi(i)),
     $  ceks(i)/(i*csi(i)*1.5d0*cikb)
	  endif
	  
	  if (csif(i).ne.0) then
	    write (fnum2,10) i,csif(i),REAL(cgsfa(i))/REAL(csif(i)),
     $	REAL(ceksf(i))/REAL(i*csif(i)),
     $  REAL(ceksf(i))/REAL(i*csif(i)*1.5d0*cikb)
	  endif
	enddo
10	format (2(3X,I7),3X,F7.3,3X,F17.15,3X,F10.5)

	return
	end
*-------------------------------------------------------------------------------
