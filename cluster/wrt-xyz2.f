*===============================================================================
*	WRITE COORDINATE FILE
*===============================================================================
* file: wrt-xyz2.f
*-------------------------------------------------------------------------------
* This routine writes a coordinate file that is readable by visualization
* programs. The routine works for both box and periodic boundary
* conditions and switches.
* The atoms are filtered with the cluster sizes for easier visulaization
* in SciAn. Only the biggest cluster is shown!
* this routine uses the sign of the timestep as a switch for displaying
* the monomers: step<0 : don't show monomers, >0 show 'em all
*-------------------------------------------------------------------------------

	subroutine wrtXYZ(step,fnum)

	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'boxpp.inc'
	include 'cvtpp.inc'
	include 'timep.inc'

	integer step,fnum
	integer i,monflag
	double precision hx,hy,hz,tx,ty,tz
		
*-------------------------------------------------------------------------------
* step		: time step number
* fnum		: file number for output
*...............................................................................
* i		: standard loop variable
* hx,hy,hz	: temporary variables for box center coordinates
* tx,ty,tz	: temporary variables for coordinates
* monflag	: switch on/off monomers
*		  1: Monomers off, 0: Monomers on
*-------------------------------------------------------------------------------


	monflag=0
	if (step.lt.0) then
	  step=-step
	  monflag=1
	endif

	call stoddard()

	hx=boxx/2
	hy=boxy/2
	hz=boxz/2

	monflag=0
	do i=1,natom
	  if (csi(i).gt.0) monflag=i
	enddo

	if (bcor.eq.'-'.or.step.ne.btim) then
	  write (fnum,*) monflag*csi(monflag)
	  write (fnum,*) step
	else
	  write (fnum,*) monflag*csi(monflag)+8
	  write (fnum,*) step
	  write (fnum,20) 'H ',-hx,-hy,-hz
	  write (fnum,20) 'H ',hx,-hy,-hz
	  write (fnum,20) 'H ',-hx,hy,-hz
	  write (fnum,20) 'H ',hx,hy,-hz
	  write (fnum,20) 'H ',-hx,-hy,hz
	  write (fnum,20) 'H ',hx,-hy,hz
	  write (fnum,20) 'H ',-hx,hy,hz
	  write (fnum,20) 'H ',hx,hy,hz
	endif

	do i=1,natom
	  if (cs(cl(i)).eq.monflag) then
	    tx=x(i)
	    ty=y(i)
	    tz=z(i)
	    if (bperx.eq.'X') tx=dmod(dmod(tx,boxx)+boxx,boxx)
	    if (bpery.eq.'X') ty=dmod(dmod(ty,boxy)+boxy,boxy)
	    if (bperz.eq.'X') tz=dmod(dmod(tz,boxz)+boxz,boxz)
	    write (fnum,20) ats,tx-hx,ty-hy,tz-hz
	  endif
	enddo

20	format (A2,3(2X,F10.4))
	return
	end
*-------------------------------------------------------------------------------
