*===============================================================================
*	WRITE COORDINATE FILE
*===============================================================================
* file: wrt-xyz3.f
*-------------------------------------------------------------------------------
* This routine writes a coordinate file that is readable by visualization
* programs. The routine works for both box and periodic boundary
* conditions and switches.
* The atoms are filtered with the cluster sizes for easier visulaization
* in SciAn. Only clusters beyond a certain size are shown!
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
	integer i,gas,minsize
	double precision hx,hy,hz,tx,ty,tz

*-------------------------------------------------------------------------------
* step		: time step number
* fnum		: file number for output
*...............................................................................
* i		: standard loop variable
* hx,hy,hz	: temporary variables for box center coordinates
* tx,ty,tz	: temporary variables for coordinates
* gas		: number of atoms not displayd
* minsize	: minimum size for display
*-------------------------------------------------------------------------------


	minsize=16
	step=abs(step)

	call stoddard()

	hx=boxx/2
	hy=boxy/2
	hz=boxz/2

	gas=0
	do i=1,minsize-1
	  gas=gas+i*csi(i)
	enddo

	if (bcor.eq.'-'.or.step.ne.btim) then
	  write (fnum,*) natom-gas
	  write (fnum,*) step
	else
	  write (fnum,*) natom-gas+8
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
	  if (cs(cl(i)).ge.minsize) then
	    tx=x(i)
	    ty=y(i)
	    tz=z(i)
	    if (bperx.eq.'X') tx=dmod(dmod(tx,boxx)+boxx,boxx)
	    if (bpery.eq.'X') ty=dmod(dmod(ty,boxy)+boxy,boxy)
	    if (bperz.eq.'X') tz=dmod(dmod(tz,boxz)+boxz,boxz)
	    write (fnum,20) ats,tx-hx,ty-hy,tz-hz
	  endif
	enddo

20	format (A2,3(2X,F10.4),2(2X,I6))
	return
	end
*-------------------------------------------------------------------------------
