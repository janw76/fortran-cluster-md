*===============================================================================
*	WRITE COORDINATE FILE
*===============================================================================
* file	  : wrt-xyz-m.f
* changed : 20-aug-2002
*-------------------------------------------------------------------------------
* This routine writes a coordinate file that is readable by visualization
* programs. The routine works for both box and periodic boundary
* conditions and switches.
* This routine can use multiple atomic flavours.
*-------------------------------------------------------------------------------

	subroutine wrtXYZ(step,fnum)

	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'boxpp.inc'
	include 'timep.inc'

	integer step,fnum,awant
	integer i
	double precision hx,hy,hz,tx,ty,tz
	character*2 aslst(mxa),asnew

	common /SPLIT/ aslst,awant,asnew
	
*-------------------------------------------------------------------------------
* step		: time step number
* fnum		: file number for output
*...............................................................................
* i		: standard loop variable
* hx,hy,hz	: temporary variables for box center coordinates
* tx,ty,tz	: temporary variables coordinates
*-------------------------------------------------------------------------------

	hx=boxx/2
	hy=boxy/2
	hz=boxz/2

	step=abs(step)
	
	if (bcor.ne.'A'.and.step.ne.btim) then
	  write (fnum,*) apnum(awant)
	  write (fnum,*) step
	elseif (bcor.eq.'A'.or.(bcor.eq.'X'.and.step.eq.btim)) then
	  write (fnum,*) apnum(awant)+8
	  write (fnum,*) step
	  write (fnum,10) 'H ',-hx,-hy,-hz
	  write (fnum,10) 'H ',hx,-hy,-hz
	  write (fnum,10) 'H ',-hx,hy,-hz
	  write (fnum,10) 'H ',hx,hy,-hz
	  write (fnum,10) 'H ',-hx,-hy,hz
	  write (fnum,10) 'H ',hx,-hy,hz
	  write (fnum,10) 'H ',-hx,hy,hz
	  write (fnum,10) 'H ',hx,hy,hz
	endif

	do i=1,natom
	
*	write(*,*) aslst(i)
*	write(*,*) apsym(awant)
	

	  if (apsym(awant).eq.aslst(i)) then
	     tx=x(i)
	     ty=y(i)
	     tz=z(i)
	     if (bperx.eq.'X') tx=dmod(dmod(tx,boxx)+boxx,boxx)
	     if (bpery.eq.'X') ty=dmod(dmod(ty,boxy)+boxy,boxy)
	     if (bperz.eq.'X') tz=dmod(dmod(tz,boxz)+boxz,boxz)
	     if (asnew.eq.'') then
		write (fnum,10) aslst(i),tx-hx,ty-hy,tz-hz
	     else
		write (fnum,10) asnew,tx-hx,ty-hy,tz-hz
 	     endif
	  endif
	enddo
10	format (A2,3(2X,F10.4))

	return
	end
*-------------------------------------------------------------------------------
