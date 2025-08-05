*===============================================================================
*	WRITE COORDINATE FILE
*===============================================================================
* file: wrt-xyz-m-vmd.f
*-------------------------------------------------------------------------------
* This routine writes a coordinate file that is readable by visualization
* programs. The routine works for both box and periodic boundary
* conditions and switches.
* This routine can use multiple atomic flavours.
* 
* 14:11 12.12.2008, JW
* - changed to write out monomners, Stillinger clusters and ten-Wolde frenkel
*
* changes by Jan Wedekind, 28/06/2005:
* when box corners are NOT shown the output file begins with
* N = number of atoms and
* step = time step
* thus, the file is readable by the VMD visualization software! 
*-------------------------------------------------------------------------------

	subroutine wrtXYZ(step,fnum)

	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'boxpp.inc'
	include 'timep.inc'
	include 'filep.inc'

	integer step,fnum
	integer i,k, count, frenkel, stillinger
	double precision hx,hy,hz,tx,ty,tz
	
*-------------------------------------------------------------------------------
* step		: time step number
* fnum		: file number for output
*...............................................................................
* i		: standard loop variable
* hx,hy,hz	: temporary variables for box center coordinates
* tx,ty,tz	: temporary variables coordinates
*-------------------------------------------------------------------------------

	open (66,file=fxave)
	open (67,file=fxrat)
	hx=boxx/2
	hy=boxy/2
	hz=boxz/2

	count = 0 
	frenkel = 0
	stillinger = 0

	if (step.lt.0) step=-step

	if (bcor.ne.'A'.and.step.ne.btim) then
	  write (fnum,*) natom
	  write (fnum,*) step
	elseif (bcor.eq.'A'.or.(bcor.eq.'X'.and.step.eq.btim)) then
	  write (fnum,*) natom+8
	  write (fnum,10) 'H ',-hx,-hy,-hz
	  write (fnum,10) 'H ',hx,-hy,-hz
	  write (fnum,10) 'H ',-hx,hy,-hz
	  write (fnum,10) 'H ',hx,hy,-hz
	  write (fnum,10) 'H ',-hx,-hy,hz
	  write (fnum,10) 'H ',hx,-hy,hz
	  write (fnum,10) 'H ',-hx,hy,hz
	  write (fnum,10) 'H ',hx,hy,hz
	elseif (step.eq.btim) then
	  write (fnum,*) natom
	  write (fnum,*) step
	endif

	write (66,*) natom
	write (66,*) step
	write (67,*) natom
	write (67,*) step
	
*** We now make three loops: 
*** 1.) All monomers (Stillinger 1, tWF 0)
*** 2.) All atoms of Stillinger clusters > 1
*** 3.) All liquid monomers acc. to tWF
	do i=1,natom
		tx=x(i)
		ty=y(i)
		tz=z(i)
		if (bperx.eq.'X') tx=dmod(dmod(tx,boxx)+boxx,boxx)
		if (bpery.eq.'X') ty=dmod(dmod(ty,boxy)+boxy,boxy)
		if (bperz.eq.'X') tz=dmod(dmod(tz,boxz)+boxz,boxz)
		write (fnum,10) apsym(api(i)),tx-hx,ty-hy,tz-hz
		count = count + 1
	enddo	
	
	do i=1,natom
		if (apsymi(i).eq.'S') then
			tx=x(i)
			ty=y(i)
			tz=z(i)
			if (bperx.eq.'X') tx=dmod(dmod(tx,boxx)+boxx,boxx)
			if (bpery.eq.'X') ty=dmod(dmod(ty,boxy)+boxy,boxy)
			if (bperz.eq.'X') tz=dmod(dmod(tz,boxz)+boxz,boxz)
			write (66,10) 'S ',tx-hx,ty-hy,tz-hz
			stillinger = stillinger + 1
			k = i
		endif
	enddo
	if (stillinger.gt.0) then
		do i = 1,natom-stillinger
			tx=x(k)
			ty=y(k)
			tz=z(k)
			if (bperx.eq.'X') tx=dmod(dmod(tx,boxx)+boxx,boxx)
			if (bpery.eq.'X') ty=dmod(dmod(ty,boxy)+boxy,boxy)
			if (bperz.eq.'X') tz=dmod(dmod(tz,boxz)+boxz,boxz)
			write (66,10) 'S ',tx-hx,ty-hy,tz-hz
		enddo
	elseif (stillinger.eq.0) then
		do i=1,natom
			write (66,10) 'S ',4*hx,4*hy,4*hz
		enddo
	endif
		
	do i=1,natom
		if (apsymi(i).eq.'F') then
			tx=x(i)
			ty=y(i)
			tz=z(i)
			if (bperx.eq.'X') tx=dmod(dmod(tx,boxx)+boxx,boxx)
			if (bpery.eq.'X') ty=dmod(dmod(ty,boxy)+boxy,boxy)
			if (bperz.eq.'X') tz=dmod(dmod(tz,boxz)+boxz,boxz)
			write (67,10) 'F ',tx-hx,ty-hy,tz-hz
			frenkel = frenkel + 1
			k = i
		endif
	enddo
	if (frenkel.gt.0) then
		do i = 1,natom-frenkel
			tx=x(k)
			ty=y(k)
			tz=z(k)
			if (bperx.eq.'X') tx=dmod(dmod(tx,boxx)+boxx,boxx)
			if (bpery.eq.'X') ty=dmod(dmod(ty,boxy)+boxy,boxy)
			if (bperz.eq.'X') tz=dmod(dmod(tz,boxz)+boxz,boxz)
			write (67,10) 'F ',tx-hx,ty-hy,tz-hz
		enddo
	elseif (frenkel.eq.0) then
		do i=1,natom
			write (67,10) 'F ',4*hx,4*hy,4*hz
		enddo
	endif
	

	
10	format (A2,3(2X,F10.4))

*	write (*,*) 
*	write (*,*) 'Count:', count

	return
	end
*-------------------------------------------------------------------------------
