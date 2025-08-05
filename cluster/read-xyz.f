*===============================================================================
*	READ XYZ-FILE
*===============================================================================
* This routine reads the positions of the particles from the XYZ-file.
* The box corners are autmatically skipped if present.
*-------------------------------------------------------------------------------

	subroutine readXYZ(step,fnum)

	implicit none
	
	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'boxpp.inc'
	include 'timep.inc'

	integer step,fnum
	integer i,nat
	double precision hx,hy,hz
	
*-------------------------------------------------------------------------------
* step		: time step number
* fnum		: file number for input
*...............................................................................
* i		: standard loop variable
* hx		: temporary variable for box center in x
* hy		: temporary variable for box center in y
* hz		: temporary variable for box center in z
*-------------------------------------------------------------------------------

	hx=boxx/2
	hy=boxy/2
	hz=boxz/2

	read (fnum,*) nat
	read (fnum,*) step

	if (nat.gt.natom) then
	  do i=1,nat-natom
	    read(fnum,10) ats,x(1),y(1),z(1)
	  enddo
	endif

	natom=min0(nat,natom)

	do i=1,natom
	  read(fnum,10) ats,x(i),y(i),z(i)
	  x(i)=x(i)+hx
	  y(i)=y(i)+hy
	  z(i)=z(i)+hz
	enddo
10	format (A2,3(2X,F10.4))

	return
	end
*-------------------------------------------------------------------------------
