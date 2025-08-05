*------------------------------------------------------------------------------
* Programm zum Ausduennen von Datenfiles
*------------------------------------------------------------------------------

	program thinout

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'filep.inc'

	character*40 fout,cskip
	integer i,j,step,maxi,skip,wrflag

	i=2
	call getarg(i,fout)
	i=3
	call getarg(i,cskip)
	read (cskip,*) skip

	call readdata

	open (1,file=fxyz)
	open (2,file=fout)

	maxi=ntim/ftxyz
	wrflag=0

	do i=1,maxi
	  call readXYZ(step,1)
	  if (wrflag.eq.0) then
	    wrflag=skip
	    call wrtXYZ(step,2)
	  endif
	  wrflag=wrflag-1
	enddo

	close(1)

	end
