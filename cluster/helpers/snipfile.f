*------------------------------------------------------------------------------
* Programm zum Ausschneiden von Teilen von Datenfiles
*------------------------------------------------------------------------------

	program thinout

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'filep.inc'

	character*40 fout,cskip
	integer i,step,maxi,stconf,econf

	i=2
	call getarg(i,fout)
	i=3
	call getarg(i,cskip)
	read (cskip,*) stconf
	i=4
	call getarg(i,cskip)
	read (cskip,*) econf

	call readdata

	open (1,file=fxyz)
	open (2,file=fout)

	maxi=ntim/ftxyz

	do i=1,maxi
	  call readXYZ(step,1)
	  if (i.ge.stconf.and.i.le.econf) then
	    call wrtXYZ(step,2)
	  endif
	enddo

	close(1)

	end
