*------------------------------------------------------------------------------
* Programm zur Berechnung des Druckes im System
*------------------------------------------------------------------------------

	program press

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'filep.inc'

	character*40 fout
	integer i,j,step,maxi,nat,ind

	i=2
	call getarg(i,fout)

	call readdata

	open (1,file=fxyz)
	open (2,file=fout)

	maxi=ntim/ftxyz

	do i=1,maxi
	  call readXYZ(step,1)
	  call AccWall()
	  call wrtAVE(step,2)
	enddo

	close(1)
	close(2)

	end
