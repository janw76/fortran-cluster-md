*------------------------------------------------------------------------------
* Programm zum Trennen von xyz-Datenfiles
*------------------------------------------------------------------------------

	program xyzsplitter

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'filep.inc'

	character*40 outfile
	integer i,j,step

	call readdata

	open (1,file='vis.xyz')

	do j=1,5
	  outfile='vis-'//char(j+96)//'.xyz'
	  open (2,file=outfile)
	  do i=1,1000
	    call readXYZ(step,1)
	    call wrtXYZ(step,2)
	  enddo
	  close(2)
	enddo

	close(1)

	end
