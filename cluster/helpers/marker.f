*------------------------------------------------------------------------------
* Programm zum Markieren eines Atoms
*------------------------------------------------------------------------------

	program thinout

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'filep.inc'

	character*40 fout
	character*2 sym1,sym2
	integer i,j,step,maxi,nat,ind

	i=2
	call getarg(i,fout)
	i=3
	call getarg(i,sym1)
	i=4
	call getarg(i,sym2)

	call readdata

	hx=boxx/2
	hy=boxy/2
	hz=boxz/2

	open (1,file=fxyz)
	open (2,file=fout)

	maxi=ntim/ftxyz

	apsym(1)=sym1
	apsym(2)=sym2

	do i=1,maxi
	  call readXYZ(step,1)
	  do j=1,natom
	    api(j)=1
	  enddo
	  api(1)=2
	  call wrtXYZ(step,2)
	enddo

	close(1)

	end
