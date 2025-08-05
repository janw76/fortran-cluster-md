*------------------------------------------------------------------------------
* Programm zur Ausgabe einer Trajektorie eines Atoms
*------------------------------------------------------------------------------

	program thinout

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'filep.inc'

	character*40 fout
	integer i,j,step,maxi,an

	i=3
	call getarg(i,fout)
	read (fout,*) an
	i=2
	call getarg(i,fout)

	call readdata

	if(an.lt.1.or.an.gt.natom) an=1

	open (1,file=fxyz)
	open (2,file=fvel)
	open (3,file=fout)

	maxi=ntim/ftxyz

	do i=1,maxi
	  call readXYZ(step,1)
	  call readVEL(step,2)
	  write (3,10) step,x(an),y(an),z(an),vx(an),vy(an),vz(an)
	enddo

10	format (I10,3(2X,F10.4),3(3X,F15.10))

	close(1)
	close(2)
	close(3)

	end
