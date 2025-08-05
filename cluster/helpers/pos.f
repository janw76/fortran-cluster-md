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
	include 'boxpp.inc'

	character*40 fout
	integer i,j,step,maxi,nat
	double precision lx,ly,lz

	i=2
	call getarg(i,fout)

	call readdata

	open (1,file=fxyz)
	open (2,file=fout)

	maxi=ntim/ftxyz

	do i=1,maxi
	  lx=0.0d0
	  ly=0.0d0
	  lz=0.0d0
	  call readXYZ(step,1)
	  do j=2,natom
	    lx=lx+dcos(4.0d0*pi*(x(i)-x(1))/(atr*1.5d0))
	    ly=ly+dcos(4.0d0*pi*(y(i)-y(1))/(atr*1.5d0))
	    lz=lz+dcos(4.0d0*pi*(z(i)-z(1))/(atr*1.5d0))
	  enddo
	  lx=lx/dble(natom)
	  ly=ly/dble(natom)
	  lz=lz/dble(natom)
	  write (2,10) i,(lx+ly+lz)/3.0d0,lx,ly,lz
	enddo
10	format (I10,4(2X,F10.5))
	close(1)

	end
