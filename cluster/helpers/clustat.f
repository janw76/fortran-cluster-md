*===============================================================================*       BAR-HISTOGRAM
*-------------------------------------------------------------------------------* file: bars.f
*-------------------------------------------------------------------------------
* call: clustat <3d-file> <outfile>
*-------------------------------------------------------------------------------

	program clustat
	
	implicit none

	include 'const.inc'
	include 'atomp.inc'
	include 'timep.inc'
	include 'filep.inc'
	include 'cvtpp.inc'

	integer i,j,ts,atsum,cg,cn
	integer szlo,szhi
	character*40 fout
	double precision tc,ek,db
	double precision h

	i=2
	call getarg(i,fout)

	call readdata()

	open (1,file=fcsi)
	open (2,file=fout)

	do i=1,natom
	  csi(i)=0
	enddo

	do i=1,ntim/ftcsi
	  read(1,*) ts
	  atsum=0
10	  read(1,11) cg,cn,ek,tc
	  if (i.ge.500) then
	    csi(cg)=csi(cg)+cn
	  endif
	  atsum=atsum+cn*cg
	  if (atsum.ne.natom) goto 10
	enddo
	do i=1,natom
	  write (2,21) i,real(csi(i))/500.0,real(csi(i))
	enddo

11      format (2(3X,I7),3X,F17.15,3X,F10.5)
21	format (I5,2(2X,F10.4))

	close (1)
	close (2)
	end

