*===============================================================================*       BAR-HISTOGRAM
*-------------------------------------------------------------------------------* file: bars.f
*-------------------------------------------------------------------------------
* call: bars <3d-file> <start-size> <endsize> <outfile>
*-------------------------------------------------------------------------------
* This program analyzes the csi-file of a simulation and produces a file for
* visualization by SciAn (NFF format). A bar-histogram is generated which
* shows the cluster size distribution vs. time.
*-------------------------------------------------------------------------------
	program bars
	
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
	read(fout,*) szlo

	i=3
	call getarg(i,fout)
	read(fout,*) szhi

	i=4
	call getarg(i,fout)

	call readdata()

	open (1,file=fcsi)
	open (2,file=fout)

	db=0.25/2.0

	write (2,*) 'f 1.0 1.0 1.0 0 0 0 0 0'
	do i=1,ntim/ftcsi
	  read(1,*) ts
	  write (2,*) 't ',i-1
	  atsum=0
10	  read(1,11) cg,cn,ek,tc
	  csi(cg)=cn
	  atsum=atsum+cn*cg
	  if (atsum.ne.natom) goto 10
	  do j=szlo,szhi
	    h=float(csi(j))
*	    if (csi(j).gt.0) then
*	      h=log10(float(csi(j)))
*	    else
*	      h=0.0
*	    endif
	    write (2,*) 'p 4'
	    write (2,21) 30.0*(-db+j)+30,0.0,0.0
	    write (2,21) 30.0*(+db+j)+30,0.0,0.0
	    write (2,21) 30.0*(+db+j)+30,0.0+h,0.0
	    write (2,21) 30.0*(-db+j)+30,0.0+h,0.0
	  enddo
	  write (2,*)
	enddo
11      format (2(3X,I7),3X,F17.15,3X,F10.5)
21	format (3(2X,F10.4))

	close (1)
	close (2)
	end

