	program mesher
	
	implicit none

	integer i,j,ts,atsum,cg,cn
	character*40 fout
	double precision tc,ek

	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'timep.inc'
	include 'filep.inc'
	include 'cvtpp.inc'

	i=2
	call getarg(i,fout)

	call readdata

	open (1,file=fcsi)
	open (2,file=fout)

*	write (2,*) 'RANK 2'
*	write (2,*) 'DIMENSIONS',200,1000
*	write (2,*) 'BOUNDS',0,1000.0,0,1000.0
*	write (2,*) 'NAME Cluster'
*	write (2,*) 'SCALAR'
*	write (2,*) 'DATA'

	write (2,*) 'NAME Cluster'
	write (2,*) 'RANK 2'
	write (2,*) 'DIMENSIONS 200 200'
	write (2,*) 'VECTOR 3'
	write (2,*) 'GRID'
	write (2,*) 'DATA'

	do i=1,200
	  do j=1,natom
	    csi(j)=0
	  enddo
	  read(1,*) ts
	  atsum=0
10	  read(1,11) cg,cn,ek,tc
11      format (2(3X,I7),3X,F17.15,3X,F10.5)
	  atsum=atsum+cn*cg
	  csi(cg)=cn
	  if (atsum.ne.natom) goto 10
	  do j=1,200
	    write (2,'(3(2X,I5),$)') i,j,min(250,csi(j))
*	    write (2,'(I5,$)') csi(j)
*	if (csi(j).gt.0) then
*	    write (2,'(1X,F8.5,$)') log(real(csi(j)))
*	else 
*	    write (2,'(1X,F8.5,$)') 0.0
*	endif
	  enddo
	  write (2,*)
	enddo

	close (1)
	close (2)
	end

