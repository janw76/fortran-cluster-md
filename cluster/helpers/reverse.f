*------------------------------------------------------------------------------
* Programm zum invertieren eines xyz-files
*------------------------------------------------------------------------------

	program reverse

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'filep.inc'

	integer mxs
	parameter (mxs=2000)
	character*40 fout
	integer i,j,maxi,step
	double precision xx(mxs,mxa),yy(mxs,mxa),zz(mxs,mxa)
	integer nat(mxs)

	i=2
	call getarg(i,fout)

	call readdata

	open (1,file=fxyz)
	open (2,file=fout)

	maxi=ntim/ftxyz
	if (maxi.gt.mxs) then
	  write(*,*) 'Too many configurations!'
	  write(*,*) 'To be read:',maxi
	  write(*,*) 'Maximum   :',mxs
	  stop
	endif

	do i=1,mxs
	  do j=1,mxa
	    xx(i,j)=0.0d0
	    yy(i,j)=0.0d0
	    zz(i,j)=0.0d0
	    nat(i)=0.0d0
	  enddo
	enddo

	do i=1,maxi
	  write (*,*) 'Reading step ',i,' of ',maxi
	  call readXYZ(step,1)
	  nat(i)=natom
	  do j=1,natom
	    xx(i,j)=x(j)
	    yy(i,j)=y(j)
	    zz(i,j)=z(j)
	  enddo
	enddo

	do i=maxi,1,-1
	  write (*,*) 'Writing step ',i,' of ',maxi
	  natom=nat(i) 
	  do j=1,mxa
	    x(j)=xx(i,j)
	    y(j)=yy(i,j)
	    z(j)=zz(i,j)
	  enddo
	  call wrtXYZ(maxi-i,2)
	enddo

	close(1)
	close(2)

	end
