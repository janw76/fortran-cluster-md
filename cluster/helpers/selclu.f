*------------------------------------------------------------------------------
* Programm zur 'Farbigen' Darstellung der Cluster
*------------------------------------------------------------------------------

	program selclu

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'filep.inc'
	include 'cvtpp.inc'

	character*40 fout,fout2,cskip
	integer i,j,step,maxi,skip,wrflag,mflag

	i=2
	call getarg(i,fout)
	i=3
	call getarg(i,fout2)
	i=4
	call getarg(i,cskip)
	read (cskip,*) skip
	i=5
	call getarg(i,cskip)
	read(cskip,*) mflag
	if (mflag.eq.0) mflag=-1	

	call readdata

	open (1,file=fxyz)
	open (2,file=fout)
	open (3,file=fout2)

	maxi=ntim/ftxyz
	wrflag=0

	do i=1,maxi
	  call readXYZ(step,1)
	  if (wrflag.eq.0) then
	    wrflag=skip
	    call wrtXYZ(step*mflag,2)
	    write (3,*) natom
	    write (3,*) step
	    do j=1,natom
	      write (3,'(A2,2X,I4,2X,F7.3,2X,I4)') 
     $	            'H ',j,sqrt(real(csi(j))),0
	    enddo
	  endif
	  wrflag=wrflag-1
	  write (*,*) 'Step ',i,' of ',maxi
	enddo

	close(1)
	close(2)

	end
