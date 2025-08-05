*------------------------------------------------------------------------------
* Programm zur 'Farbigen' Darstelung des Groessten Clusters im System
*------------------------------------------------------------------------------
* file: bigcl.f
* date: 05-feb-2001
*------------------------------------------------------------------------------

	program bigcl

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'timep.inc'
	include 'filep.inc'
	include 'cvtpp.inc'

	character*40 fout,cskip
	integer i,j,step,maxi,bcl,ind,skip,wrflag

	i=2
	call getarg(i,fout)
	i=3
	call getarg(i,cskip)
	read (cskip,*) skip

	call readdata

	open (1,file=fxyz)
	open (2,file=fout)

	maxi=ntim/ftxyz
	wrflag=0

	do i=1,maxi
	  write (*,*) 'Reading config:',i
	  call readXYZ(step,1)
	  if (wrflag.eq.0) then
	    call stoddard()
	    bcl=0
	    do j=1,natom
	      if (cs(j).gt.bcl) then
	        bcl=cs(j)
	        ind=j
	      endif
	    enddo
	    do j=1,natom
	      if (cl(j).eq.ind) then
	        api(j)=2
	      else
	        api(j)=1
	      endif
	    enddo

	    apsym(1)='S'
	    apsym(2)='N'
	    wrflag=skip
	    call wrtXYZ(step,2)
	  endif
	  wrflag=wrflag-1
	enddo

	close(1)
	close(2)

	end
