*------------------------------------------------------------------------------
* Programm zur 'Farbigen' Darstellung der Cluster
*------------------------------------------------------------------------------

	program selclu

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'timep.inc'
	include 'filep.inc'
	include 'cvtpp.inc'

	character*40 fout,fout2,cskip
	integer i,j,step,maxi,bcl,ind,skip,wrflag
	integer sum,lst
	integer apiorg(mxa),apilast(mxa)

	i=2
	call getarg(i,fout)
	i=3
	call getarg(i,fout2)
	i=4
	call getarg(i,cskip)
	read (cskip,*) skip

	call readdata

	open (1,file=fxyz)
	open (2,file=fout)
	open (3,file=fout2)

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
	    if (i.eq.1) then
	      do j=1,natom
	        apiorg(j)=api(j)
	        apilast(j)=api(j)
              enddo
	    endif

	    apsym(1)='S'
	    apsym(2)='N'
	    sum=0
	    lst=0
	    do j=1,natom
	      if (api(j).ne.apiorg(j).and.cl(j).ne.ind) sum=sum+1
	      if (api(j).ne.apilast(j).and.cl(j).ne.ind) lst=lst+1
	      apilast(j)=api(j)
	      api(j)=apiorg(j)
	    enddo
	    wrflag=skip
	    call wrtXYZ(step,2)
	    write (3,'(I9,3(2X,I5))') step,lst,sum,bcl
	  endif
	  wrflag=wrflag-1
	enddo

	close(1)
	close(2)

	end
