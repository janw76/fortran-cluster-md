	program sizer
	
	implicit none

	integer i,ts,atsum,cg,cn,cd,ocn,ocg
	character*40 fout
	character*80 readln
	double precision tc,ek,eks,oek,otc

	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'timep.inc'
	include 'filep.inc'

	i=2
	call getarg(i,fout)

	call readdata()

	open (1,file=fcsi)
	open (2,file=fout)

	cd=0

	read(1,'(A)') readln
	do i=1,ntim/ftcsi
	  read (readln,*) ts
	  atsum=0
	  eks=0.0d0
10	  read(1,'(A)') readln
	  read(readln,11) cg,cn,ek,tc
11      format (2(3X,I7),3X,F17.15,3X,F10.5)	
	  if (cn.ne.0) then
	    atsum=atsum+cn*cg
	    eks=eks+ek*cn*cg
	    ocn=cn
	    ocg=cg
	    oek=ek
	    otc=tc
	    goto 10 
	  else
	    eks=eks-oek*ocn*ocg
	    cd=cd+ocg
	    if (ocg.eq.atsum) then
	      write (2,21) ts,ocg,otc,atsum-ocg,0,real(cd)/i
	    else
	      write (2,21) ts,ocg,otc,atsum-ocg,
     $        eks/((atsum-ocg)*1.5d0*cikb),real(cd)/i
	    endif
	  endif
21	format (I9,2X,2(2X,I6,2X,F7.2),2X,F10.3)
	enddo

	write (*,*) 'Average Cluster Size:',real(cd)/real(ntim/ftcsi)

	close (1)
	close (2)
	end

