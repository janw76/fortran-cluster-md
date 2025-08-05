	program sizer
	
	implicit none

	integer i,ts,atsum,cg,cn,cd
	character*40 fout
	double precision tc,ek,eks

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

	do i=1,ntim/ftcsi
	  read(1,*) ts
	  atsum=0
	  eks=0.0d0
10	  read(1,11) cg,cn,ek,tc
11      format (2(3X,I7),3X,F17.15,3X,F10.5)
	  atsum=atsum+cn*cg
	  eks=eks+ek*cn*cg
	  if (atsum.ne.natom) goto 10
	  eks=eks-ek*cn*cg
	  cd=cd+cg
	  if (cg.eq.natom) then
	    write (2,21) ts,cg,tc,natom-cg,0,real(cd)/i
	  else
	    write (2,21) ts,cg,tc,natom-cg,
     $      eks/((natom-cg)*1.5d0*cikb),real(cd)/i
	  endif
21	format (I9,2X,2(2X,I6,2X,F7.2),2X,F10.3)
	enddo

	write (*,*) 'Average Cluster Size:',real(cd)/real(ntim/ftcsi)

	close (1)
	close (2)
	end

