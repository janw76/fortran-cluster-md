*===============================================================================
*	HELPER PRORGAM diffusivity.f
*-------------------------------------------------------------------------------
* file: diffusivity.f
* date: 06-sep-2005 Jan Wedekind 
*-------------------------------------------------------------------------------
* uses the output file of SIZER to get the the "diffusivity" of the critical
* size from the top of the barrier. n* must be known before.
*-------------------------------------------------------------------------------

	program diffusivity
	
	implicit none

	integer i,ts,cs,ngas,startsize,samples,ios
	integer tnull,tnew,ntemp,SWITCH,count(10000),totcount(10000)
	double precision tc,tg,avgcs,dnsq,avgdn(10000)
	character*40 fsizer,fout,dummy
	character*80 readln
	

*-------------------------------------------------------------------------------
* i		: loop variable
* ts		: time step
* cs		: cluster size
* ngas		: number of gas atoms
* tc		: kin. temperature of cluster
* tg		: kin. temperature of gas
* avgcs		: average cluster size
*-------------------------------------------------------------------------------

	
	call getarg(1,fsizer)
	call getarg(2,fout)
	call getarg(3,dummy)
	   read (dummy,*) startsize
	call getarg(4,dummy)
	   read (dummy,*) samples
	call getarg(5,dummy)
	   read (dummy,*) SWITCH

	open (1,file=fsizer)
	open (2,file=fout)

	if (SWITCH.eq.1) THEN

	cs = 0
	i = 0
	dnsq = 0.0d0
	tnull = 0

10	read(1,'(A)',IOSTAT=ios) readln
 	if (ios.ne.0) goto 70
	
	read(readln,11) ts,cs,tc,ngas,tg,avgcs
11	format (I9,2X,2(2X,I6,2X,F7.2),2X,F10.3)

	if (cs.eq.startsize) then
	   
	  tnull = ts
*** 	  begin do in if	
	  do i = 1,samples
	   read(1,'(A)',IOSTAT=ios) readln
	   if (ios.ne.0) goto 70
	   read(readln,11) ts,cs,tc,ngas,tg,avgcs	   
	   tnew = ts - tnull
	   ntemp = (cs-startsize) * (cs-startsize)
	   dnsq = dble(ntemp)
12	   format (I9,2X,F10.3)	
	   write (2,12) tnew/1000,dnsq
	  enddo	
	endif
	goto 10

70	close (1)	
	close (2)

*********************************
	ELSEIF (SWITCH.EQ.2) THEN
	
	do i = 1,samples
	   avgdn(i) = 0
	enddo

	count = 0

30	continue
	read(1,'(A)',IOSTAT=ios) readln
 	if (ios.ne.0) goto 40
	read(readln,12) tnew,dnsq
		
	avgdn(tnew) = avgdn(tnew) + dnsq
	count(tnew) = count(tnew) + 1

	goto 30

40	continue

	do i = 1,samples
	   avgdn(i) = avgdn(i) / dble(count(i))
	   write (2,12) i,avgdn(i)
	enddo

20	write (*,*) 'diffusivity taken for: ',fsizer
	close (1)	
	close (2)

***
	ELSEIF (SWITCH.EQ.3) THEN
***

	do i = 1,10
	   avgdn(i) = 0
	   count(i) = 0
	   totcount(i) = 0
	enddo
	
	samples = 0
	
50	continue
	read(1,'(A)',IOSTAT=ios) readln
 	if (ios.ne.0) goto 60
	read(readln,13) cs,tc,i
13	format (I3,2X,F7.4,2X,I15) 	
	
	avgdn(cs) = avgdn(cs) + tc
	count(cs) = count(cs) + 1
	totcount(cs) = totcount(cs) + i

	goto 50

60	continue

	do i = 2,8
	   write (2,13) i,avgdn(i)/dble(count(i)),totcount(i)
	enddo
	write (*,*) 'full meantemp written for: ',fsizer
	close (1)	
	close (2)	
	
	ENDIF
	end

*-------------------------------------------------------------------------------
