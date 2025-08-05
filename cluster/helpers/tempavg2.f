*===============================================================================
*	HELPER PRORGAM temperature average
*-------------------------------------------------------------------------------
* file: tempavg.f
* date: 13-sep-2005 Jan Wedekind 
*-------------------------------------------------------------------------------
* this calculates the mean temperature of a cluster at a given size i.
* it can be applied on the sizer.f output files and takes ALL timesteps, when
* the cluster has the size i, takes the temperatures and averages!
* limeted to cluster sizes <=4096
*-------------------------------------------------------------------------------

	program tempavg
	
	implicit none

	integer i,ts,cs,ngas,startsize,endsize,refsize,fpt,ios,total
	double precision tc,tg,avgcs,endtime
	character*40 fsizer,fout,dummy
	character*80 readln
	integer count(4096)
	double precision tott(4096),mtemp(4096),tottsq(4096),std(4096)
	double precision gtot(4096),gmean(4096),gtotsq(4096),gstd(4096)
	double precision occ(4096)

*-------------------------------------------------------------------------------
* i		: loop variable
* ts		: time step
* cs		: cluster size
* ngas		: number of gas atoms
* tc		: kin. temperature of cluster
* tg		: kin. temperature of gas
* avgcs		: average cluster size
* std()		: standard deviation
* endtime 	: time to reach endsize for the first time
* refsize	: reference size for endtime
* fpt		: first passage-time through startsize
* total		: total number of cluster sampled
* occ		: 
*-------------------------------------------------------------------------------

	
	call getarg(1,fsizer)
	call getarg(2,fout)
	call getarg(3,dummy)
	   read (dummy,*) startsize
	call getarg(4,dummy)
	   read (dummy,*) endsize
	call getarg(5,dummy)
	   read (dummy,*) refsize

	open (1,file=fsizer)
	open (2,file=fout)

	cs=0
	fpt=0
		
	endtime=0
	do i=1,4096
	   mtemp(i) = 0
	   tott(i)  = 0
	   count(i) = 0
	   tottsq(i)= 0
	   std(i)   = 0
	   occ(i)   = 0
	   gtot(i)  = 0
	   gmean(i) = 0
	   gtotsq(i)= 0
	   gstd(i)  = 0
	enddo
	 
10	read(1,'(A)',IOSTAT=ios) readln
	if (ios.ne.0) goto 20

	read(readln,11) ts,cs,tc,ngas,tg,avgcs
11	format (I9,2X,2(2X,I6,2X,F7.2),2X,F10.3)

	tott(cs)   = tott(cs)   + tc
	tottsq(cs) = tottsq(cs) + tc*tc
	gtot(cs)   = gtot(cs)   + tg
	
	count(cs) = count(cs) + 1

	if (cs.ge.startsize) then
	   if (fpt.eq.0) fpt=ts
	endif
	if (cs.ge.refsize) then
	   if (endtime.eq.0) endtime=(real(ts)-real(fpt))/500000
	endif	

	goto 10

20	continue
	
	do i=1,4096
	   total = total + count(i)
	enddo

	do i=startsize,endsize
	   if (count(i).ne.0) then
	      mtemp(i) = tott(i)/DBLE(count(i))
	      std(i)   = tottsq(i)/DBLE(count(i)) - mtemp(i)
	      gmean(i) = gtot(i)/DBLE(count(i))
	      
	write(2,12) i,mtemp(i),std(i),gmean(i),endtime,count(i),total
12	format (I5,2X,4(F12.4,2X),I6,2X,I16)
	   endif
	enddo

	write (*,*) 'average temperatures written to: ',fout	
	
	close (1)	
	close (2)

	end

*-------------------------------------------------------------------------------
