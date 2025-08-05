*===============================================================================
*	HELPER PRORGAM mean passage time
*-------------------------------------------------------------------------------
* file: mpt.f
* date: 13-sep-2005 Jan Wedekind 
*-------------------------------------------------------------------------------
* this calculates the mean passage times from any dump file of passage-time.f
* can be used in a batch script. Up to cluster sizes of 4096.
* CAREFUL: times are in ns assuming a timestep of 2 fs!!!
*-------------------------------------------------------------------------------

	program meanpt
	
	implicit none

	integer i,ts,cs,ngas,startsize,endsize,ios
	double precision tc,tg,avgcs
	character*40 fsizer,fout,dummy
	character*80 readln
	double precision pt(4096),ptsq(4096),temp(4096)
	double precision tempsq(4096),cgs,acgs(4096)
	double precision mpt(4096),mtemp(4096),gtemp(4096),gtempsq(4096)
	double precision stda(4096),stdb(4096),stdc(4096),mgt(4096)
	integer count(4096)

*-------------------------------------------------------------------------------
* i		: loop variable
* ts		: time step
* cs		: cluster size
* ngas		: number of gas atoms
* tc		: kin. temperature of cluster
* tg		: kin. temperature of gas
* avgcs		: average cluster size
* pt		: passage-times
* temp:		: temperatures
* mtemp:	: mean temperature
* _sq		: squared value (for std. dev.)
* count(i)	: counting numbers of cluster of size i
* stda/b/c	: standard errors
*-------------------------------------------------------------------------------
	
	call getarg(1,fsizer)
	call getarg(2,fout)
	call getarg(3,dummy)
	   read (dummy,*) startsize
	call getarg(4,dummy)
	   read (dummy,*) endsize

	open (1,file=fsizer)
	open (2,file=fout)

	cs=0

	do i=1,4096
	   pt(i)=0
	   ptsq(i)=0
	   temp(i)=0
	   tempsq(i)=0
	   gtemp(i)=0
	   gtempsq(i)=0
	   acgs = 0
	   count(i)=0
	   mpt(i)=0
	   mtemp(i)=0
	   mgt(i)=0
	   stda(i)=0
	   stdb(i)=0
	   stdc(i)=0
	enddo	

	open (1,file=fsizer)

10	read(1,'(A)',IOSTAT=ios) readln

	if (ios.ne.0) goto 20

	read(readln,11) ts,cs,tc,ngas,tg,avgcs,cgs
11	format (I9,2X,2(2X,I6,2X,F7.2),2X,F10.3,F7.3)
	
	pt(cs)=pt(cs)+DBLE(ts)/500000
	ptsq(cs)=ptsq(cs)+(DBLE(ts)/500000)*(DBLE(ts)/500000)
	temp(cs)=temp(cs)+tc
	tempsq(cs)=tempsq(cs)+tc*tc
	gtemp(cs)=gtemp(cs)+tg
	gtempsq(cs)=gtempsq(cs)+tg*tg
	acgs(cs) = acgs(cs) + cgs
	count(cs)=count(cs)+1
	
	goto 10

20	continue
	
	do i=startsize,endsize

	   mpt(i)=(pt(i)/count(i))
	   stda(i)=ptsq(i)/count(i)-mpt(i)*mpt(i)
	   stda(i)=SQRT(stda(i)/count(i))
	   mtemp(i)=temp(i)/count(i)
	   stdb(i)=tempsq(i)/count(i)-mtemp(i)*mtemp(i)
	   stdb(i)=SQRT(stdb(i)/count(i))
	   mgt(i)=gtemp(i)/count(i)
	   stdc(i)=gtempsq(i)/count(i)-mgt(i)*mgt(i)
	   stdc(i)=SQRT(stdc(i)/count(i))
	   acgs(i)=acgs(i)/count(i)

	write (2,12) mpt(i),stda(i),i,mtemp(i),stdb(i),mgt(i),
     $	    stdc(i),acgs(i),count(i)
12	format (F15.4,2X,F12.6,2X,I6,2X,5(F9.2,2X),I6)
	write (*,12) mpt(i),stda(i),i,mtemp(i),stdb(i),mgt(i),
     $	    stdc(i),acgs(i),count(i)

	enddo

	write (*,*) 'mean-passage-times calculated and written to: ',fout	

	close (1)
	close (2)

	end

*-------------------------------------------------------------------------------
