*===============================================================================
*	HELPER PRORGAM temperature average
*-------------------------------------------------------------------------------
* file: tempavg.f
* date: 13-sep-2005 Jan Wedekind 
*-------------------------------------------------------------------------------
* this calculates the mean temperature of a cluster at a given size i.
* it can be applied on the sizer.f output files and takes ALL timesteps, when
* the cluster has the size i, takes the temperatures and averages!
* limeted to cluster sizes <=1000
*-------------------------------------------------------------------------------

	program tempavg
	
	implicit none

	integer i,ts,cs,ngas,startsize,endsize,lines,io
	double precision tc,tg,avgcs
	character*40 fsizer,fout,dummy
	character*80 readln
	double precision tct(1000),mtemp(1000),temp(1000),clcount(1000)
	double precision stda(100),stdb(100)
	

*-------------------------------------------------------------------------------
* i		: loop variables
* ts		: time step
* cs		: cluster size
* ngas		: number of gas atoms
* tc		: kin. temperature of cluster
* tg		: kin. temperature of gas
* tct		: total sum of temperatures of cluster of size (i)
* mtemp()	: mean temperature of cluster of size (i)
* clcount()	: number of clusters of size (i)
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
	lines=0

	do i=1,1000
	   tct(i)=0
	   mtemp(i)=0
	   temp(i)=0
	   clcount(i)=0
	   stda(i)=0
	   stdb(i)=0
	enddo	

10	read(1,'(A)',IOSTAT=io) readln
	if (io.ne.0) then
	   goto 20
	else
	   read(readln,11) ts,cs,tc,ngas,tg,avgcs
11	format (I9,2X,2(2X,I6,2X,F7.2),2X,F10.3)
	   tct(cs)=tct(cs)+tc
	   clcount(cs)=clcount(cs)+1
	   goto 10
	endif

20	continue	   	

	do i=startsize,endsize
	   mtemp(i)=tct(i)/clcount(i)
	   write (2,12) i,mtemp(i),clcount(i) 
12	format (I5,2X,F4.4)
	enddo

	write (*,*) 'average cluster temperatures written to: ',fout	
	write (2,11) -1,0,0,0,0,0
	
	close (1)	
	close (2)

	end

*-------------------------------------------------------------------------------
