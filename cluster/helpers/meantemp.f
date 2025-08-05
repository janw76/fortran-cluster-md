*===============================================================================
*	HELPER PRORGAM temperature average
*-------------------------------------------------------------------------------
* file: meantemp.f
* date: 13-sep-2005 Jan Wedekind 
*-------------------------------------------------------------------------------
* this calculates the mean temperature of a cluster at a given size i
* of ALL simulation runs.
* can be used with perl script to analyse tempavg output files.
*-------------------------------------------------------------------------------

	program meantemp
	
	implicit none

	integer i,ts,cs,startsize,endsize,num,ios,total
	double precision tc,err,endtime,occ
	double precision tg,errg
	character*40 fin,fout,dummy
	character*80 readln
	integer count(4096),k(4096)
	double precision tott(4096),mtemp(4096),toterr(4096),std(4096)
	double precision gtot(4096),gmean(4096),gtoterr(4096),gstd(4096)
	double precision totocc(4096)


*-------------------------------------------------------------------------------
* i		: loop variable
* ts		: time step
* cs		: cluster size
* ngas		: number of gas atoms
* tc		: kin. temperature of cluster
* tg		: kin. temperature of gas
* avgcs		: average cluster size
* std()		: standard deviation
* err()		: standard error
*-------------------------------------------------------------------------------

	
	call getarg(1,fin)
	call getarg(2,fout)
	call getarg(3,dummy)
	   read (dummy,*) startsize
	call getarg(4,dummy)
	   read (dummy,*) endsize

	open (1,file=fin)
	open (2,file=fout)
	
	do i=1,4096
	   mtemp(i) = 0
	   tott(i)  = 0
	   count(i) = 0
	   toterr(i)= 0
	   std(i)   = 0
	   totocc(i)= 0
	   k(i)     = 0
	   gtot(i)  = 0
	   gmean(i) = 0
	   gtoterr(i) = 0
	   gstd(i)  = 0
	enddo
	 
10	read(1,'(A)',IOSTAT=ios) readln
	
	if (ios.ne.0) goto 20
	
	read(readln,12) cs,tc,err,tg,endtime,num,total
12	format (I5,2X,4(F12.4,2X),I6,2X,I16)
	
	tott(cs)    = tott(cs)    + tc/err
	toterr(cs)  = toterr(cs)  + 1/err
	gtot(cs)    = gtot(cs)    + tg
	count(cs)   = count(cs)   + num
		
	k(cs) = k(cs) + 1

	goto 10

20	continue

	do i=startsize,endsize
	   if (k(i).ne.0) then
	      mtemp(i) = tott(i)/toterr(i)
	      std(i)   = SQRT(toterr(i)*k(i))
	      std(i)   = 1/std(i)
	      gmean(i) = gtot(i)/DBLE(k(i))

	   write(2,22) i,mtemp(i),std(i),gmean(i),k(i)
22	format (I5,2X,3(F12.4,2X),I12,2X,F12.4)
   	   write(*,22) i,mtemp(i),std(i),gmean(i),k(i)
	   endif
	enddo

	close (1)	
	close (2)
	write (*,*) 'Ok - ',fout
	end


*-----------------------------------------------------------------------
