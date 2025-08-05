*===============================================================================
*	HELPER PRORGAM passage-time.f
*-------------------------------------------------------------------------------
* file: passage-time.f
* date: 06-sep-2005 Jan Wedekind 
*-------------------------------------------------------------------------------
* uses the output file of SIZER to get the timestep when the largest cluster
* passes through a specific size for the first time.
* the mean passage time of many individual simulations at same conditions can
* then be used to evaluate the rate, barrier height and critical cluster size
*-------------------------------------------------------------------------------
* 12:50 26.02.2008: probab gets the probability distribution of the LARGEST

	program passagetime
	
	implicit none

	integer i,j,ts,cs,ngas,startsize,endsize,warmup,ios,flag
	integer count(0:4096)
	double precision tc,tg,avgcs,cgs
	character*40 fsizer,fout,dummy,probab
	character*80 readln	
    
*-------------------------------------------------------------------------------
* i,j		: loop variable
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
	   read (dummy,*) endsize
	call getarg(5,dummy)
	   read (dummy,*) warmup
 
	probab = 'probab.largest'
	open (1,file=fsizer)
	open (2,file=fout)
	open (3,file=probab)
	do i = 0,4096
	   count(i) = 0
	enddo	  
	cs=0
	i=startsize
	flag = -1
	
	
10	read(1,'(A)',IOSTAT=ios) readln
 
	if (ios.ne.0) goto 20
	
	read(readln,11) ts,cs,tc,ngas,tg,avgcs,cgs
11	format (I9,2X,2(2X,I6,2X,F7.2),2X,F10.3,F7.3)

	if (flag.ne.1) Then 
	   if (cs.ge.endsize) flag = 1
	   if (ts.ge.warmup) count(cs) = count(cs) + 1
	endif
    
	if (cs.ge.i) then
	   do j = i,cs 
	      write (2,11) ts,j,tc,ngas,tg,avgcs,cgs
	   enddo
	   i = cs + 1
	   goto 10
	else
	   goto 10
	endif

20	if (i.lt.endsize) then
	   write (*,*)
	write (*,*) '----->    !!! Cluster did NOT reach size size/in',endsize,
     $	    '/',fsizer,'after',ts*2/1000000,'ns'
	   write (*,*)
	else
	   write (*,*) 'passage times written to: ',fout
	endif
	
	do i = 0,endsize
	  if (count(i).ne.0) then
	    write (3,'(I5,2X,I12)') i,count(i)  
	  endif
	enddo
	
	close (1)	
	close (2)
	close (3)
	end

*-------------------------------------------------------------------------------
