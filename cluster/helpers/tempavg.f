*===============================================================================
*	HELPER PRORGAM temperature average
*-------------------------------------------------------------------------------
* file: tempavg.f
* date: 13-sep-2005 Jan Wedekind 
*-------------------------------------------------------------------------------
* it's actually in order to get the histograms, so the peak temperatures
* it just writes out the info for one specific size n
*-------------------------------------------------------------------------------

	program tempavg
	
	implicit none

	integer i,j,ts,cs,ngas,startsize,endsize,ios
	double precision tc,tg,avgcs
	character*40 fsizer,dummy
	character*80 readln
	character*40 filename,ext


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
	if (fsizer.eq.'help') call helptext()
	call getarg(2,dummy)
	   read (dummy,*) startsize
	call getarg(3,dummy)
	   read (dummy,*) endsize

	ext = '.nnn'
	do i = startsize,endsize
 	   write(filename,'(I4,A4)') i,ext
	   j = i + 10
	   open(j,file=filename,access='append')
	enddo	  

	open (1,file=fsizer)

	cs=0
	
10	read(1,'(A)',IOSTAT=ios) readln
	if (ios.ne.0) goto 20

	read(readln,11) ts,cs,tc,ngas,tg,avgcs
11	format (I9,2X,2(2X,I6,2X,F7.2),2X,F10.3)

	if ((cs.ge.startsize).and.(cs.le.endsize)) then
	   i = cs + 10
	   write(i,'(A)') readln
	endif

	goto 10
	
20	continue

	do i=startsize,endsize
	   j = i + 10
 	   close(j)
	enddo

	write (*,*) 'values written for: ',fsizer	
	
	close (1)	

	end

*-------------------------------------------------------------------------------
*** SUBROUTINES
*** HELP ROUTINE *********************************************

	subroutine helptext()
	
	implicit none
	
	write (*,*) 'Here be HELP!'
	
	stop
	end

*-------------------------------------------------------------------------------