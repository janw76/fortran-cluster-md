*===============================================================================
*	HELPER PRORGAM temperature histograms
*-------------------------------------------------------------------------------
* file: t-hist.f
* date: 13-sep-2005 Jan Wedekind 
*-------------------------------------------------------------------------------


	program thist
	
	implicit none

	integer i,ts,cs,ngas,ios,sum
	double precision tc,tg,avgcs,norm,P
	character*40 fsizer,fout
	character*80 readln
	integer bin(4096),max

	
	call getarg(1,fsizer)
	call getarg(2,fout)

	open (1,file=fsizer)
	open (2,file=fout)

	do i=1,4096
	   bin(i) = 0	  
	enddo
	 
10	read(1,'(A)',IOSTAT=ios) readln
	if (ios.ne.0) goto 20

	read(readln,11) ts,cs,tc,ngas,tg,avgcs
11	format (I9,2X,2(2X,I6,2X,F7.2),2X,F10.3)

	bin(int(tc)) = bin(int(tc)) + 1

	goto 10

20	continue
	
	max = 0
	sum = 0
	do i=1,300
	   if (bin(i).gt.max) max = bin(i)
	   sum = sum + bin(i)
	enddo

	do i=1,300
	   norm = REAL(bin(i))/REAL(max)
	   P = REAL(bin(i))/REAL(sum)
	   write (2,'(I4,2X,I16,2X,F18.8,2X,F18.8)') i,bin(i),norm,P
	enddo

	write (*,*) 'temperature histogram written to: ',fout	
	
	close (1)	
	close (2)

	end

*-------------------------------------------------------------------------------
