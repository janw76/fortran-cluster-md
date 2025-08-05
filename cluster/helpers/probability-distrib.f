*===============================================================================
*	HELPER PRORGAM probability distribution
*-------------------------------------------------------------------------------
* file: probability-distrib.f
* date: 2007 Jan Wedekind 
*-------------------------------------------------------------------------------
* Takes the pn.*** output of sizer.f (dumped into one file
* dumped into one file by a script) and gets the total number 
*-------------------------------------------------------------------------------

	program probabilitydistrib
	
	implicit none

	integer i,cs,cn,ios
	double precision count(0:4096),SUM
	double precision PROB(0:4096)
	character*40 fin,fout
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

	call getarg(1,fin)
	if (fin.eq.'help') call helptext()
	call getarg(2,fout)
	
	do i = 0,4096
 	   count(i) = 0.0d0
	   PROB(i) = -1.0d0
	enddo	  
	
	open (1,file=fin)
	open (2,file=fout)
	
10	read(1,'(A)',IOSTAT=ios) readln
	if (ios.ne.0) goto 20
	
	read(readln,11) cs,cn
11	format (I5,2X,I12)

	count(cs) = count(cs) + DBLE(cn)
	goto 10
	
20	continue
	
	SUM = 0.0d0
	do i = 0,4096
	   SUM = SUM + count(i)
	enddo	

	do i = 0,4096
	 if (count(i).ne.0) Then
		PROB(i) = DBLE(count(i)) / DBLE(SUM)
		write (2,'(I5,2X,F15.9,2X,F16.0)') i,PROB(i),count(i)
	 endif
	enddo

	write (*,*) 'probability distribution written for: ',fin
	
	close (1)	

	end
	
	
*** SUBROUTINES
*** HELP ROUTINE *********************************************

	subroutine helptext()
	
	implicit none
	
	write (*,*) 'Probability distribution routine HELP'
	write (*,*) 'call: "probability-distrib infile outfile:'

	stop
	end

*-------------------------------------------------------------------------------
