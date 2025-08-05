	program sizer
	
*********************************************************************	
* Extensively rewritten by Jan Wedekind to meet
* the requirements of new subroutines
* 04/10/2007 (JW)
*   - added subroutine for probability counting
*   - added help subroutine
*  
*********************************************************************

	implicit none

	integer i,ts,atsum,cg,cn,cd,ocn,ocg,ios,warmup
	integer count(0:4096),nmax,flag
	character*40 fout,fpn,dummy
	character*80 readln
	double precision tc,ek,eks,oek,otc
	double precision avgcgs,oavgcgs

	integer SWITCH

	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'timep.inc'
	include 'filep.inc'

	call getarg(1,dummy)
	if (dummy.eq.'help') call helptext()
	call getarg(2,fout)
	call getarg(3,dummy)
	   read (dummy,*) SWITCH
	call getarg(4,dummy)
	   read (dummy,*) nmax
	call getarg(5,dummy)
	   read (dummy,*) warmup

	call readdata()

	if (SWITCH.eq.2) fcsi = 'fsi'
	write (*,*) 'reading: ',fcsi
	
	open (1,file=fcsi)
	open (2,file=fout)

	do i = 0,4096
	   count(i) = 0
	enddo	  
	
	flag = -1
	
	cd = 0
10	format (I9,2X,2(2X,I6,2X,F7.2),2X,F10.3,F7.3)

*** SWITCH = 1
*** Standard sizer for Stillinger definition

	if (SWITCH.eq.1) THEN
	
	fpn = 'pn.stillinger'
	open (3,file=fpn)
	
	read(1,'(A)',IOSTAT=ios) readln
	if (ios.ne.0) goto 30
		
31	continue

	read (readln,*,IOSTAT=ios) ts
	  atsum=0
	  eks=0.0d0
21	  read(1,'(A)',IOSTAT=ios) readln
	  if (ios.ne.0) goto 30
	  read(readln,11) cg,cn,avgcgs,ek,tc
11      format (2(3X,I7),3X,F7.3,3X,F17.15,3X,F10.5)
	
	if (cn.ne.0) then
	    atsum=atsum+cn*cg
	    eks=eks+ek*cn*cg
	    ocn=cn
	    ocg=cg
	    oavgcgs = avgcgs
	    oek=ek
	    otc=tc
	    
	    if (flag.ne.1) Then 
	      if (cg.ge.nmax) flag = 1
	      if (ts.ge.warmup) count(cg) = count(cg) + cn
	    endif	    
	    
	    goto 21
	else
	    eks=eks-oek*ocn*ocg
	    cd=cd+ocg
	    if (ocg.eq.atsum) then
	      write (2,10) ts,ocg,otc,atsum-ocg,0,real(cd)/i,oavgcgs	      
	    else
	      write (2,10) ts,ocg,otc,atsum-ocg,
     $        eks/((atsum-ocg)*1.5d0*cikb),real(cd)/i,oavgcgs
	    endif
	endif

	goto 31


*** SWITCH = 2
*** Sizer for Frenkel file (SWITCH = 2)

	ELSEIF (SWITCH.eq.2) THEN
	
	fpn = 'pn.twf'
	open (3,file=fpn)
	
	read(1,'(A)',IOSTAT=ios) readln
	if (ios.ne.0) goto 30
		
32	continue

	read (readln,*,IOSTAT=ios) ts
	  atsum=0
	  eks=0.0d0
22	  read(1,'(A)',IOSTAT=ios) readln
	  if (ios.ne.0) goto 30
	  read(readln,12) cg,cn,avgcgs,ek,tc
12      format (2(3X,I7),3X,F7.3,3X,F17.15,3X,F10.5)	

	if (cn.ne.0) then  		
	    if (cg.eq.0) Then
	      atsum=atsum+cn
	    endif
	    atsum=atsum+cn*cg
	    eks=eks+ek*cn*cg
	    ocn=cn
	    ocg=cg
	    oavgcgs = avgcgs
	    oek=ek
	    otc=tc

	    if (flag.ne.1) Then 
	      if (cg.ge.nmax) flag = 1
	      if (ts.ge.warmup) count(cg) = count(cg) + cn
	    endif	
       	    
	    goto 22
	else
	    eks=eks-oek*ocn*ocg
	    cd=cd+ocg
	    if (ocg.eq.atsum) Then
	      write (2,10) ts,ocg,otc,atsum-ocg,0,real(cd)/i,oavgcgs
	    else
	      write (2,10) ts,ocg,otc,atsum-ocg,
     $        eks/((atsum-ocg)*1.5d0*cikb),real(cd)/i,oavgcgs
	    endif
	endif

	goto 32


*** SWITCH = 3
*** Old sizer (without the carrier gas coloumn) (SWITCH = 3)

	ELSEIF (SWITCH.eq.3) THEN
	oavgcgs = 0
	flag = -1
	fpn = 'pn.stillinger.old'
	open (3,file=fpn)
	
	read(1,'(A)',IOSTAT=ios) readln
	if (ios.ne.0) goto 30	
33	continue

	read (readln,*,IOSTAT=ios) ts

	atsum=0
	eks=0.0d0
23	read(1,'(A)',IOSTAT=ios) readln
	if (ios.ne.0) goto 30
	read(readln,13) cg,cn,ek,tc
13    format (2(3X,I7),3X,F17.15,3X,F10.5)	
	
	
	if (cn.ne.0) then
	   if (flag.ne.1) Then 
	     if (cg.ge.nmax) flag = 1
	     if (ts.ge.500000) count(cg) = count(cg) + cn
	   endif
	   atsum=atsum+cn*cg
	   eks=eks+ek*cn*cg
	   ocn=cn
	   ocg=cg
	   oek=ek
	   otc=tc
	   goto 23 
	else
	   eks=eks-oek*ocn*ocg
	   cd=cd+ocg
	   if (ocg.eq.atsum) then
	     write (2,10) ts,ocg,otc,atsum-ocg,0,real(cd)/i,oavgcgs
	   else
	     write (2,10) ts,ocg,otc,atsum-ocg,
     $        eks/((atsum-ocg)*1.5d0*cikb),real(cd)/i,oavgcgs
	    endif
	endif
	
	goto 33

	ELSE
	   write (*,*) 'Incorrect switch!'
	ENDIF

*** FINAL 
*** write-out

30	continue

	if (ios.ne.0) then	
	write (2,10) ts,ocg,otc,atsum-ocg,
     $        eks/((atsum-ocg)*1.5d0*cikb),real(cd)/i,oavgcgs
	endif

	do i = 0,nmax
	  if (count(i).ne.0) write (3,'(I5,2X,I12)') i,count(i)  
	enddo
	
	write (*,*) '-->> DONE << sizer wrote to:   ',fout
	write (*,*) '                        and:   ',fpn

	close (1)
	close (2)
	close (3)

	end





*** SUBROUTINES *************************************************************
*** Probabilities -----------------------------------------------------------
	subroutine probability(flag,ocg,nmax,ts,warmup,count)
	
	implicit none

	integer flag, ocg, nmax, ts, warmup,count(*)

	if (flag.ne.1) Then 
	   if (ocg.ge.nmax) flag = 1
	   if (ts.ge.warmup) count(ocg) = count(ocg) + 1
	endif
    
	return
	end
*** Yasuoka ----------------------------------------------------------------
	subroutine yasuoka()

	implicit none
	
	

	return
	end
	
	
*** HELP ROUTINE -----------------------------------------------------------

	subroutine helptext()
	
	implicit none
	
	write (*,*) 'Sizer Routine HELP'
	write (*,*) 'call: "sizer file.3d outputfile SWITCH NMAX WARMUP'
	write (*,*) 'SWITCH:'
	write (*,*) '   1: Stillinger'
	write (*,*) '   2: ten Wolde/Frenkel'
	write (*,*) '   3: Old CSI files (w/o carrier gas)'
	write (*,*) 'NMAX:'
	write (*,*) '   maximum size probability-counting'
	write (*,*) 'WARMUP:'
	write (*,*) '   warmup for probability-counting'
	stop
	end

*****************************************************************************
