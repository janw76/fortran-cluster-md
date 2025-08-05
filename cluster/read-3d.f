*===============================================================================
*	READ DATA FROM INPUT FILE
*===============================================================================
* file:	read-3d.f
* date: 3-jan-2001
*-------------------------------------------------------------------------------
* This routine reads the relevant data form the input (3d)-file
* change 15-jul-2005: updated to 3.d File version 1.12
*-------------------------------------------------------------------------------

	subroutine readdata()
	
	implicit none
	
	include 'const.inc'
	include 'filep.inc'
	include 'boxpp.inc'
	include 'timep.inc'
	include 'atomp.inc'
	include 'energ.inc'
	include 'dmmpp.inc'
	include 'cvtpp.inc'	
	include 'atomc.inc'

	character*80 dummy,infile
	character*4 frenkel,fpt
	integer i,j,plen,fspc
	real ver,fver

*-------------------------------------------------------------------------------
* i,j		: standard loop variables
* dummy		: dummy variable for skipping division lines
* infile 	: name of the input file, read from command line
* fnam		: temporary variable for filnames
* plen		: length of pathname
* ver		: minimum version number readable by this routine
* fver		: version of the input file (read from file)
*-------------------------------------------------------------------------------

	ver=1.02
	i=1
	dummy=''
	infile=''
	frenkel='fsi'
	fpt='fpt'

	call getarg(i,infile)
	if (infile.eq.dummy) then
		write(*,*) 'Please supply an input file name!'
		stop
	else
		write(*,*) 'Reading data from file: ',infile
	endif
	
	open (1,file=infile)
	
	read (1,'(A)')	dummy
	read (1,'(41X,F8.6)') fver
	if (fver.lt.ver) then
		write(*,*) 'Inputfile is unreadable due to age!'
		write(*,*) 'File version number: ',fver
		write(*,*) 'Minimum Version    : ',ver
		stop
	endif
	read (1,'(A)') dummy
	read (1,'(41X,A)') fpath
	plen=fspc(fpath)
	read (1,'(41X,A)') fxyz
	fxxyz=fpath
	fxxyz(plen:)=fxyz
        read (1,'(41X,I10)') ftxyz
	
	read (1,'(41X,A)') fvel
	fxvel=fpath
	fxvel(plen:)=fvel
	read (1,'(41X,I10)') ftvel
        read (1,'(41X,A)') facc
	fxacc=fpath
	fxacc(plen:)=facc
        read (1,'(41X,I10)') ftacc
	read (1,'(41X,A)') frst
	fxrst=fpath
	fxrst(plen:)=frst
        read (1,'(41X,I10)') ftrst
	read (1,'(41X,A)') fcsi
	fxcsi=fpath
	fxcsi(plen:)=fcsi
	ffpt = fpt
	fxfpt = fpath
	fxfpt(plen:)=fpt
	fcsif = frenkel
	fxcsif = fpath
	fxcsif(plen:) = fcsif
	
        read (1,'(41X,I10)') ftcsi
	read (1,'(41X,A)') fave
	fxave=fpath
	fxave(plen:)=fave
        read (1,'(41X,I10)') ftave
	if (fver.ge.1.06) then
	  read (1,'(41X,A)') frat
	  fxrat=fpath
	  fxrat(plen:)=frat
          read (1,'(41X,I10)') ftrat
	else
	  ftrat=0
	endif
	cstep=-1
	ctrat=ftrat
	cen=0
	ccn=0
	read (1,'(A)') dummy
	read (1,'(41X,F14.5)') dt
	read (1,'(40X,I10)') ntim
	read (1,'(41X,I10)') btim
	read (1,'(41X,I10)') sctim 
	read (1,'(A)') dummy
	read (1,'(41X,F14.5)') boxx
	read (1,'(41X,F14.5)') boxy
	read (1,'(41X,F14.5)') boxz
	if (fver.ge.1.08) then
	  read (1,'(41X,A)') bperx
	  read (1,'(41X,A)') bpery
	  read (1,'(41X,A)') bperz
	else
	  bperx='X'	  
	  bpery='X'	  
	  bperz='X'	  
	endif
	read (1,'(41X,F14.5)') bpot
	read (1,'(41X,F14.5)') beng
	read (1,'(41X,A)') bcor
	if (fver.ge.1.05) then
	  read (1,'(41X,A)') bshr
	  read (1,'(41X,F14.5)') boxxn
	  read (1,'(41X,F14.5)') boxyn
	  read (1,'(41X,F14.5)') boxzn
	  read (1,'(41X,F14.5)') boxxs
	  read (1,'(41X,F14.5)') boxys
	  read (1,'(41X,F14.5)') boxzs
	  read (1,'(41X,I10)') boxt
	else
	  boxxn=boxx
	  boxyn=boxy
	  boxzn=boxz
	  boxxs=0.0
	  boxys=0.0
	  boxzs=0.0
	  bshr='-'
	  boxt=1000
	endif
	boxx=boxx*10
	boxy=boxy*10
	boxz=boxz*10
	boxxn=boxxn*10
	boxyn=boxyn*10
	boxzn=boxzn*10
	boxxs=abs(boxxs*10)
	boxys=abs(boxys*10)
	boxzs=abs(boxzs*10)
	beng=beng*cikb
	read (1,'(A)') dummy
	if (fver.lt.1.10) read (1,'(41X,I10)') atsorts
	read (1,'(41X,A)') frctrl
	read (1,'(41X,A)') finp
	fxinp=fpath
	fxinp(plen:)=finp
	read (1,'(41X,A)') ethst

	if (fver.ge.1.2) then
	  read (1,'(41X,F14.5)') Q
	endif

	read (1,'(41X,A)') rndctrl
	if (fver.ge.1.07) then
	  read (1,'(41X,A)') eclth
          read (1,'(41X,I10)') ecltim
	else
	  eclth='-'
	  ecltim=1000
	endif
	if (fver.ge.1.10) then
	  read (1,'(A)') dummy
	  read (1,'(41X,F14.5)') att
	  read (1,'(41X,F14.5)') atskin
	  if (fver.ge.1.30) then
	    read (1,'(41X,I10)') atnenum
	  endif 
	  read (1,'(41X,F14.5)') atrcut
	  read (1,'(41X,F14.5)') atrskin
	  read (1,'(41X,I10)') nbtim
	
	  if (fver.ge.1.11) then
	    read (1,'(41X,I10)') chgNtim
	  else
	    chgNtim=0
	  endif
	  if (fver.ge.1.12) then
	    read (1,'(A)') dummy
	    read (1,'(41X,A)') clnuctrl
	    read (1,'(41X,I10)') csover
	    read (1,'(41X,I10)') clwait
	  endif
	  read (1,'(A)') dummy
	  read (1,'(41X,I10)') atsorts
	endif

	read (1,'(A)') dummy
	if (fver.lt.1.10) then
	  read (1,'(41X,A)') ats
	  read (1,'(41X,I10)') natom
	  read (1,'(41X,F14.5)') atm
	  read (1,'(41X,F14.5)') atr
	  read (1,'(41X,F14.5)') ate
	  read (1,'(41X,F14.5)') atpr
	  read (1,'(41X,F14.5)') atpa
	  read (1,'(41X,F14.5)') att
	  read (1,'(41X,F14.5)') atskin
	  read (1,'(41X,F14.5)') atrcut
          read (1,'(41X,F14.5)') atrskin
          read (1,'(41X,I10)') nbtim
	  atm=atm*cimu
	  atr=atr*10.0d0
	  ate=ate*cikb
	endif

	if (fver.ge.1.10) then
	  natom=1
	  do i=1,atsorts
	    read (1,'(41X,A)') apsym(i)
	    read (1,'(41X,I10)') apnum(i)
	    read (1,'(41X,F14.5)') apmas(i)
	    read (1,'(41X,F14.5)') apsig(i)
	    read (1,'(41X,F14.5)') apeps(i)
	    read (1,'(41X,F14.5)') appr(i)
	    read (1,'(41X,F14.5)') appa(i)
	    read (1,'(41X,A)') apth(i)
	    if (fver.ge.1.11) then
	      read (1,'(41X,A)') apNch(i)
	      read (1,'(41X,I10)') apnnum(i)
	    else
	      apNch(i)='-'
	      apnnum(i)=apnum(i)
	    endif
	    read (1,'(A)') dummy
	    apmas(i)=apmas(i)*cimu
	    apsig(i)=apsig(i)*10.0d0
	    apeps(i)=apeps(i)*cikb
	    do j=natom,natom+apnum(i)-1
	      api(j)=i
		  apsymi(j) = apsym(i)
	    enddo
	    natom=natom+apnum(i)
	  enddo
	  Q = Q*apmas(1)*apsig(1)*apsig(1)
	  natom=natom-1	  	
	endif

	atskin=atskin*atr
	atrcut=atrcut*atr
	atrskin=atrskin*atr

	do i=0,255
	  dmmh(i)=(atrskin-atrcut)/real(nbtim)
	enddo
	dmmst=0
	dmmi=0
	dmmsu=256.0d0*(atrskin-atrcut)/real(nbtim)
	dmms=0.0d0
	ewallsum=0.0d0
	fwallsum=0.0d0
	
	clookfor = 1
	clookforf= 1
	cbiggest = 0
	cbiggestf= 0
	
	do i=1,natom
	  x0(i)=0.0d0
	  x1(i)=0.0d0
	  x2(i)=0.0d0
	  x3(i)=0.0d0
	  x4(i)=0.0d0
	  x5(i)=0.0d0
	  y0(i)=0.0d0
	  y1(i)=0.0d0
	  y2(i)=0.0d0
	  y3(i)=0.0d0
	  y4(i)=0.0d0
	  y5(i)=0.0d0
	  z0(i)=0.0d0
	  z1(i)=0.0d0
	  z2(i)=0.0d0
	  z3(i)=0.0d0
	  z4(i)=0.0d0
	  z5(i)=0.0d0
	  cfpt(i)=-1
	  cfptf(i)=-1
	enddo

	return
	end

*-------------------------------------------------------------------------------

        integer function fspc(a) 
        implicit none

        character*(*) a
        integer i,p1

        p1=0

        do i=1,len(a) 
          if (a(i:i).eq.' '.and.p1.eq.0) p1=i
        enddo
        fspc=p1
	return
        end
*-------------------------------------------------------------------------------
