*===============================================================================
*
*		C L U S T E R
*
*===============================================================================
* program for simulation of small lennard-Jones atomic clusters
* Author	: Stephan Wonczak
* Version 	: serial version for single processors
* Created	: 9.7.1997
*-------------------------------------------------------------------------------
* file		: main.f
* last changed	: 07-may-2006 Jan Wedekind
* change	: added call chk_nucleation(i)
* change	: added Nose-Hoover
*
* NOSE HOOVER ONLY WORKS WITH ONE SORT OF ATOMS!!!
* 11/10/07 (JW):
*   - added first passage time routine
*
*-------------------------------------------------------------------------------

	program cluster
	
	implicit none

	include 'const.inc'
	include 'atomc.inc'
	include 'filep.inc'
	include 'timep.inc'
	include 'energ.inc'
	include 'boxpp.inc'
	include 'atomp.inc'
	include 'cvtpp.inc'

	integer i

*-------------------------------------------------------------------------------
* i,j		: loop variable
*-------------------------------------------------------------------------------

	call kpsInit(0)
	call readdata()
	
 	if (ftxyz.ne.0) open (61,file=fxxyz)
        if (ftvel.ne.0) open (62,file=fxvel) 
        if (ftacc.ne.0) open (63,file=fxacc) 
        if (ftcsi.ne.0) open (65,file=fxcsi) 
        if (ftave.ne.0) open (66,file=fxave) 
        if (ftrat.ne.0) open (67,file=fxrat)
	if (ftcsi.ne.0) open (68,file=fxcsif)
	if (ftcsi.ne.0) open (69,file=fxfpt)

*** Screen output of simulation details
	write (*,'(A,I6)'  ) 'Total number of atoms N: ',natom
        write (*,'(A,F9.2)') 'Box volume [nm^3]      : ',
     $    boxx*boxy*boxz/1000d0
        write (*,'(A,F9.2)') 'System temperature [K] : ',att
        write (*,'(A,F9.2)') 'Simulation time [ns]   : ',dt*ntim/1d6
        write (*,'(A,A)')    'use thermostat?        : ',ethst
	write (*,*)	
	write (*,*) 'Sort Sym Number Mass [u] Thermostat'
        do i=1,atsorts
	  write (*,'(X,I4,2X,A2,1X,I6,2X,F5.1,8X,A1)') i,apsym(i),apnum(i),
     $      apmas(i)/cimu,apth(i)
        enddo
	write (*,*)
 	if (clnuctrl.eq.'X') then
	  write (*,'(A,I3)') 
     $      'Simulation will be aborted if a cluster exceeds n=',csover
	  write (*,*)
        endif

*** Initialization
	if (frctrl.eq.'X') then
	  call readRST(i,64)
	  if (ethst.eq.'N') then
	  write (*,'(A,F14.5)') 'Q (internal): ',Q
	  write (*,'(A,F14.5)') 'Q (reduced) : ',Q/(apmas(1)*apsig(1)*apsig(1))
	    write (*,*)
	    call stopdead(-1)
	  endif
	  call stopdead(-1)
	else
	  call initconf
	  i=btim
	endif
	
*** For VMD: after init or read we set the last two atoms of species 1 to the
*** atom symbols "S" and "F" respectively. This order is demanded by VMD and
*** must be present in the first time step in the XYZ file!
*	write (*,*) apnum(1)-1, apnum(1)
*	apsymi(apnum(1)-1) = 'S'
*	apsymi(apnum(1)) = 'F'

*** Initial values of Nose-Hoover coefficients
	s = 0
	PSI = 0
	sumv2 = 2d0*ekin/apmas(1)		
	if (ethst.eq.'N') then
	 PSI = (apmas(1)*sumv2 - 3.0d0*natom*cikb*att) / Q
	else
	 PSI = 0
	 s = 0
	endif

	H = ekin + epot + PSI*PSI*Q/2 + 3*natom*cikb*att*s
	
*	if (ftxyz.ne.0) then 
*	    if (mod(i,ftxyz).eq.0) call wrtXYZ(i,61)
*	endif

	call stoddard(-1)

	write (*,*)

*** start of simulation	loop
10	continue

	  if (ftrat.ne.0) call rates(i)

	  if (ftxyz.ne.0) then 
	    if (mod(i,ftxyz).eq.0) call wrtXYZ(i,61)
	  endif
	  if (ftvel.ne.0) then 
	    if (mod(i,ftvel).eq.0) call wrtVEL(i,62)
	  endif
	  if (ftacc.ne.0) then
	    if (mod(i,ftacc).eq.0) call wrtACC(i,63)
	  endif
	  if (ftrst.ne.0) then 
	    if (mod(i,ftrst).eq.0) call wrtRST(i,64) 
	  endif
	  if (ftcsi.ne.0) then 
	    if (mod(i,ftcsi).eq.0) call wrtCSI(i,65)
	  endif
	  if (ftave.ne.0) then 
	    if (mod(i,ftave).eq.0) call wrtAVE(i,66)
	  endif

*** screen ouput during simulation
	  if (mod(i,sctim).eq.0) then
	    H = ekin + epot + PSI*PSI*Q/2 + 3*natom*cikb*att*s
	    write (*,'(I5,2X,4(F18.10,2X),A40)')
     $	  i/sctim,2d0*ekin/(3*natom*cikb),ekin,epot,H,fcsi
	  endif

	  if (ethst.eq.'X') call thermostat(i)

	  if (bshr.eq.'X') then
	    if (mod(i,boxt).eq.0) call resize(i)
	  endif
	  if (chgNctl.eq.'X') then
	    if (mod(i,chgNtim).eq.0) call chgNatom(i)
	  endif

*** check for a cluster in the system
	  if (clnuctrl.eq.'X') then 
	    if (mod(i,ftcsi).eq.0) then
	      call chk_nucleation(i)
	      if (i.lt.0) goto 20
	    endif
	  endif
	
	  call integrator(i)

	  i=i+1

	if (i.le.ntim+btim) goto 10

20	continue

*** write out first passage times
	
*	call wrtFPT()

	end

*-------------------------------------------------------------------------------
