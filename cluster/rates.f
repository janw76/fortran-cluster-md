*===============================================================================
*	CONDENSATION AND EVAPORATION RATES 
*-------------------------------------------------------------------------------
* file: rates.f
* date: 31-jan-2001
*-------------------------------------------------------------------------------
* In this routine the condensation and evaporation rates of the biggest
* cluster are calculated and shown in a statistic outputfile.
* CAUTION: Only use in systems with ONE big cluster !
*-------------------------------------------------------------------------------

	subroutine rates(step)
	
	implicit none

	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'cvtpp.inc'

	integer step
	integer ocs(mxa),ocl(mxa)
	double precision ocek(mxa)
	integer mcs,mcln,omcs,omcln,en,cn
	integer i
	double precision mct,omct

*-------------------------------------------------------------------------------
* i		: standard loop variable
* step		: time step number
* ocs()		: old size of cluster N
* ocl()		: old clusterlabel of atom
* ocek()	: old kinetic energy of cluster N
* msc		: size of max. cluster
* mcln		: clusterlabel of max. cluster
* omsc          : size of old max. cluster
* omcln         : clusterlabel of old max. cluster
* en		: evaporation rate
* cn 		: condensation rate
* mct		: temperature of maximal sized cluster
* omct		: temperature of old maximal sized cluster
*-------------------------------------------------------------------------------

	do i=1,natom
	  ocs(i)=cs(i)
	  ocl(i)=cl(i)
	  ocek(i)=cek(i)
	enddo

	if (cstep.ne.step) then
	  call stoddard(step)
	else
	  write(*,*) '*** WARNING *** No old Stoddard-list !!'
	endif

	mcs=0
	omcs=0
	do i=1,natom
	  if (ocs(i).gt.omcs) then
	    omcs=ocs(i)
            omcln=i
	  endif
	  if (cs(i).gt.mcs) then
	    mcs=cs(i)
	    mcln=i
	  endif
	enddo
	
	cn=0
	en=0
	do i=1,natom
	  if (cl(i).eq.mcln .and. ocl(i).ne.omcln) cn=cn+1
	  if (cl(i).ne.mcln .and. ocl(i).eq.omcln) en=en+1
	enddo

	if (mcs.ne.omcs) write (*,1) step,omcs,mcs,omcln,mcln,cn,en
1	format (I12,4(2X,I5),2(2X,I2))

	mct=cek(mcln)/(dble(mcs)*1.5d0*cikb)
	mcts=mcts+mct
	omct=ocek(omcln)/(dble(omcs)*1.5d0*cikb)
	omcts=omcts+omct

	cen=cen+en
	ccn=ccn+cn
	if (ctrat.ne.1) then
	  if (mod(step,ctrat).eq.0) then
	    write(67,10) step,omcs,mcs,
     $  mcts/real(ctrat),omcts/real(ctrat),ccn,cen,mcln,omcln
	    ccn=0
	    cen=0
	    mcts=0.0d0
	  endif
	else
	  if ((en+cn).ne.0) then 
	    write(67,10) step,omcs,mcs,mct,omct,cn,en,mcln,omcln
	  endif
	endif

10	format (I12,2(1X,I6),2(1X,F10.4),2(1X,I3),2(1X,I6))

	return
	end
*-------------------------------------------------------------------------------
