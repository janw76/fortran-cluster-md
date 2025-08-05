*===============================================================================
*	CHECKING IF A CLUSTER NUCLEATED IN THE SYSTEM
*-------------------------------------------------------------------------------
* file: chk-nucleation.f
* date: 14-jul-2005 Jan Wedekind 
*-------------------------------------------------------------------------------
* Checking, if there exists a cluster>"csover" for at least "clwait" timesteps.
* If so, this information can be used to automatically stop the simulation
* because the remaining simulation steps are not interesting in nucleation
* experiments. Values "csover", "clwait" controlled through the 3d Input file
* This routine uses the same timer as wrt-csi.f! If this one is set to 0 
* this module will NOT work!
*-------------------------------------------------------------------------------

	subroutine chk_nucleation(step)
	
	implicit none

	include 'const.inc'
	include 'atomp.inc'
	include 'cvtpp.inc'
	
	integer i,step,tim	

*-------------------------------------------------------------------------------
* i		: loop variable
* step		: Time step number
*...............................................................................
* tim		: temporary timestep when cluster>csover is found
*-------------------------------------------------------------------------------

	
 	tim=0
  	
	do i=csover,natom
	  if (csi(i).gt.0) then
      	    tim=step
	    goto 10
    	  endif
  	enddo

10  	if (tim.gt.0) then
	  if (clbrkt.eq.0) then
            clbrkt=tim
	    else
            if ((tim-clbrkt).ge.clwait) then
	      write(*,*) 
	      write(*,*) 'Found a cluster larger than: ',csover,'!'
	      write(*,*) '---> ABORTING simulation at timestep: ',step,'!'
	      write(*,*) 
	      step=-step
            endif
      	  endif
	else 
	  clbrkt=0
	endif
	
	end

*-------------------------------------------------------------------------------
