*===============================================================================
*	WRITE AVERAGES
*===============================================================================
* file: wrt-ave.f
*-------------------------------------------------------------------------------
* This routine writes several interesting properties - file format:
* 1: Time step
* 2: total kinetic energy [10^-17 J]
* 3: total potential energy [10^-17J]
* 4: total system energy [10^-17J]
* 5: wall interaction energy [10^-17 J]
* 6: wall interaction energy, sum over last n timesteps
* 7: force on walls
* 8: force on walls, sum over last n timesteps
* 9: System Temperature
*10: pressure, based on 7 [pa]
*11: pressure, based on 8 [pa]
*-------------------------------------------------------------------------------

	subroutine wrtAVE(step,fnum)

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'energ.inc'
	include 'filep.inc'
	include 'boxpp.inc'

	integer step,fnum
	
*-------------------------------------------------------------------------------
* step		: Schrittnummer
* fnum		: Filenumber fuer Ausgabe
*-------------------------------------------------------------------------------

	write (fnum,10) step,ekin,epot,ekin+epot,
     $  ewall,ewallsum,
     $  fwall,fwallsum,
     $	ekin/(natom*1.5d0*cikb),
     $  1.0d13*fwall/(2.0d0*(boxx*boxy+boxy*boxz+boxx*boxz)),
     $  1.0d13*fwallsum/(2.0d0*ftave*(boxx*boxy+boxy*boxz+boxx*boxz))
10 	format (I15,7(3X,F17.14),3X,F10.5,2(3X,F13.4)) 
	
	ewallsum=0.0d0
	fwallsum=0.0d0

	return
	end
*-------------------------------------------------------------------------------
