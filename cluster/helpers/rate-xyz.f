*------------------------------------------------------------------------------
* Programm zum erstellen eines Rates-Files aus dem xyz-file einer
* Simulation
*------------------------------------------------------------------------------

	program rater

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'timep.inc'
	include 'filep.inc'

	integer i,maxi,step

	call readdata

	open (1,file=fxyz)
	open (67,file=frat)

	maxi=ntim/ftxyz

	do i=1,maxi
	  write (*,*) 'Reading config:',i
	  call readXYZ(step,1)
	  if (i.eq.1) call stoddard(step-1)
	  call rates(step)
	enddo

	close(1)
	close(67)

	end
