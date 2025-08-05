*------------------------------------------------------------------------------
* Programm zum Testen der Stoddard-Routine!
*------------------------------------------------------------------------------

	program stotest

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'filep.inc'
	include 'cvtpp.inc'

	integer i,s

	call readdata
	open (1,file=fxyz)
	open (2,file='stotest.csi')

	do i=1,20
	  s=i
	  call readXYZ(s,1)
	  call neighbor(-1)
	  call wrtCSI(s,2)
	enddo 

	close (1)
	close (2)

	end
