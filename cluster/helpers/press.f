*------------------------------------------------------------------------------
* Programm zur Berechnung des Druckes im System
*------------------------------------------------------------------------------

	program press

	implicit none
	
	include 'timep.inc'
	include 'filep.inc'

	character*40 fout,avc
	integer i,maxi,avs
	double precision w1,w2,w3,w4,w5,w6,w7,w8
	double precision w9,w10,w11,psum

	fout=''

	i=2
	call getarg(i,fout)
	i=3
	call getarg(i,avc)
	read (avc,*) avs

	call readdata

	open (1,file=fave)
	open (2,file=fout)

	maxi=ntim/ftxyz
	psum=0.0d0

	do i=1,maxi
	  read(1,10) w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11
	  psum=psum+w11
	  if (mod(i+2*avs-1,avs).eq.0) then
	    write (2,20) w1,psum/real(avs)
	    psum=0.0d0
	  endif
	enddo

10 	format (I15,7(3X,F17.15),3X,F10.5,2(3X,F13.4))
20 	format (I15,3X,F11.4)

	close(1)
	close(2)

	end
