*------------------------------------------------------------------------------
* Programm zum Selektieren einer Atomsorte aus den xyz-Datenfiles
* SYNTAX: selats *.3d 1[2] outfile
* geaendert: 18.11.2002  Marc Hamacher 
*------------------------------------------------------------------------------

	program selats

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'filep.inc'

	character*2  aslst(mxa)
	character*40 fout,asnew
	integer i,step,maxi,awant
	
	common /SPLIT/ aslst,awant,asnew

	i=2
	call getarg(i,fout)
	read (fout,*) awant

	i=3
	call getarg(i,fout)

	i=4
	call getarg(i,asnew)

	call readdata

	open (1,file=fxyz)
	open (2,file=fout)

	maxi=ntim/ftxyz + 1
	
	write(*,*) "You want atom's with symbol "// apsym(awant)
	
        do i=1,maxi
       	  call readXYZ(step,1)
	  call wrtXYZ(step,2)	  
        enddo

	close(1)
	close(2)

	end
