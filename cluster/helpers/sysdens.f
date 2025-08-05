*------------------------------------------------------------------------------
* Programm zur Berechnung der Dichteverteilung im System
* Das System wird in kleine Teile unterteilt; es wird eine Statistik
* der Dichten erstellt
*------------------------------------------------------------------------------

	program sysdens

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'boxpp.inc'
	include 'filep.inc'
	include 'cvtpp.inc'

	character*40 fout,dummy
	integer i,j,step,maxi,ind,disc,divs
	double precision hx,hy,hz
	double precision rho,rhoAV,rhoSC,bstSC
	integer nb,ix,iy,iz
	parameter (nb=10)
	integer bx(0:nb*nb*nb),bst(0:mxa)

	dummy=''
	fout=''

	i=2
	call getarg(i,fout)
	i=3
	call getarg(i,dummy)
	read (dummy,*) disc
	i=4
	call getarg(i,dummy)
	read (dummy,*) divs
	if (divs.gt.nb) then
	  write (*,*) 'Only ',nb,' Divisions allowed'
	  stop
	endif

	call readdata

	hx=boxx/2
	hy=boxy/2
	hz=boxz/2

	open (1,file=fxyz)
	open (2,file=fout)

	maxi=ntim/ftxyz
	do j=0,nb*nb*nb
	  bx(j)=0
	enddo
	do j=0,mxa
	  bst(j)=0
	enddo

	do i=1,maxi
	  call readXYZ(step,1)

	if (i.gt.disc) then

	  do j=1,natom
	    ix=min0(int(x(j)*real(divs)/boxx),divs-1)
	    iy=min0(int(y(j)*real(divs)/boxy),divs-1)
	    iz=min0(int(z(j)*real(divs)/boxz),divs-1)
	    ind=ix+divs*iy+divs*divs*iz
	    if (ix.lt.0.or.iy.lt.0.or.iz.lt.0) then
	      write (*,*) 'eeek! Too Low!', ix,iy,iz
	      ind=divs*divs*divs-1
	    endif
	    bx(ind)=bx(ind)+1
	  enddo

	  do j=0,divs*divs*divs-1
	    bst(bx(j))=bst(bx(j))+1
	    bx(j)=0
	  enddo
	endif
	enddo

	bstSC=1.0d0/(real(maxi-disc)*real(divs)**4.5d0)
	rhoAV=real(natom)*40.0d0*cimu/real(boxx*boxy*boxz)
	rhoSC=dsqrt(1.0d0/real(divs*divs*divs))

	ind=0
	do i=0,natom
	  ind=ind+bst(i)
	  rho=real(i)*40.0d0*cimu*real(divs*divs*divs)/
     $	      real(boxx*boxy*boxz)
*	  rho=rhoAV+(rho-rhoAV)*rhoSC
	  write (2,20) i,rho,real(bst(i))*bstSC,
     $                 real(ind)*bstSC,bst(i),ind
	enddo

20	format (I7,3(2X,F15.5),2(2X,I7))
	close(1)
	close(2)

	end
