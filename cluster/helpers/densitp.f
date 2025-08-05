*------------------------------------------------------------------------------
* Programm zur Berechnung des Dichteprofils des groessten Clusters
*------------------------------------------------------------------------------

	program density

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'boxpp.inc'
	include 'filep.inc'
	include 'cvtpp.inc'

	character*40 fout,dummy
	integer i,j,k,step,maxi,mcs,disc
	double precision hx,hy,hz,cx,cy,cz,dx,dy,dz
	double precision rstep,rrat,msum(0:20000),ks
	double precision ibx,iby,ibz

	i=2
	call getarg(i,fout)
	i=3
	call getarg(i,dummy)
	read (dummy,*) disc
	i=4
	call getarg(i,dummy)
	read (dummy,*) rstep

	call readdata

	open (1,file=fxyz)
	open (2,file=fout)
*	open (3,file='cmass.xyz')

	maxi=ntim/ftxyz

	do i=0,20000
	  msum(i)=0
	enddo

	ibx=1.0d0/boxx
	iby=1.0d0/boxy
	ibz=1.0d0/boxz

	do i=1,maxi
	  call readXYZ(step,1)

	if (i.gt.disc) then

	  call stoddard()

	  mcs=1
	  k=1
	  do j=1,cnum
	    if (cs(j).gt.k) then
	      mcs=j
	      k=cs(j)
	    endif
	  enddo

	  cx=0.0d0
	  cy=0.0d0
	  cz=0.0d0

	  k=1
	  do j=1,natom
	    if (cl(j).eq.mcs) then
	      if (k.eq.1) then
		k=0
		hx=x(j)
		hy=y(j)
		hz=z(j)
	      endif
	      dx=x(j)-hx
	      dy=y(j)-hy
	      dz=z(j)-hz
	      cx=cx+dx+boxx-boxx*int(1.5d0+dx*ibx)
	      cy=cy+dy+boxy-boxy*int(1.5d0+dy*iby)
	      cz=cz+dz+boxz-boxz*int(1.5d0+dz*ibz)
	    endif
	  enddo
	  cx=dmod(dmod(hx+cx/real(cs(mcs)-1),boxx)+boxx,boxx)
	  cy=dmod(dmod(hy+cy/real(cs(mcs)-1),boxy)+boxy,boxy)
	  cz=dmod(dmod(hz+cz/real(cs(mcs)-1),boxz)+boxz,boxz)

*	write (3,*) 2
*	write (3,*) step
*	write (3,1) 'N',cx-boxx/2,cy-boxy/2,cz-boxz/2
*1	format (A2,3(2X,F10.4))

	  do j=1,natom
	    dx=cx-x(j)
	    dy=cy-y(j)
	    dz=cz-z(j)
	    dx=dx+boxx-boxx*int(1.5d0+dx*ibx)
	    dy=dy+boxy-boxy*int(1.5d0+dy*iby)
	    dz=dz+boxz-boxz*int(1.5d0+dz*ibz)
	    rrat=dsqrt(dx*dx+dy*dy+dz*dz)
	    k=int(rrat/rstep)
	    if (k.lt.20000) msum(k)=msum(k)+atm
	  enddo
	endif
	write (*,*) i,' of',maxi

	enddo

	cx=0.0d0
	do j=0,20000
	  cx=cx+msum(j)
	  cy=real(maxi-disc)*4.0d-3*3.1415027*
     $	  ((rstep*real(j))**3)/3.0d0
          ks=real(maxi-disc)*4.0d-3*3.1415927*
     $    ((rstep*real(j))**3-(rstep*real(j-1))**3)/3.0d0
	  write (2,20) j,rstep*j,int(msum(j)),msum(j)/ks,int(cx),cx/cy
	enddo

*------
* Output:
* 1: i
* 2: r [10^-10 m]
* 3: m(r) [10^-27 kg]
* 4: m(r)/r^2dr [kg m^-3]
* 5: sum(m(r)) [10^-27 kg]
* 6: sum(m(r))/V [kg m^-3]


20 	format (I4,2X,F11.5,2X,I10,2X,F11.5,2X,I10,F11.5)

	close(1)
	close(2)

	end
