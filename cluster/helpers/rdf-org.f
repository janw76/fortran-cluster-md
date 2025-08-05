*------------------------------------------------------------------------------
* Programm zur Berechnung der Paarverteilungsfunktion im System
*------------------------------------------------------------------------------

	program rdf

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'boxpp.inc'
	include 'filep.inc'

	integer mxst
	parameter (mxst=10000)
	character*40 fout,fstr,dummy
	integer i,j,k,maxi,step
	integer ind,disc
	double precision tx,ty,tz,dx,dy,dz
	double precision bf1,bf2,ibx,iby,ibz
	double precision rij,rstep,irstep
	integer gr(0:mxst)
	double precision pofq(0:mxst),dbq(0:mxst)
	double precision pi

	pi=4.0d0*datan(1.0d0)

	dummy=''

	i=2
	call getarg(i,fout)
	i=3
	call getarg(i,fstr)
	i=4
	call getarg(i,dummy)
	read (dummy,*) disc
	i=5
	call getarg(i,dummy)
	read (dummy,*) rstep

	irstep=1.0d0/rstep

	call readdata

	bf1=2.0d0+int(atrskin/atr)
	bf2=bf1+0.5d0

	ibx=1.0d0/boxx
	iby=1.0d0/boxy
	ibz=1.0d0/boxz

	open (1,file=fxyz)
	open (2,file=fout)

	maxi=ntim/ftxyz

	do i=0,mxst
	  gr(i)=0
	enddo

	do i=1,maxi
	  write (*,*) 'Reading Config:',i
	  call readXYZ(step,1)
	  if (i.gt.disc) then
	    do j=1,natom
	      tx=x(j)
	      ty=y(j)
	      tz=z(j)
	      do k=j+1,natom
	        dx=tx-x(k)
	        dy=ty-y(k)
	        dz=tz-z(k)
	        dx=dx+boxx*(bf1-int(bf2+dx*ibx))
                dy=dy+boxy*(bf1-int(bf2+dy*iby))
                dz=dz+boxz*(bf1-int(bf2+dz*ibz))
	        rij=dsqrt(dx*dx+dy*dy+dz*dz)
	        ind=int(rij*irstep)
	        if (ind.le.mxst) then
	          gr(ind)=gr(ind)+2
	        else
*		  write (*,*) 'Eeek...',rij
	        endif
	      enddo
	    enddo
	  endif
	enddo

	gr(0)=gr(0)+natom

	dx=4.0d0*pi*real(natom)/(boxx*boxy*boxz*3.0d0)
*	dx=real(natom-600)/(35.0d0**3)
	do i=0,mxst
	  dz=real(i)*rstep
	  dy=dx*(dz**3-(dz-rstep)**3)
	  write (2,10) dz,gr(i)/real(maxi-disc),dy,
     $                real(gr(i))/(real(maxi-disc)*real(natom)*dy)
	enddo

10	format(F10.4,2X,I8,2X,F10.4,2X,F10.4)

	do i=0,mxst
	  dbq(i)=10.0d0**(-10.0d0+dble(i)*0.02d0)
	enddo

	do i=1,mxst
	  if (gr(i).gt.0) then
	    do j=1,mxst
	      pofq(j)=pofq(j)+dble(gr(i))*dsin(dbq(j)*dble(i))/dble(i)
	    enddo
	  endif
	enddo

	do i=1,mxst
	  pofq(i)=(pofq(i)+dble(natom))/(dbq(i)*dble(natom)**2)
	  if (pofq(i).lt.0) pofq(i)=1.0d-10
	enddo

	open (3,file=fstr)

	do i=1,mxst
	  write (3,1) i,dbq(i),log10(dbq(i)),pofq(i),log10(pofq(i))
	enddo
1       format (I4,'  ',2(F15.10,'  ',F15.8,'  '))

	close(1)
	close(2)
	close(3)
	end
