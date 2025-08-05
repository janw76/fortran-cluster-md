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
	integer i,j,k,step,maxi,nat,ind,mcs,disc
	double precision cx,cy,cz
	double precision rstep,rrs,rrsi,rrat,msum(0:2000),ks

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

	maxi=ntim/ftxyz

	do i=0,2000
	  msum(i)=0
	enddo

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
	  cx=0
	  cy=0
	  cz=0
	  do j=1,natom
	    if (cl(j).eq.mcs) then
	      cx=cx+x(j)
	      cy=cy+y(j)
	      cz=cz+z(j)
	    endif
	  enddo
	  cx=cx/real(cs(mcs))
	  cy=cy/real(cs(mcs))
	  cz=cz/real(cs(mcs))

	  do j=1,natom
	    rrat=dsqrt((x(k)-cx)*(x(k)-cx)+
     $		       (y(k)-cy)*(y(k)-cy)+
     $		       (z(k)-cz)*(z(k)-cz))
	    k=int(rrat/rstep)
	    if (k.lt.2000) msum(k)=msum(k)+atm
	  enddo

*	  write (2,'(I7,$)') step
*	  do j=1,500
*	    rrs=(rstep*real(j))**2
*	    rrsi=(rstep*real(j-1))**2
*	    msum(j)=0.0d0
*	    do k=1,natom
*	      if (cl(k).eq.mcs) then
*		rrat=(x(k)-cx)**2+(y(k)-cy)**2+(z(k)-cz)**2
*		if (rrat.le.rrs .and. rrat.gt.rrsi) then 
*		  msum(j)=msum(j)+atm
*		endif
*	      endif
*	    enddo
*	    ks=4.0d0*3.1415927*
*     $	    ((rstep*real(j))**3-(rstep*real(j-1))**3)/3.0d0
*	    write (2,'(1X,F10.3,1X,F10.3,$)')
*     $	    msum(j)/ks,msum(j)/atm
*	    write (2,'(1X,F10.3),$)') msum(j)/ks
*	  enddo
*	  write (2,*)
	endif

*	write (*,'(A,$)') '.'
*	if (mod(i,maxi/10).eq.0) write (*,'(I5,$)') i
	write (*,*) i,' of',maxi

	enddo

	cx=0.0d0
	do j=0,2000
	  cx=cx+msum(j)
	  cy=real(maxi-disc)*4.0d-3*3.1415027*
     $	  ((rstep*real(j))**3)/3.0d0
          ks=real(maxi-disc)*4.0d-3*3.1415927*
     $    ((rstep*real(j))**3-(rstep*real(j-1))**3)/3.0d0
	  write (2,20) rstep*j,int(msum(j)),msum(j)/ks,int(cx),cx/cy
	enddo

*------
* Output:
* 1: r
* 2: m(r)
* 3: m(r)/r^2dr
* 4: sum(m(r))
* 5: sum(m(r))/V


20 	format (F11.5,2X,I10,2X,F11.5,2X,I10,F11.5)

	close(1)
	close(2)

	end
