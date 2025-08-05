*===============================================================================
*	DETERMINATION OF CLUSTER SIZES
*-------------------------------------------------------------------------------
* file: stoddard.f
*-------------------------------------------------------------------------------
* Stoddard-Routine, see Allen/Tildesley, S.288
*-------------------------------------------------------------------------------
	subroutine stoddard(step)

	implicit none

	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'cvtpp.inc'
	
	integer step
	integer i,j,k
	integer flag
	double precision dx,dy,dz
	double precision rij2,rijc2
	
*-------------------------------------------------------------------------------
* step		: Time step number
*...............................................................................
* i,j,k		: standard loop variables
* flag		: control flag for calculations
* dx,dy,dz	: components of distance between two atoms
* rij2		: squared distance of two atoms
* rijc2		: squared maximum distance of two atoms (cluster criterion)
*-------------------------------------------------------------------------------

	cstep=step

	rijc2=atskin*atskin

	do i=1,natom
	  cl(i)=0
	  cs(i)=0
	  cc(i)=0
	  csi(i)=0
	  cek(i)=0
	  ceks(i)=0
	enddo

	cnum=1
	do i=1,natom
	  if (cl(i).eq.0) then
	    cl(i)=cnum
10	    flag=0
	    do j=i,natom
	      if (cc(j).eq.0.and.cl(j).eq.cnum) then
	        cc(j)=1
	        cs(cnum)=cs(cnum)+1
		cek(cnum)=cek(cnum)+
     $	  0.5d0*atm*(vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j))
	        do k=i,natom
	          if (cl(k).eq.0) then
	            dx=x(k)-x(j)
	            dy=y(k)-y(j)
	            dz=z(k)-z(j)
	            rij2=dx*dx+dy*dy+dz*dz
	            if (rij2.le.rijc2) then
	              cl(k)=cnum
	              flag=1
	            endif
	          endif
	        enddo
	      endif
	    enddo
	    if (flag.eq.1) goto 10
	    cnum=cnum+1
	  endif
	enddo
	cnum=cnum-1

	do i=1,cnum
	  csi(cs(i))=csi(cs(i))+1
	  ceks(cs(i))=ceks(cs(i))+cek(i)
	enddo
	
	return
	end
*-------------------------------------------------------------------------------
