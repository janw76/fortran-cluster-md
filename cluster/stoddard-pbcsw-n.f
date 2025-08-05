*===============================================================================
*	DETERMINATION OF CLUSTER SIZES
*-------------------------------------------------------------------------------
* file: stoddard-pbcsw-n.f
*-------------------------------------------------------------------------------
* Stoddard-Routine, see Allen/Tildesley, S.288
* This routine uses the neighborlist to save calculation time!
* additionaly, the periodic boundaries can be switched off.
*-------------------------------------------------------------------------------
	subroutine stoddard(step)
	
	implicit none

	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'cvtpp.inc'
	include 'boxpp.inc'

	integer step
	integer i,j,k
	integer tcnum,ind
	integer cvtnt(mxnlist),cvtii(mxa+1)
	integer cvtn(mxnlist),cvti(mxa+1),cvtl(mxa+1)
	double precision dx,dy,dz
	double precision rij2,rijc2
	double precision bf1,bf2
	double precision ibx,iby,ibz
	double precision tx,ty,tz

*-------------------------------------------------------------------------------
* step          : Time step number
*...............................................................................
* i,j,k		: standard loop variable
* tcnum		: temporary cluster number
* ind		: Index for neighborlist
* dx,dy,dz	: components of distance between two atoms
* rij2		: squared distance of two atoms
* rijc2		: squared maximum distance of two atoms (cluster criterion)
* bf1,bf2	: box-factors for minimum image convention
* ibx,iby,ibz	: inverse box lengths in X,Y,Z-directions
* tx,ty,tz	: temporary variables for position of atom i
* dx,dy,dz	: components for distance of two atoms
* cvtn()	: *Full* neighborlist with atskin (with back references)
* cvtnt()	: neighborlist with atskin (*without* back references)
* cvti()	: indices for atoms in cvtn
* cvtii()	: counters for sub-neighborlists
* cvtl()	: length of sub-neighborlists
*-------------------------------------------------------------------------------

	cstep=step

	rijc2=atskin*atskin
	bf1=int(atrskin/atr)+2.0d0
	bf2=bf1+0.5d0
	ibx=1.0d0/boxx
	iby=1.0d0/boxy
	ibz=1.0d0/boxz	
	
	do i=1,natom
	  cl(i)=0
	  cc(i)=0
	  cs(i)=0
	  csi(i)=0
	  cek(i)=0
	  ceks(i)=0
	  cvtl(i)=0
	enddo

	ind=1
	do i=1,natom
	  cvtnt(ind)=-i
	  ind=ind+1
	  if (i.eq.natom) goto 110
	  tx=x(i)
	  ty=y(i)
	  tz=z(i)
	  k=atnidx(i)
100	  k=k+1
	  j=atnlist(k)
	  if (j.lt.0) goto 110
	    dx=tx-x(j)
	    dy=ty-y(j)
	    dz=tz-z(j)
	    if (bperx.eq.'X') dx=dx+boxx*(bf1-int(bf2+dx*ibx))
	    if (bpery.eq.'X') dy=dy+boxy*(bf1-int(bf2+dy*iby))
	    if (bperz.eq.'X') dz=dz+boxz*(bf1-int(bf2+dz*ibz))
	    rij2=dx*dx+dy*dy+dz*dz

	    if (rij2.le.rijc2) then
	      cvtnt(ind)=j
	      cvtl(i)=cvtl(i)+1
	      cvtl(j)=cvtl(j)+1
	      ind=ind+1
	    endif
	  goto 100
110	  continue
	enddo
	cvtnt(ind)=-(natom+1)

	ind=1
	do i=1,natom
	  cvti(i)=ind
	  cvtn(ind)=-i
	  ind=ind+cvtl(i)+1
	  cvtii(i)=cvti(i)
	enddo
	cvtn(ind)=-(natom+1)

	ind=1
510	i=-cvtnt(ind)
	ind=ind+1
	if (i.eq.natom+1) goto 500
520	j=cvtnt(ind)
	if (j.lt.0) goto 510
	ind=ind+1
	  cvtii(i)=cvtii(i)+1
	  cvtii(j)=cvtii(j)+1
	  cvtn(cvtii(j))=i
	  cvtn(cvtii(i))=j
	goto 520
500	continue

	cnum=0
	do k=1,natom
	  cnum=cnum+1
	  j=k
	  if (cc(j).eq.0) then
200	    cl(j)=cnum
	    cc(j)=1
	    cs(cnum)=cs(cnum)+1
	    cek(cnum)=cek(cnum)+
     $	      0.5d0*atm*(vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j))

	    ind=cvti(j)+1
	    do i=ind,ind+cvtl(j)-1
	      cl(cvtn(i))=cnum
	    enddo

	    do j=k+1,natom
	      if(cc(j).eq.0.and.cl(j).eq.cnum) goto 200
	    enddo
	  endif
	enddo	

	do i=1,cnum
	  csi(cs(i))=csi(cs(i))+1
	  ceks(cs(i))=ceks(cs(i))+cek(i)
	enddo
	
	return
	end
*-------------------------------------------------------------------------------
