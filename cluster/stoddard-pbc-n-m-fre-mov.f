*===============================================================================
*	DETERMINATION OF CLUSTER SIZES
*-------------------------------------------------------------------------------
*  This routine just writes out the liquid monomers!
*
*
*
* file: stoddard-pbc-n.f
* date: 08-feb-2002
*-------------------------------------------------------------------------------
* Stoddard-Routine, see Allen/Tildesley, S.288
* This routine uses the neighborlist to save calculation time!
* additionally this routine handles multiple atomic sorts
*-------------------------------------------------------------------------------
	subroutine stoddard(step)

	implicit none

	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'cvtpp.inc'
	include 'boxpp.inc'
	include 'timep.inc'

	integer step
	integer i,j,k,m,catoms
	integer ind,lind
	integer nliqmon
	integer cvtnt(mxnlist),cvtii(mxa+1)
	integer cvtn(mxnlist),cvti(mxa+1),cvtl(mxa+1),liqmon(natom)
	double precision dx,dy,dz
	double precision rij2,rijc2
	double precision bf1,bf2
	double precision ibx,iby,ibz
	double precision tx,ty,tz

*-------------------------------------------------------------------------------
* step          : Time step number
*...............................................................................
* i,j,k		: standard loop variable
* ind,lind	: Index for neighborlist
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
	  clf(i)=0
	  ccf(i)=0
	  csf(i)=0
	  csif(i)=0
	  cekf(i)=0
	  ceksf(i)=0
	  cvtl(i)=0
	  cgs(i)=0
	  cgsf(i)=0
	  cgsa(i)=0
	  cgsfa(i)=0
	enddo
	csif(0) = 0
	cgsf(0) = 0
	cgsfa(0)= 0.0d0

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
	    dx=dx+boxx*(bf1-int(bf2+dx*ibx))
	    dy=dy+boxy*(bf1-int(bf2+dy*iby))
	    dz=dz+boxz*(bf1-int(bf2+dz*ibz))
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
	if (i.eq.natom) goto 500
520	j=cvtnt(ind)
	if (j.lt.0) goto 510
	ind=ind+1
	  cvtii(i)=cvtii(i)+1
	  cvtii(j)=cvtii(j)+1
	  cvtn(cvtii(j))=i
	  cvtn(cvtii(i))=j
	goto 520
500	continue

*** Building the cluster statistics

	lind = 1
	cnum = 0
	nliqmon = 0
	do k=1,natom

***  This part filters out the liquid monomers
	  if (cvtl(k).ge.atnenum) then
	    liqmon(lind) = k
	    nliqmon = nliqmon + 1
	    lind = lind + 1
	  Else
	    cekf(0) = cekf(0)+
     $	      0.5d0*apmas(api(k))*
     $	      (vx(k)*vx(k)+vy(k)*vy(k)+vz(k)*vz(k))
	  endif
	  cnum=cnum+1
	  j=k
	  if (cc(j).eq.0) then
200	    cl(j)=cnum
	    cc(j)=1
	    if (api(j).eq.1) Then
	       cs(cnum) = cs(cnum) + 1
	       cek(cnum)=cek(cnum)+
     $	      0.5d0*apmas(api(j))*
     $	      (vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j))
	    else
	       cgs(cnum) = cgs(cnum) + 1
	    endif

	    ind=cvti(j)+1
	    do i=ind,ind+cvtl(j)-1
	      cl(cvtn(i))=cnum
	    enddo

	    do j=k+1,natom
	      if(cc(j).eq.0.and.cl(j).eq.cnum) goto 200
	    enddo
	  endif
	enddo



*** Do the same for Frenkel definition
	cnumf = 0
	do k=1,nliqmon
	  cnumf = cnumf + 1
	  j=k
	  m = liqmon(j)
	  if (ccf(m).eq.0) then
300	    clf(m) = cnumf
	    ccf(m) = 1
	    if (api(m).eq.1) Then
	       csf(cnumf) = csf(cnumf) + 1
	       cekf(cnumf) = cekf(cnumf)+
     $	      0.5d0*apmas(api(m))*
     $	      (vx(m)*vx(m)+vy(m)*vy(m)+vz(m)*vz(m))
	    else
	       cgsf(cnumf) = cgsf(cnumf) + 1
	    endif

	    ind = cvti(m) + 1
	    do i = ind,ind+cvtl(m)-1
	      clf(cvtn(i)) = cnumf
	    enddo

	    do j = k+1,nliqmon
	      m = liqmon(j)
	      if(ccf(m).eq.0.and.clf(m).eq.cnumf) goto 300
	    enddo
	  endif
	enddo

*** Do the statistics

	j = 0
	do i=1,cnum
	  csi(cs(i)) = csi(cs(i)) + 1
	  ceks(cs(i)) = ceks(cs(i)) + cek(i)
	  cgsa(cs(i)) = cgsa(cs(i)) + cgs(i)
	  if (cs(i).gt.j) j = i
	enddo

	catoms = 0
	do i=1,cnumf
	  csif(csf(i)) = csif(csf(i)) + 1
	  ceksf(csf(i)) = ceksf(csf(i)) + cekf(i)
	  cgsfa(csf(i)) = cgsfa(csf(i)) + cgsf(i)
	  catoms = catoms + csf(i)
	enddo
	csif(0) = apnum(1) - catoms


***	First reset all atom labels to the default one
	do i=1,apnum(1)
	   apsymi(i) = apsym(api(1))
	enddo

***	Ok - now mark all the atoms belonging to Stillinger clusters >=2 with the 
*** symbol "S"
*** Remember: all atoms have a clustelabel (cl() or clf()) which is the number 
***(not the size! that's cs()!) of the cluster they belong to. 
	k = 0
	do i=1,apnum(1)
	   j = cl(i)
	   if (cs(j).gt.5) Then
		apsymi(i) = 'S'
		k = k + 1
	   endif
	enddo
*	write (*,*) 'Still:', k
*** Finally, we loop again over all particles and now mark those which are
*** liguid monomers according to ten-Wolde-Frenkel as "F"
	k = 0
	do i=1,nliqmon
	   j = liqmon(i)
	   apsymi(j) = 'F'
	   k = k + 1
	enddo
*	write (*,*) 'twf:', k

88	format (A2,3(2X,F10.4))

99	continue
	return
	end
*-------------------------------------------------------------------------------
