*===============================================================================
*	ACCELERATION ATOM-ATOM (GPU Accelerated)
*===============================================================================
* file: for-a-lj-n-pbc-m-gpu.f
*-------------------------------------------------------------------------------
* GPU-accelerated version using OpenACC directives
* Experimental - for future GPU support on Apple Silicon
*-------------------------------------------------------------------------------

	subroutine AccAtom()
	
	implicit none

	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'energ.inc'
	include 'boxpp.inc'

	integer i,j,k,ij
	double precision atr2,rij2,rij6
	double precision fij,esum
	double precision dx,dy,dz,tx,ty,tz
	double precision c111,c112,c113
	double precision c221,c222,c223
	double precision c121,c122,c123
	double precision eps11,eps22,eps12
	double precision bf1,bf2
	double precision fsumx(mxa),fsumy(mxa),fsumz(mxa)
	double precision ibx,iby,ibz,iatm

*-------------------------------------------------------------------------------
* GPU acceleration using OpenACC:
* - Parallel loops on GPU
* - Optimized memory transfers
* - Vector operations
*-------------------------------------------------------------------------------

	atr2=atrcut*atrcut
	c111=apsig(1)*apsig(1)*apsig(1)*apsig(1)*apsig(1)*apsig(1)
	c112=c111*c111
	c113=0.5d0*c111
	c221=apsig(2)*apsig(2)*apsig(2)*apsig(2)*apsig(2)*apsig(2)
	c222=c221*c221
	c223=0.5d0*c221
	c121=0.5d0*(apsig(1)+apsig(2))
	c121=c121*c121*c121*c121*c121*c121
	c122=c121*c121
	c123=0.5d0*c121

	eps11=apeps(1)*apeps(1)
	eps22=apeps(2)*apeps(2)
	eps12=apeps(1)*apeps(2)

	bf1=int(atrskin/atr)+2.0d0
	bf2=bf1+0.5d0
	ibx=1.0d0/boxx
	iby=1.0d0/boxy
	ibz=1.0d0/boxz
	esum=0.0d0

	! Copy data to GPU and initialize force arrays
!$ACC DATA COPYIN(x(1:natom),y(1:natom),z(1:natom),atnlist,api(1:natom)) &
!$ACC      COPYOUT(fsumx(1:natom),fsumy(1:natom),fsumz(1:natom)) &
!$ACC      COPY(esum)

!$ACC KERNELS
!$ACC LOOP INDEPENDENT
	do i=1,natom
	  fsumx(i)=0.0d0
	  fsumy(i)=0.0d0
	  fsumz(i)=0.0d0
	enddo
!$ACC END KERNELS

	k=1
	
	! GPU-accelerated force calculation
!$ACC PARALLEL LOOP REDUCTION(+:esum) &
!$ACC PRIVATE(i,j,k,tx,ty,tz,dx,dy,dz,rij2,rij6,fij,ij)
10	do while (k.le.mxnlist)
	  i=-atnlist(k)
	  if (i.eq.natom) goto 30

	  k=k+1
	  tx=x(i)
	  ty=y(i)
  	  tz=z(i)

20	  j=atnlist(k)
	  if (j.gt.0) then
	    k=k+1

	    dx=tx-x(j)
	    dy=ty-y(j)
	    dz=tz-z(j)
	    dx=dx+boxx*(bf1-int(bf2+dx*ibx))
	    dy=dy+boxy*(bf1-int(bf2+dy*iby))
	    dz=dz+boxz*(bf1-int(bf2+dz*ibz))

	    rij2=dx*dx+dy*dy+dz*dz

	    if (rij2.le.atr2) then
	      ij=api(i)*10+api(j)
	      rij2=1.0d0/rij2
	      rij6=rij2*rij2*rij2
	      
	      if (ij.eq.11) then
	        fij=rij2*rij6*(c112*rij6-c113)*eps11
	        esum=esum+rij6*(c112*rij6-c111)*eps11
	      else if (ij.eq.22) then
	        fij=rij2*rij6*c222*eps22*12.0d0
	        esum=esum+rij6*c222*eps22
	      else
	        fij=rij2*rij6*c122*eps12*12.0d0
	        esum=esum+rij6*c122*eps12
	      endif
	      
!$ACC ATOMIC
	      fsumx(i)=fsumx(i)+fij*dx
!$ACC ATOMIC
	      fsumx(j)=fsumx(j)-fij*dx
!$ACC ATOMIC
	      fsumy(i)=fsumy(i)+fij*dy
!$ACC ATOMIC
	      fsumy(j)=fsumy(j)-fij*dy
!$ACC ATOMIC
	      fsumz(i)=fsumz(i)+fij*dz
!$ACC ATOMIC
	      fsumz(j)=fsumz(j)-fij*dz
	    endif
	    goto 20
	  else
	    goto 10
	  endif
	enddo
30	continue

!$ACC END DATA

	epot=esum*4.0d0
	
	! Convert forces to accelerations
	do i=1,natom
	  iatm=1.0d0/apmas(api(i))
	  ax(i)=ax(i)+48.0d0*fsumx(i)*iatm
	  ay(i)=ay(i)+48.0d0*fsumy(i)*iatm
	  az(i)=az(i)+48.0d0*fsumz(i)*iatm
	enddo

	return
	end
*-------------------------------------------------------------------------------
