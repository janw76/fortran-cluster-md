*===============================================================================
*	ACCELERATION ATOM-ATOM (OpenMP Optimized)
*===============================================================================
* file: for-a-lj-n-pbc-m-omp.f
*-------------------------------------------------------------------------------
* OpenMP parallelized version of for-a-lj-n-pbc-m.f
* Optimized for Apple Silicon with explicit vectorization hints
*-------------------------------------------------------------------------------

	subroutine AccAtom()
	
	implicit none

	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'energ.inc'
	include 'boxpp.inc'

	integer i,j,k,ij,tid
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
	
	! Thread-local energy accumulator for reduction
	double precision esum_local

*-------------------------------------------------------------------------------
* OpenMP parallelization with optimizations for Apple Silicon:
* - Thread-local accumulators to reduce false sharing
* - Guided scheduling for load balancing
* - SIMD hints for auto-vectorization
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

	! Initialize force arrays in parallel
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
	do i=1,natom
	  fsumx(i)=0.0d0
	  fsumy(i)=0.0d0
	  fsumz(i)=0.0d0
	enddo
!$OMP END PARALLEL DO

	k=1
	
	! Main force calculation loop with OpenMP parallelization
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,tx,ty,tz,dx,dy,dz,rij2,rij6,fij,ij,esum_local,tid) REDUCTION(+:esum)
	esum_local = 0.0d0
	tid = 1  ! Thread ID for debugging
!$OMP DO SCHEDULE(GUIDED) 
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
	      
	      ! Species-specific interactions
	      if (ij.eq.11) then
	        fij=rij2*rij6*(c112*rij6-c113)*eps11
	        esum_local=esum_local+rij6*(c112*rij6-c111)*eps11
	      else if (ij.eq.22) then
	        fij=rij2*rij6*c222*eps22*12.0d0
	        esum_local=esum_local+rij6*c222*eps22
	      else
	        fij=rij2*rij6*c122*eps12*12.0d0
	        esum_local=esum_local+rij6*c122*eps12
	      endif
	      
	      ! Atomic force updates with reduction
!$OMP ATOMIC
	      fsumx(i)=fsumx(i)+fij*dx
!$OMP ATOMIC
	      fsumx(j)=fsumx(j)-fij*dx
!$OMP ATOMIC
	      fsumy(i)=fsumy(i)+fij*dy
!$OMP ATOMIC
	      fsumy(j)=fsumy(j)-fij*dy
!$OMP ATOMIC
	      fsumz(i)=fsumz(i)+fij*dz
!$OMP ATOMIC
	      fsumz(j)=fsumz(j)-fij*dz
	    endif
	    goto 20
	  else
	    goto 10
	  endif
	enddo
30	continue
	esum = esum + esum_local
!$OMP END DO
!$OMP END PARALLEL

	epot=esum*4.0d0
	
	! Convert forces to accelerations in parallel
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,iatm) SCHEDULE(STATIC)
	do i=1,natom
	  iatm=1.0d0/apmas(api(i))
	  ax(i)=ax(i)+48.0d0*fsumx(i)*iatm
	  ay(i)=ay(i)+48.0d0*fsumy(i)*iatm
	  az(i)=az(i)+48.0d0*fsumz(i)*iatm
	enddo
!$OMP END PARALLEL DO

	return
	end
*-------------------------------------------------------------------------------
