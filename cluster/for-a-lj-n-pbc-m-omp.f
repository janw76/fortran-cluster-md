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

	integer i,j,k,ij
	integer istart(mxa),iend(mxa)
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
* OpenMP parallelization: the force walk is parallelized over atom i (see
* Pass 1/Pass 2 below); Newton's-third-law updates to fsum*(j) and the
* i/j collision case are protected with ATOMIC, esum via REDUCTION.
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

	eps11=apeps(1)
	eps22=apeps(2)
	eps12=dsqrt(apeps(1)*apeps(2))

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

	! Pass 1 (serial, cheap): walk the packed neighbor list once and record
	! each atom's segment [istart(i),iend(i)] in atnlist. atnlist is built
	! with a -i marker for every atom 1..natom (even ones with zero
	! neighbors), followed by its positive-j neighbor entries; the walk
	! must stop as soon as the marker for atom natom is reached without
	! reading into its (always-empty, possibly stale) segment, exactly as
	! the original single-threaded walk did. This keeps the neighbor-list
	! traversal itself serial (it is an inherent pointer-chase, not a valid
	! OMP worksharing loop) while making the actual force evaluation below
	! embarrassingly parallel over i.
	k=1
10	i=-atnlist(k)
	if (i.eq.natom) goto 15
	k=k+1
	istart(i)=k
20	j=atnlist(k)
	if (j.gt.0) then
	  k=k+1
	  goto 20
	endif
	iend(i)=k-1
	goto 10
15	continue

	! Pass 2 (parallel over i): each thread evaluates the pairs of one
	! atom i (segment istart(i)..iend(i), all j>i by construction). Forces
	! obey Newton's third law (fsum*(i) and fsum*(j) both updated per
	! pair), and different threads can legitimately target the same array
	! slot: atom i's own accumulation here races with another thread's
	! j-side update where j equals this i (from a smaller i' processed by
	! that thread). All six fsum* updates are therefore atomic, applied
	! pair-by-pair in exactly the same order as the original serial code
	! (no private per-i batching) so that OMP_NUM_THREADS=1 reproduces the
	! pre-parallelization floating-point result bit-for-bit; esum uses an
	! OMP reduction, which for a single thread is likewise an exact-order
	! accumulation into a zero-initialized accumulator.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,tx,ty,tz,dx,dy,dz,rij2,
!$OMP& rij6,ij,fij) SCHEDULE(STATIC) REDUCTION(+:esum)
	do i=1,natom-1
	  tx=x(i)
	  ty=y(i)
	  tz=z(i)
	  do k=istart(i),iend(i)
	    j=atnlist(k)

	    dx=tx-x(j)
	    dy=ty-y(j)
	    dz=tz-z(j)
	    dx=dx+boxx*(bf1-int(bf2+dx*ibx))
	    dy=dy+boxy*(bf1-int(bf2+dy*iby))
	    dz=dz+boxz*(bf1-int(bf2+dz*ibz))

	    rij2=dx*dx+dy*dy+dz*dz

	    if (rij2.le.atr2) then
	      ij=api(i)+api(j)
	      rij2=1.0d0/rij2
	      rij6=rij2*rij2*rij2

	      ! Species-specific interactions
	      if (ij.eq.2) then
	        fij=eps11*48.0d0*rij2*rij6*(c112*rij6-c113)
	        esum=esum+eps11*4.0d0*rij6*(c112*rij6-c111)
	      else if (ij.eq.4) then
	        fij=eps22*48.0d0*rij2*rij6*(c222*rij6-c223)
	        esum=esum+eps22*4.0d0*rij6*(c222*rij6-c221)
	      else if (ij.eq.3) then
	        fij=eps12*48.0d0*rij2*rij6*(c122*rij6-c123)
	        esum=esum+eps12*4.0d0*rij6*(c122*rij6-c121)
	      else
	        write (*,*) 'Nanu???',i,j
	      endif

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
	  enddo
	enddo
!$OMP END PARALLEL DO

	epot=esum

	! Convert forces to accelerations in parallel
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,iatm) SCHEDULE(STATIC)
	do i=1,natom
	  iatm=1.0d0/apmas(api(i))
	  ax(i)=ax(i)+fsumx(i)*iatm
	  ay(i)=ay(i)+fsumy(i)*iatm
	  az(i)=az(i)+fsumz(i)*iatm
	enddo
!$OMP END PARALLEL DO

	return
	end
*-------------------------------------------------------------------------------
