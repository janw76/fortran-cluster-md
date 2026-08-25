*===============================================================================
* 	NEIGHBOR LIST - PBC (OpenMP Optimized)
*===============================================================================
* file: neighbor-pbc-mav-omp.f
*-------------------------------------------------------------------------------
* OpenMP parallelized neighbor list construction
* Optimized for Apple Silicon with memory-efficient algorithms
*-------------------------------------------------------------------------------

	subroutine neighbor(step)

	implicit none

	include 'const.inc'
	include 'atomc.inc'
	include 'atomp.inc'
	include 'boxpp.inc'
	include 'timep.inc'
	include 'dmmpp.inc'

	integer step
	integer i,j,ind,chunk_size
	double precision tx,ty,tz,dx,dy,dz,r2,rs2
	double precision ibx,iby,ibz
	
	! Thread-private variables for neighbor list construction
	integer, allocatable :: local_neighbors(:,:)
	integer, allocatable :: local_counts(:)
	integer max_neighbors_per_atom
	parameter (max_neighbors_per_atom = 1024)

*-------------------------------------------------------------------------------
* Apple Silicon optimizations:
* - Parallel boundary condition wrapping
* - Thread-local neighbor lists to avoid contention
* - Cache-friendly memory access patterns
*-------------------------------------------------------------------------------

	if (step.eq.0) then
	  write (*,*) 'Determining neighbors for the first time ...'
	else
	  write (*,'(A,I10,A)') 'Step ',step,': Updating neighborlist'
	endif

	! Boundary condition wrapping in parallel
	ibx=1.0d0/boxx
	iby=1.0d0/boxy
	ibz=1.0d0/boxz

	! Serial parity: skin-exceeded fatal check, dmmav bookkeeping and
	! adaptive atrskin update (ported from neighbor-pbc-mav.f, same
	! placement relative to ibx/iby/ibz and rs2 so behavior matches).
	if (dmms.gt.(atrskin-atrcut)) then
	  write (*,*) '*** FATAL ERROR! ***  Skin radius exceeded!!'
	  stop
	endif

	dmmav=real(dmmst)*dmmsu/256.0d0

	if (step.ne.btim.and.step.gt.0) then
	  atrskin=atrcut+(atrskin-atrcut)*dble(nbtim)/dble(dmmst)
	  atrskin=min(boxx/2.0d0,atrskin)
	  atrskin=min(boxy/2.0d0,atrskin)
	  atrskin=min(boxz/2.0d0,atrskin)
	endif

	if (step.lt.0) step=-step
	
	rs2=atrskin*atrskin

	! Parallel boundary wrapping
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
	do i=1,natom
	  x(i)=dmod(dmod(x(i),boxx)+boxx,boxx)
	  y(i)=dmod(dmod(y(i),boxy)+boxy,boxy)
	  z(i)=dmod(dmod(z(i),boxz)+boxz,boxz)
	enddo
!$OMP END PARALLEL DO

	! Initialize neighbor list
	ind=1

	! Allocate shared per-atom scratch storage (each atom i is written
	! only by the single iteration that owns it, so no race despite the
	! arrays being shared across threads). Must be allocated/zeroed
	! single-threaded, before the parallel region, so all threads see
	! the same storage.
	allocate(local_neighbors(max_neighbors_per_atom, natom))
	allocate(local_counts(natom))
	local_counts = 0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,tx,ty,tz,dx,dy,dz,r2)

	! Parallel neighbor detection with spatial decomposition
!$OMP DO SCHEDULE(DYNAMIC,32)
	do i=1,natom
	  tx=x(i)
	  ty=y(i)
	  tz=z(i)
	  
	  do j=i+1,natom
	    dx=tx-x(j)
	    dy=ty-y(j)
	    dz=tz-z(j)
	    
	    ! Minimum image convention
	    dx=dx+boxx-boxx*int(1.5d0+dx*ibx)
	    dy=dy+boxy-boxy*int(1.5d0+dy*iby)
	    dz=dz+boxz-boxz*int(1.5d0+dz*ibz)

	    r2=dx*dx+dy*dy+dz*dz
	    if (r2.le.rs2) then
	      local_counts(i) = local_counts(i) + 1
	      if (local_counts(i).le.max_neighbors_per_atom) then
	        local_neighbors(local_counts(i), i) = j
	      endif
	    endif
	  enddo
	enddo
!$OMP END DO

	! Merge thread-local results into global neighbor list
!$OMP SINGLE
	ind = 1
	do i=1,natom
	  atnlist(ind) = -i
	  atnidx(i) = ind
	  ind = ind + 1
	  
	  if (local_counts(i).gt.max_neighbors_per_atom) then
	    write (*,*) '*** FATAL ERROR *** '//
     $                  'Per-atom neighbor list is full!'
	    stop 1
	  endif

	  do j=1,local_counts(i)
	    if (ind.lt.mxnlist) then
	      atnlist(ind) = local_neighbors(j,i)
	      ind = ind + 1
	    else
	      write (*,*) '*** FATAL ERROR *** '//
     $                    'Neighborlist is nearly full!'
	      stop
	    endif
	  enddo
	enddo
!$OMP END SINGLE

!$OMP END PARALLEL

	deallocate(local_neighbors)
	deallocate(local_counts)

	! Update displacement tracking
	dmmax=max(dmms,dmmax)
	if (step.gt.0) then
	  write (*,10) step,ind-1-natom,dmms,dmmav,dmmax,
     $	             atrskin-atrcut, atrskin/atr
	endif
	dmms=0.0d0
	dmmst=0
10	format (2(I10,2X),5(F8.4,2X))

	return
	end

*-------------------------------------------------------------------------------
