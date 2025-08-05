	program veltest

	implicit none

	include 'timep.inc'
	include 'filep.inc'

	integer i,j
	integer steps,ind,nat
	double precision vx,vy,vz
	integer stat(-10000:10000)
	integer vstat(0:10000)
	character*40 ofile,par

	i=2
	call getarg(i,ofile)

	call readdata()

	steps=ntim/ftvel
	write (*,*) ofile
	write (*,*) steps
	

	do i=0,10000
	  stat(i)=0
	  stat(-i)=0
	  vstat(i)=0
	enddo

	open (1,file=fvel)

	do i=1,steps
	  read(1,*) nat
	  read(1,*) ind
	  do j=1,nat
	    read (1,'(A2,3(3X,F15.10))') par,vx,vy,vz
	    vx=1.0d5*vx
	    vy=1.0d5*vy
	    vz=1.0d5*vz
            ind=int(vx)
	    if (ind.gt.10000) write (*,*) 'OOOOPPPPPSSS!  X'
	    stat(ind)=stat(ind)+1
	    ind=int(vy)
	    if (ind.gt.10000) write (*,*) 'OOOOPPPPPSSS!  Y'
	    stat(ind)=stat(ind)+1
	    ind=int(vz)
	    if (ind.gt.10000) write (*,*) 'OOOOPPPPPSSS!  Z'
	    stat(ind)=stat(ind)+1
	    ind=int(dsqrt(vx*vx+vy*vy+vz*vz))
	    if (ind.gt.10000) write (*,*) 'Hurza!'
	    vstat(ind)=vstat(ind)+1
	  enddo
	enddo
	
	close(1)

	open (1,file=ofile)
	do i=0,10000
	  write (1,*) i,stat(i),(dble(vstat(i))/dble((nat*steps)))
	enddo

	close(1)

	end

