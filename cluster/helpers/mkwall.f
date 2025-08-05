	program makewall

	implicit none

	integer i,j
	integer sy,dy,sx
	real r,x,y,ap
	character*40 dummy

	i=1
	call getarg(i,dummy)
	read(dummy,*) r
	i=2
	call getarg(i,dummy)
	read(dummy,*) dy
	i=3
	call getarg(i,dummy)
	read(dummy,*) sy

	sx=2*sy

	open (1,file='wall.nff')

	write (1,*) 'f 1.0 1.0 1.0 0 0 0 0 0'

	do i=-sy,sy,dy
	  if (i*i.ge.r*r) then
	    write (1,*) 'pl 2'
	    write (1,'(3(1X,I5))') 0,i,-sy
	    write (1,'(3(1X,I5))') 0,i, sy
	    write (1,*) 'pl 2'
	    write (1,'(3(1X,I5))') 0,-sy,i
	    write (1,'(3(1X,I5))') 0, sy,i
	  else
	    ap=sqrt(abs(r*r-real(i*i)))
	    write (1,*) 'pl 2'
	    write (1,'(3(1X,I5))') 0,i,-sy
	    write (1,'(3(1X,F10.5))') 0,real(i),-ap
	    write (1,*) 'pl 2'
	    write (1,'(3(1X,F10.5))') 0,real(i), ap
	    write (1,'(3(1X,I5))') 0,i, sy
	    write (1,*) 'pl 2'
	    write (1,'(3(1X,I5))') 0,-sy,i
	    write (1,'(3(1X,F10.5))') 0,-ap,real(i)
	    write (1,*) 'pl 2'
	    write (1,'(3(1X,F10.5))') 0, ap,real(i)
	    write (1,'(3(1X,I5))') 0, sy,i
	  endif
	enddo

	if (r.gt.0.0) then
	  write (1,*) 'pl 31'
	  do i=30,0,-1
	    ap=real(i)*8.0*atan(1.0)/30.0
	    x=r*sin(ap)
	    y=r*cos(ap)
	    write (1,'(3(1X,F10.5))') 0,x,y
	  enddo
	endif

*	goto 10

	do i=-sy,sy,dy
	  write(1,*) 'pl 2'
	  write (1,'(3(1X,I5))') -sx,-sy,i
	  write (1,'(3(1X,I5))')  sx,-sy,i
	  write(1,*) 'pl 2'
	  write (1,'(3(1X,I5))') -sx,i,-sy
	  write (1,'(3(1X,I5))')  sx,i,-sy
	  write(1,*) 'pl 2'
	  write (1,'(3(1X,I5))') -sx,i,-sy
	  write (1,'(3(1X,I5))') -sx,i, sy
	  write(1,*) 'pl 2'
	  write (1,'(3(1X,I5))') -sx,-sy,i
	  write (1,'(3(1X,I5))') -sx, sy,i
	enddo

	do i=-sx,sx,dy
	  write(1,*) 'pl 2'
	  write (1,'(3(1X,I5))') i,-sy,-sy
	  write (1,'(3(1X,I5))') i, sy,-sy
	  write(1,*) 'pl 2'
	  write (1,'(3(1X,I5))') i,-sy, sy
	  write (1,'(3(1X,I5))') i,-sy,-sy
	enddo

	write (1,*) 'f 0.5 0.5 0.5 0 0 0 0 0'
	write (1,*) 'p 4'
	write (1,'(3(1X,I5))') -sx,-sy,-sy
	write (1,'(3(1X,I5))')  sx,-sy,-sy
	write (1,'(3(1X,I5))')  sx, sy,-sy
	write (1,'(3(1X,I5))') -sx, sy,-sy
	write (1,*) 'p 4'
	write (1,'(3(1X,I5))') -sx,-sy,-sy
	write (1,'(3(1X,I5))') -sx,-sy, sy
	write (1,'(3(1X,I5))') -sx, sy, sy
	write (1,'(3(1X,I5))') -sx, sy,-sy
	write (1,*) 'p 4'
	write (1,'(3(1X,I5))') -sx,-sy,-sy
	write (1,'(3(1X,I5))')  sx,-sy,-sy
	write (1,'(3(1X,I5))')  sx,-sy, sy
	write (1,'(3(1X,I5))') -sx,-sy, sy

	write (1,*) 'f 1.0 1.0 1.0 0 0 0 0 0'
	write (1,*) 'pl 5'
	write (1,'(3(1X,I5))') -sx,-sy,-sy
	write (1,'(3(1X,I5))')  sx,-sy,-sy
	write (1,'(3(1X,I5))')  sx, sy,-sy
	write (1,'(3(1X,I5))') -sx, sy,-sy
	write (1,'(3(1X,I5))') -sx,-sy,-sy
	write (1,*) 'pl 5'
	write (1,'(3(1X,I5))') -sx,-sy,-sy
	write (1,'(3(1X,I5))') -sx, sy,-sy
	write (1,'(3(1X,I5))') -sx, sy, sy
	write (1,'(3(1X,I5))') -sx,-sy, sy
	write (1,'(3(1X,I5))') -sx,-sy,-sy
	write (1,*) 'pl 5'
	write (1,'(3(1X,I5))') -sx,-sy,-sy
	write (1,'(3(1X,I5))')  sx,-sy,-sy
	write (1,'(3(1X,I5))')  sx,-sy, sy
	write (1,'(3(1X,I5))') -sx,-sy, sy
	write (1,'(3(1X,I5))') -sx,-sy,-sy

10	continue

	end

