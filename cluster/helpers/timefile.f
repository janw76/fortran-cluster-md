	program Timefile

	include 'const.inc'
	include 'timep.inc'
	include 'filep.inc'

	integer i,maxi
	character*40 fout
	double precision t

	call readdata

	i=2
	call getarg(i,fout)

	open (1,file=fout)

	maxi=ntim/ftxyz

	write (*,*) dt,ftxyz

	do i=0,maxi-1
	  write (1,*) 'NAME TimeFile'
	  write (1,'(A,I10)') ' TIME',i
	  write (1,*) 'RANK 2'
	  write (1,*) 'DIMENSIONS 2 2'
	  write (1,*) 'BOUNDS 0 1000000000 0 1000000000'
	  write (1,*) 'SCALAR'
	  write (1,*) 'DATA'
	  t=real(i)*dt*real(ftxyz)+real(btim)*dt
	  write (1,'(2(2X,F18.8))') t/1000.0,t
	  write (1,'(2(2X,F18.8))') t/1000000.0,real(i)
	  write (1,*) 'END'
	  write (1,*)
	enddo

* 0,0: Time in picoseconds
* 1,0: Time in femtoseconds
* 1,0: Time in nanoseconds
* 1,0: Time in frames

	close (1)

	end
