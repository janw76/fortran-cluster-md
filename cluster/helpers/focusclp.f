*------------------------------------------------------------------------------
* Programm zum Focussieren auf einen Cluster
*------------------------------------------------------------------------------

	program focusclp

	implicit none
	
	include 'const.inc'
	include 'atomp.inc'
	include 'atomc.inc'
	include 'timep.inc'
	include 'boxpp.inc'
	include 'filep.inc'
	include 'cvtpp.inc'

	character*40 fout,dummy
	integer i,j,k,step,maxi,mcs,disc,nnat
	double precision hx,hy,hz,cx,cy,cz,dx,dy,dz,hbx,hby,hbz
	double precision sbox,sbpx,sbpy,sbpz,sbmx,sbmy,sbmz
	double precision ibx,iby,ibz

	i=2
	call getarg(i,fout)
	i=3
	call getarg(i,dummy)
	read (dummy,*) disc
	i=4
	call getarg(i,dummy)
	read (dummy,*) sbox
	sbox=sbox*5.0d0

	call readdata

	open (1,file=fxyz)
	open (2,file=fout)

	maxi=ntim/ftxyz

	ibx=1.0d0/boxx
	iby=1.0d0/boxy
	ibz=1.0d0/boxz

	hbx=boxx/2.0d0
	hby=boxy/2.0d0
	hbz=boxz/2.0d0

	do i=1,maxi
	  call readXYZ(step,1)

	if (i.gt.disc) then

	  call stoddard()

	  mcs=1
	  k=1
	  do j=1,cnum
	    if (cs(j).gt.k) then
	      mcs=j
	      k=cs(j)
	    endif
	  enddo

	  cx=0.0d0
	  cy=0.0d0
	  cz=0.0d0

	  k=1
	  do j=1,natom
	    if (cl(j).eq.mcs) then
	      if (k.eq.1) then
		k=0
		hx=x(j)
		hy=y(j)
		hz=z(j)
	      endif
	      dx=x(j)-hx
	      dy=y(j)-hy
	      dz=z(j)-hz
	      cx=cx+dx+boxx-boxx*int(1.5d0+dx*ibx)
	      cy=cy+dy+boxy-boxy*int(1.5d0+dy*iby)
	      cz=cz+dz+boxz-boxz*int(1.5d0+dz*ibz)
	    endif
	  enddo
	  cx=dmod(dmod(hx+cx/real(cs(mcs)-1),boxx)+boxx,boxx)
	  cy=dmod(dmod(hy+cy/real(cs(mcs)-1),boxy)+boxy,boxy)
	  cz=dmod(dmod(hz+cz/real(cs(mcs)-1),boxz)+boxz,boxz)

	  dx=cx-hbx
	  dy=cy-hby
	  dz=cz-hbz

	write (*,*) cx,cy,cz
	write (*,*) dx,dy,dz


	  sbpx=hbx+sbox
	  sbpy=hby+sbox
	  sbpz=hbz+sbox
	  sbmx=hbx-sbox
	  sbmy=hby-sbox
	  sbmz=hbz-sbox
	  nnat=0
	  do j=1,natom
	    x(j)=dmod(dmod(x(j)-dx,boxx)+boxx,boxx)	    
	    y(j)=dmod(dmod(y(j)-dy,boxy)+boxy,boxy)	    
	    z(j)=dmod(dmod(z(j)-dz,boxz)+boxz,boxz)	    
	    if (x(j).ge.sbmx .and. x(j).le.sbpx .and. 
     $	        y(j).ge.sbmy .and. y(j).le.sbpy .and.
     $	        z(j).ge.sbmz .and. z(j).le.sbpz) then
	      nnat=nnat+1
	    endif 
	  enddo

	  if (bcor.ne.'A'.and.step.ne.btim) then
            write (2,*) nnat 
            write (2,*) step
          elseif (bcor.eq.'A'.or.(bcor.eq.'X'.and.step.eq.btim)) then
            write (2,*) nnat+8
            write (2,*) step 
	    write (2,10) 'H ',sbpx-hbx,sbpy-hby,sbpz-hbz
	    write (2,10) 'H ',sbmx-hbx,sbpy-hby,sbpz-hbz
	    write (2,10) 'H ',sbpx-hbx,sbmy-hby,sbpz-hbz
	    write (2,10) 'H ',sbmx-hbx,sbmy-hby,sbpz-hbz
	    write (2,10) 'H ',sbpx-hbx,sbpy-hby,sbmz-hbz
	    write (2,10) 'H ',sbmx-hbx,sbpy-hby,sbmz-hbz
	    write (2,10) 'H ',sbpx-hbx,sbmy-hby,sbmz-hbz
	    write (2,10) 'H ',sbmx-hbx,sbmy-hby,sbmz-hbz
          endif 

	  do j=1,natom
	    if (x(j).ge.sbmx .and. x(j).le.sbpx .and. 
     $	        y(j).ge.sbmy .and. y(j).le.sbpy .and.
     $	        z(j).ge.sbmz .and. z(j).le.sbpz) then
	      write (2,10) apsym(api(j)),x(j)-hbx,y(j)-hby,z(j)-hbz	  
	    endif
	  enddo
10        format (A2,3(2X,F10.4))

	endif

	  write (*,*) i,' of',maxi

	enddo

	close(1)
	close(2)

	end

