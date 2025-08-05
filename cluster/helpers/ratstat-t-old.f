*===============================================================================
*	CONDENSATION AND EVAPORATION RATES ... Statistics
*-------------------------------------------------------------------------------
* file: ratstat-t.f
* date: 23-mar-2001
* call: ratstat <input-file> <threshold>
*-------------------------------------------------------------------------------
* This program generates various statistics from the .rat-file of a
* Simulation. Multiple .rat-files of a single run can be concatenated and
* processed in one go!
* The threshold gives the minimum number of occurences of a cluster number
* to be included in the life-statistics
*-------------------------------------------------------------------------------
* This program uses *BIG* fields, set 'ulimit -s 500000' (bash), 
* or limit stacksize 500000 (sh)
*-------------------------------------------------------------------------------


	program ratstat
	
	implicit none

	integer lblen,mxc,mts,mxal,mxau,mxa
	double precision btl,btu,dts
	parameter (lblen=3000,mxc=20,mxal=70,mxau=900,mxa=1000)
	parameter (btl=78.0d0,btu=88.0d0,dts=1.0d0)
	parameter (mts=int(real(btu-btl)/dts))

	logical flag
	integer i,j,io,ind
	character*40 fout,fin,fout2
	character*40 arg
	character*100 inline
	integer step,omcs,mcs,cn,en,mcln,omcln
	integer dtime,steps,thresh
	double precision mct,omct
	integer grstat(mxa),grn(mxa),ostep
	integer grtstat(mxa,0:mts),grtn(mxa,0:mts)
	integer slbstat(mxal:mxau,0:lblen,0:mts)
	integer elbstat(mxal:mxau,0:lblen,0:mts)
	integer clbstat(mxal:mxau,0:lblen,0:mts)
	integer ssum,sssum,esum,essum,csum,cssum
	real ssm,csm,esm

	i=1
	call getarg(i,fin)

	if (fin(1:4).eq.'init') then
	  write (*,*) 'Writing save file'
	  write (*,*) 'Size:', lblen+1,mts+1
	  do io=0,mts
	    do i=mxal,mxau
	      do j=0,lblen
	        slbstat(i,j,io)=0
	        elbstat(i,j,io)=0
	        clbstat(i,j,io)=0
	      enddo
	    enddo
	  enddo
	  do i=1,mxa
	    grstat(i)=0
	    grn(i)=0
	  enddo
	  do io=0,mts
	    do i=1,mxa
	      grtstat(i,io)=0
	      grtn(i,io)=0
	    enddo
	  enddo
	  open (1,file='rat-slb.dat',form='unformatted')
	  write (1) slbstat
	  close (1)
	  open (1,file='rat-elb.dat',form='unformatted')
	  write (1) elbstat
	  close (1)
	  open (1,file='rat-clb.dat',form='unformatted')
	  write (1) clbstat
	  close (1)
	  open (1,file='rat-gr.dat',form='unformatted')
	  write (1) grstat,grn,grtstat,grtn
	  close (1)
	  stop
	else
	  open (1,file='rat-slb.dat',form='unformatted')
	  read (1) slbstat
	  close (1)
	  open (1,file='rat-elb.dat',form='unformatted')
	  read (1) elbstat
	  close (1)
	  open (1,file='rat-clb.dat',form='unformatted')
	  read (1) clbstat
	  close (1)
	  open (1,file='rat-gr.dat',form='unformatted')
	  read (1) grstat,grn,grtstat,grtn
	  close (1)
	endif

	if (fin(1:3).eq.'out') then
	  goto 999
	endif

	i=2
	call getarg(i,arg)
	if (arg.eq.'') arg='0'
	read (arg,*) thresh

	open (1,file=fin)

	steps=0
	read(1,11,iostat=io) ostep,omcs,mcs,mct,omct,cn,en,mcln,omcln
	if (io.ne.0) stop
10	continue
	  read(1,'(A)',iostat=io) inline
	  read(inline,11) step,omcs,mcs,mct,omct,cn,en,mcln,omcln
	  if (io.ne.0.or.omcln.eq.0) goto 20
	  steps=steps+1
11	  format (I12,2(X,I6),2(X,F10.4),2(X,I3),2(X,I6))
	  dtime=step-ostep
	  grstat(omcs)=grstat(omcs)+dtime
	  grn(omcs)=grn(omcs)+1

	  ind=int((omct-btl+0.5d0*dts)/dts)
	  if (ind.ge.0.and.ind.le.mts) then
	    grtstat(omcs,ind)=grtstat(omcs,ind)+dtime
	    grtn(omcs,ind)=grtn(omcs,ind)+1
	  endif

	  ostep=step
	goto 10
20	continue
	close (1)

	open (1,file=fin)
	read(1,11) ostep,omcs,mcs,mct,omct,cn,en,mcln,omcln
	do i=1,steps
	  read(1,11) step,omcs,mcs,mct,omct,cn,en,mcln,omcln
	  dtime=step-ostep
	  ind=int((omct-btl+0.5d0*dts)/dts)
	if (ind.ge.0.and.ind.le.mts.
     $	   and.omcs.ge.mxal.and.omcs.le.mxau) then
	  if (dtime.le.lblen) then
	    slbstat(omcs,dtime,ind)=slbstat(omcs,dtime,ind)+1
	    slbstat(omcs,0,ind)=slbstat(omcs,0,ind)+1
	    if (en.gt.0) then
	      elbstat(omcs,dtime,ind)=elbstat(omcs,dtime,ind)+1
	      elbstat(omcs,0,ind)=elbstat(omcs,0,ind)+1
	    endif
	    if (cn.gt.0) then
	      clbstat(omcs,dtime,ind)=clbstat(omcs,dtime,ind)+1
	      clbstat(omcs,0,ind)=clbstat(omcs,0,ind)+1
	    endif
	  endif
	endif
	  ostep=step
	enddo
	close(1)

	  open (1,file='rat-slb.dat',form='unformatted')
	  write (1) slbstat
	  close (1)
	  open (1,file='rat-elb.dat',form='unformatted')
	  write (1) elbstat
	  close (1)
	  open (1,file='rat-clb.dat',form='unformatted')
	  write (1) clbstat
	  close (1)
	  open (1,file='rat-gr.dat',form='unformatted')
	  write (1) grstat,grn,grtstat,grtn

	goto 9999

999	continue

	flag=.false.
	if (fin(4:4).eq.'c') flag=.true.

	open (2,file='sizes-total.dat')
	do i=0,mxa
	  write (2,101) i,grstat(i),grn(i)
	enddo
	close (2)

	do io=0,mts
	  fout2='sizes-0000.dat'
	  write (fout2,'(A6,I4.4,A4)')
     $		      'sizes-',10*int(dble(io)*dts+btl+.01),'.dat'	
	  open (2,file=fout2)
	  do i=0,mxa
	    write (2,101) i,grtstat(i,io),grtn(i,io)
	  enddo
	  close (2)
	enddo
101	format (I5,2X,I9,2X,I7)

	do io=0,mts
	fout2='lm-0000.dat'
	write (fout2,'(A3,I4.4,A4)')
     $		      'lm-',10*int(dble(io)*dts+btl+.01),'.dat'
	open (1,file=fout2)
	do i=mxal,mxau
	  fout='life0000-00000.dat'
	  write (fout,'(A4,I4.4,A1,I5.5,A4)') 
     $		      'life',10*int(dble(io)*dts+btl+.01),'-',i,'.dat'
	  if (slbstat(i,0,io).gt.thresh) then
	    if (flag) open (2,file=fout)
	    ssum=0
	    sssum=0
	    esum=0
	    essum=0
	    csum=0
	    cssum=0
	    do j=1,lblen
	      if (flag) write (2,103)
     $		j,clbstat(i,j,io),elbstat(i,j,io),slbstat(i,j,io)
103	format (I7,3(2X,I5))
	      sssum=sssum+j*slbstat(i,j,io)
	      ssum=ssum+slbstat(i,j,io)
	      essum=essum+j*elbstat(i,j,io)
	      esum=esum+elbstat(i,j,io)
	      cssum=cssum+j*clbstat(i,j,io)
	      csum=csum+clbstat(i,j,io)
	    enddo
	    if (flag) close (2)
	    ssm=real(sssum)/real(ssum)
	    csm=real(cssum)/real(csum)
	    esm=real(essum)/real(esum)
	    write (1,104) i,csm,esm,ssm
104	format (I7,3(2X,F9.4))
	  endif
	enddo
	close (1)
	enddo
	
9999 	continue
	end
*-------------------------------------------------------------------------------
