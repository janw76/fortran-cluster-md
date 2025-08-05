*===============================================================================
*	CONDENSATION AND EVAPORATION RATES ... Statistics
*-------------------------------------------------------------------------------
* file: ratstat-t.f
* date: 15-may-2001
* call: ratstat <script-file> <threshold>
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
	integer ncycle,cycle
	parameter(ncycle=20)
	logical flag
	integer i,j,io,ind
	character*40 fout,fin,fout2,fin2
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
	integer ssn,csn,esn
	double precision ssx,ssy,ssxx,ssyy,ssxy,sa,sb,sr,stmp
	double precision esx,esy,esxx,esyy,esxy,ea,eb,er,etmp
	double precision csx,csy,csxx,csyy,csxy,ca,cb,cr,ctmp
	double precision cao,cbo,cro,eao,ebo,ero,sao,sbo,sro

	i=1
	call getarg(i,fin2)

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

	i=2
	call getarg(i,arg)
	if (arg.eq.'') arg='0'
	read (arg,*) thresh

	open (42,file=fin2)
4242	read (42,'(A)') fin
	if (fin(1:3).eq.'OUT') then
	  flag=.false.
	  if (fin(4:4).eq.'C') flag=.true.
	  close (42)
	  goto 999
	endif

	open (1,file=fin)
	write (*,*) 'Processing: ',fin

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
	goto 4242

999 	continue

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
	    sro=0.0d0
	    ero=0.0d0
	    cro=0.0d0
	    sao=-10.0d0
	    eao=-10.0d0
	    cao=-10.0d0
	    sbo=0.0d0
	    ebo=0.0d0
	    cbo=0.0d0
	    cycle=ncycle
	    if (flag) open (2,file=fout)
666	    ssum=0
	    sssum=0
	    esum=0
	    essum=0
	    csum=0
	    cssum=0
	    ssn=0
	    ssx=0.0d0
	    ssy=0.0d0
	    ssxx=0.0d0
	    ssyy=0.0d0
	    ssxy=0.0d0
	    sa=0.0d0
	    sb=0.0d0
	    esn=0
	    esx=0.0d0
	    esy=0.0d0
	    esxx=0.0d0
	    esyy=0.0d0
	    esxy=0.0d0
	    ea=0.0d0
	    eb=0.0d0
	    csn=0
	    csx=0.0d0
	    csy=0.0d0
	    csxx=0.0d0
	    csyy=0.0d0
	    csxy=0.0d0
	    ca=0.0d0
	    cb=0.0d0
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
	      if(slbstat(i,j,io).gt.0.or.j.gt.lblen) then
		ssn=ssn+1
		if (slbstat(i,j,io).eq.0) then 
		  stmp=sao+sbo*dble(j)
		else 
		  stmp=log(dble(slbstat(i,j,io)))
		endif
		ssx=ssx+dble(j)
		ssy=ssy+stmp
		ssxx=ssxx+dble(j)*dble(j)
		ssyy=ssyy+stmp*stmp
		ssxy=ssxy+dble(j)*stmp
	      endif
	      if(elbstat(i,j,io).gt.0.or.j.gt.lblen) then
		esn=esn+1
		if (elbstat(i,j,io).eq.0) then 
		  etmp=eao+ebo*dble(j)
		else 
		  etmp=log(dble(elbstat(i,j,io)))
		endif
		esx=esx+dble(j)
		esy=esy+etmp
		esxx=esxx+dble(j)*dble(j)
		esyy=esyy+etmp*etmp
		esxy=esxy+dble(j)*etmp
	      endif
	      if(clbstat(i,j,io).gt.0.or.j.gt.lblen) then
		csn=csn+1
		if (clbstat(i,j,io).eq.0) then 
		  ctmp=cao+cbo*dble(j)
		else 
		  ctmp=log(dble(clbstat(i,j,io)))
		endif
		csx=csx+dble(j)
		csy=csy+ctmp
		csxx=csxx+dble(j)*dble(j)
		csyy=csyy+ctmp*ctmp
		csxy=csxy+dble(j)*ctmp
	      endif
	    enddo
	    if (flag) close (2)
	    ssm=real(sssum)/real(ssum)
	    csm=real(cssum)/real(csum)
	    esm=real(essum)/real(esum)
	    if (ssn.gt.0) then
	      ssx=ssx/dble(ssn)
	      ssy=ssy/dble(ssn)
	      ssxx=ssxx/dble(ssn)
	      ssyy=ssyy/dble(ssn)
	      ssxy=ssxy/dble(ssn)
	      sb=(ssxy-ssx*ssy)/(ssxx-ssx*ssx)
	      sa=ssy-sb*ssx
	      sr=(ssxy-ssx*ssy)**2/((ssxx-ssx*ssx)*(ssyy-ssy*ssy))
	      sao=sa
	      sbo=sb
	      sro=sr
	      sb=-1.0d0/sb
	      sa=exp(sa)
	    endif
	    if (esn.gt.0) then
	      esx=esx/dble(esn)
	      esy=esy/dble(esn)
	      esxx=esxx/dble(esn)
	      esyy=esyy/dble(esn)
	      esxy=esxy/dble(esn)
	      eb=(esxy-esx*esy)/(esxx-esx*esx)
	      ea=esy-eb*esx
	      er=(esxy-esx*esy)**2/((esxx-esx*esx)*(esyy-esy*esy))
	      eao=ea
	      ebo=eb
	      ero=er
	      eb=-1.0d0/eb
	      ea=exp(ea)
	    endif
	    if (csn.gt.0) then
	      csx=csx/dble(csn)
	      csy=csy/dble(csn)
	      csxx=csxx/dble(csn)
	      csyy=csyy/dble(csn)
	      csxy=csxy/dble(csn)
	      cb=(csxy-csx*csy)/(csxx-csx*csx)
	      ca=csy-cb*csx
	      cr=(csxy-csx*csy)**2/((csxx-csx*csx)*(csyy-csy*csy))
	      cao=ca
	      cbo=cb
	      cro=cr
	      cb=-1.0d0/cb
	      ca=exp(ca)
	    endif
	    if (cycle.gt.0) then
	      cycle=cycle-1
	      goto 666
	    endif
	    write (1,104) i,csm,esm,ssm,cb,eb,sb,ca,ea,sa,cr,er,sr
104	format (I7,9(1X,F10.4),3(1X,f7.5))
	  endif
	enddo
	close (1)
	enddo
	
	end
*-------------------------------------------------------------------------------
