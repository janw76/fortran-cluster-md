*===============================================================================
*	CONDENSATION AND EVAPORATION RATES ... Statistics
*-------------------------------------------------------------------------------
* file: ratstat.f
* date: 15-may-2001
* call: ratstat <input-file> <threshold>
*-------------------------------------------------------------------------------
* This program generates various statistics from the .rat-file of a
* Simulation. Multiple .rat-files of a single run can be concatenated and
* processed in one go!
* The threshold gives the minimum number of occurences of a cluster number
* to be included in the life-statistics
*-------------------------------------------------------------------------------

	program ratstat
	
	implicit none

	integer lblen,mxc,mxa,ncycle
	parameter (lblen=4000,mxc=20,mxa=1000,ncycle=10)
	integer i,j,io,cycle
	character*40 fout,fin
	character*40 arg
	character*100 inline
	integer step,omcs,mcs,cn,en,mcln,omcln
	integer lob,upb,dtime,steps,nd,ndd,thresh
	integer ssn,csn,esn
	double precision ssx,ssy,ssxx,ssyy,ssxy,sa,sb,sr,stmp
	double precision esx,esy,esxx,esyy,esxy,ea,eb,er,etmp
	double precision csx,csy,csxx,csyy,csxy,ca,cb,cr,ctmp
	double precision cao,cbo,cro,eao,ebo,ero,sao,sbo,sro
	double precision mct,omct
	integer grstat(mxa),grn(mxa),ostep
	integer lbstat(0:lblen)
	integer ebstat(0:lblen)
	integer cbstat(0:lblen)
	integer slbstat(0:mxa,0:lblen)
	integer elbstat(0:mxa,0:lblen)
	integer clbstat(0:mxa,0:lblen)
	integer ssum,sssum,esum,essum,csum,cssum
	real ssm,csm,esm
	integer sstat(0:mxa,-mxc:mxc+1)

	i=1
	call getarg(i,fin)

	i=2
	call getarg(i,arg)
	if (arg.eq.'') arg='0'
	read (arg,*) thresh

	do i=1,mxa
	  grstat(i)=0
	  grn(i)=0
	enddo
	do i=1,lblen
	  lbstat(i)=0
	enddo

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
	  if (dtime.le.lblen) lbstat(dtime)=lbstat(dtime)+1
	  if (dtime.gt.lblen) lbstat(0)=lbstat(0)+1
	  if (dtime.le.lblen.and.en.gt.0) ebstat(dtime)=ebstat(dtime)+1
	  if (dtime.gt.lblen.and.en.gt.0) ebstat(0)=ebstat(0)+1
	  if (dtime.le.lblen.and.cn.gt.0) cbstat(dtime)=cbstat(dtime)+1
	  if (dtime.gt.lblen.and.cn.gt.0) cbstat(0)=cbstat(0)+1
	  ostep=step
	goto 10
20	continue
	close (1)

	lob=0
	do i=1,mxa
	  if (grstat(i).ne.0.and.lob.eq.0) lob=max(0,i-20)
	enddo
	upb=0
	do i=mxa,1,-1
	  if (grstat(i).ne.0.and.upb.eq.0) upb=min(mxa,i+20)
	enddo

	open (2,file='sizes.dat')
	do i=lob,upb
	  write (2,101) i,grstat(i),grn(i)
	enddo
	close (2)
101	format (I5,2X,I9,2X,I7)

	open (2,file='life-tot.dat')
	do i=0,lblen
	  write (2,102) i,cbstat(i),ebstat(i),lbstat(i)
	enddo
	close (2)
102	format (I7,3(2X,I5))

	do i=0,lblen
	  do j=0,mxa
	    slbstat(j,i)=0
	  enddo
	enddo
	do i=-mxc,mxc+1
	  do j=0,mxa
	    sstat(j,i)=0
	  enddo
	enddo

	nd=0
	ndd=0
	open (1,file=fin)
	read(1,11) ostep,omcs,mcs,mct,omct,cn,en,mcln,omcln
	do i=1,steps
	  read(1,11) step,omcs,mcs,mct,omct,cn,en,mcln,omcln
	  dtime=step-ostep
	  if (cn*en.gt.0) nd=nd+1
	  if (cn*en.gt.1) ndd=ndd+1
	  if (dtime.le.lblen) then
	    slbstat(omcs,dtime)=slbstat(omcs,dtime)+1
	    slbstat(omcs,0)=slbstat(omcs,0)+1
	    if (en.gt.0) then
	      elbstat(omcs,dtime)=elbstat(omcs,dtime)+1
	      elbstat(omcs,0)=elbstat(omcs,0)+1
	    endif
	    if (cn.gt.0) then
	      clbstat(omcs,dtime)=clbstat(omcs,dtime)+1
	      clbstat(omcs,0)=clbstat(omcs,0)+1
	    endif
	  endif
	  if (en.gt.0.and.en.le.mxc) then
	    sstat(omcs,-en)=sstat(omcs,-en)+1
	    sstat(omcs,mxc+1)=sstat(omcs,mxc+1)+1
	    sstat(0,-en)=sstat(0,-en)+1
	  endif
	  if (cn.gt.0.and.cn.le.mxc) then
	    sstat(omcs,cn)=sstat(omcs,cn)+1
	    sstat(0,cn)=sstat(0,cn)+1
	  endif
	  sstat(omcs,0)=sstat(omcs,0)+dtime-1
	  sstat(omcs,mxc+1)=sstat(omcs,mxc+1)+1
	  sstat(0,0)=sstat(0,0)+dtime-1
	  ostep=step
	enddo
	close(1)
	write (*,'(A,I5,2X,I4)') 'Double Events: ',nd,ndd

	open (1,file='event-tot.dat')
	do j=-mxc,mxc
	  write (1,'(I4,2X,I9)') j,sstat(0,j)
	enddo
	close (1)
	
	do i=1,mxa
	  fout='event00000.dat'
	  write (fout,'(A5,I5.5,A4)') 'event',i,'.dat'
	  if (sstat(i,mxc+1).gt.0) then
	    open (2,file=fout)
	    do j=-mxc,mxc
	      write (2,'(I4,2X,I6)') j,sstat(i,j)
	    enddo
	    close (2)
	  endif
	enddo

	open (1,file='life-mean.dat')
	do i=0,mxa
	  fout='life00000.dat'
	  write (fout,'(A4,I5.5,A4)') 'life',i,'.dat'
	  if (slbstat(i,0).gt.thresh) then
	    sro=0.0d0
	    ero=0.0d0
	    cro=0.0d0
	    sao=-10.0d0
	    eao=-10.0d0
	    cao=-2.0d0
	    sbo=0.0d0
	    ebo=0.0d0
	    cbo=0.0d0
	    cycle=ncycle
	    open (2,file=fout)
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
	      if (cycle.eq.ncycle) then
	 	write (2,103) j,clbstat(i,j),elbstat(i,j),slbstat(i,j)
	      endif
103	format (I7,3(2X,I5))
	      sssum=sssum+j*slbstat(i,j)
	      ssum=ssum+slbstat(i,j)
	      essum=essum+j*elbstat(i,j)
	      esum=esum+elbstat(i,j)
	      cssum=cssum+j*clbstat(i,j)
	      csum=csum+clbstat(i,j)
	      if(slbstat(i,j).gt.0.or.j.gt.lblen) then
		ssn=ssn+1
		if (slbstat(i,j).eq.0) then 
		  stmp=sao+sbo*dble(j)
		else 
		  stmp=log(dble(slbstat(i,j)))
		endif
		ssx=ssx+dble(j)
		ssy=ssy+stmp
		ssxx=ssxx+dble(j)*dble(j)
		ssyy=ssyy+stmp*stmp
		ssxy=ssxy+dble(j)*stmp
	      endif
	      if(elbstat(i,j).gt.0.or.j.gt.lblen) then
		esn=esn+1
		if (elbstat(i,j).eq.0) then 
		  etmp=eao+ebo*dble(j)
		else 
		  etmp=log(dble(elbstat(i,j)))
		endif
		esx=esx+dble(j)
		esy=esy+etmp
		esxx=esxx+dble(j)*dble(j)
		esyy=esyy+etmp*etmp
		esxy=esxy+dble(j)*etmp
	      endif
	      if(clbstat(i,j).gt.0.or.j.gt.lblen) then
		csn=csn+1
		if (clbstat(i,j).eq.0) then 
		  ctmp=cao+cbo*dble(j)
		else 
		  ctmp=log(dble(clbstat(i,j)))
		endif
		csx=csx+dble(j)
		csy=csy+ctmp
		csxx=csxx+dble(j)*dble(j)
		csyy=csyy+ctmp*ctmp
		csxy=csxy+dble(j)*ctmp
	      endif
	    enddo
	    close (2)
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
	
	end
*-------------------------------------------------------------------------------
