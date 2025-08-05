*===============================================================================
*	RANDOM NUMBER GENERATOR
*===============================================================================
* file: kps.f
*-------------------------------------------------------------------------------
	subroutine kpsInit(flag)
	
	implicit none

	character*24 idate
	double precision kpsRND
	double precision kpsinv,kpsinv2,dummy
	integer ir(250)
	integer ibm,jh,jm,js,k,i,flag
	integer kpsa,kpsb
	
	common /kps/ ir,kpsinv,kpsinv2,kpsa,kpsb

	kpsinv=2147483648.0**(-1)
	kpsinv2=kpsinv
	
	call fdate(idate)
	read (idate,'(11X,I2)') jh
	read (idate,'(14X,I2)') jm
	read (idate,'(17X,I2)') js
	ibm=65539*js+16807*jm+(11038+10437)*jh
	if (mod(ibm,2).eq.0) ibm=ibm+1
	ibm=mod(ibm,100000)

	if (flag.eq.1) ibm=65539

	do k=1,250
	  ibm=ibm*16807
	  ir(k)=0
	enddo

	do k=1,250
	  do i=1,32
	    ibm=ibm*16807
	    if (ibm.lt.0) ir(k)=ir(k)+1
	    ir(k)=ir(k)*2
	  enddo
	enddo
	
	kpsa=1
	kpsb=148
	
	do i=1,5000
	  dummy=kpsRND()
	enddo

	return
	end
*-------------------------------------------------------------------------------
	subroutine kpsSetMax(max)

	implicit none

	double precision kpsinv,kpsinv2
	integer ir(250)
	integer kpsa,kpsb,max
	
	common /kps/ ir,kpsinv,kpsinv2,kpsa,kpsb

	kpsinv2=kpsinv*max

	return
	end
*-------------------------------------------------------------------------------
	double precision function kpsRND()

	implicit none

	double precision kpsinv,kpsinv2
	integer kpsa,kpsb
	integer ir(250)
	
	common /kps/ ir,kpsinv,kpsinv2,kpsa,kpsb
	
	if (kpsa.eq.251) kpsa=1
	if (kpsb.eq.251) kpsb=1
	ir(kpsa)=ieor(ir(kpsa),ir(kpsb))
	kpsRND=kpsinv*iand(ir(kpsa),2147483647)
	kpsa=kpsa+1
	kpsb=kpsb+1

	return
	end
*-------------------------------------------------------------------------------
	integer function kpsINT()

	implicit none

	double precision kpsinv,kpsinv2
	integer kpsa,kpsb
	integer ir(250)
	
	common /kps/ ir,kpsinv,kpsinv2,kpsa,kpsb
	
	if (kpsa.eq.251) kpsa=1
	if (kpsb.eq.251) kpsb=1
	ir(kpsa)=ieor(ir(kpsa),ir(kpsb))
	kpsINT=iand(ir(kpsa),2147483647)
	kpsa=kpsa+1
	kpsb=kpsb+1

	return
	end
*-------------------------------------------------------------------------------
	integer function kpsINTf()

	implicit none

	double precision kpsinv,kpsinv2
	integer kpsa,kpsb
	integer ir(250)

	common /kps/ ir,kpsinv,kpsinv2,kpsa,kpsb

	if (kpsa.eq.251) kpsa=1
	if (kpsb.eq.251) kpsb=1
	ir (kpsa)=ieor(ir(kpsa),ir(kpsb))
	kpsINTf=ir(kpsa)
	kpsa=kpsa+1
	kpsb=kpsb+1

	return
	end
*-------------------------------------------------------------------------------
	integer function kpsMAX(maximum)
	
	implicit none

	double precision kpsinv,kpsinv2
	integer kpsa,kpsb,maximum
	integer ir(250)

	common /kps/ ir,kpsinv,kpsinv2,kpsa,kpsb
	
	if (kpsa.eq.251) kpsa=1
	if (kpsb.eq.251) kpsb=1
	ir(kpsa)=ieor(ir(kpsa),ir(kpsb))
	kpsMAX=int(kpsinv*maximum*iand(ir(kpsa),2147483647))
	kpsa=kpsa+1
	kpsb=kpsb+1

	return
	end
*-------------------------------------------------------------------------------
	integer function kpsMAX2()
	
	implicit none

	double precision kpsinv,kpsinv2
	integer kpsa,kpsb
	integer ir(250)

	common /kps/ ir,kpsinv,kpsinv2,kpsa,kpsb
	
	if (kpsa.eq.251) kpsa=1
	if (kpsb.eq.251) kpsb=1
	ir(kpsa)=ieor(ir(kpsa),ir(kpsb))
	kpsMAX2=int(kpsinv2*iand(ir(kpsa),2147483647))
	kpsa=kpsa+1
	kpsb=kpsb+1

	return
	end
*-------------------------------------------------------------------------------
