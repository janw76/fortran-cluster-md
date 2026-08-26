*===============================================================================
*	RANDOM NUMBER GENERATOR
*===============================================================================
* file: kps.f
*-------------------------------------------------------------------------------
	subroutine kpsInit(flag)

	implicit none

	character*24 idate
	character*32 cseed
	double precision kpsRND
	double precision kpsinv,kpsinv2,dummy
	integer ir(250)
	integer ibm,jh,jm,js,k,i,flag
	integer kpsa,kpsb
	integer seed,ios,envlen,envstat
	logical seeded

	common /kps/ ir,kpsinv,kpsinv2,kpsa,kpsb

	kpsinv=2147483648.0**(-1)
	kpsinv2=kpsinv

*** CLUSTER_SEED overrides both the wall-clock (flag=0) and fixed (flag=1)
*** seeding below, so every kpsInit caller (main.f and the init-*.f family,
*** which reseed independently based on rndctrl) picks up the same explicit
*** seed. An unset/empty variable falls through silently; a set-but-invalid
*** value warns once and falls through to the default seeding for this call.
	seeded=.false.
	cseed=' '
	call get_environment_variable('CLUSTER_SEED',cseed,envlen,envstat)
	if (envstat.eq.0 .and. cseed.ne.' ') then
	  read (cseed,*,iostat=ios) seed
	  if (ios.eq.0) then
*** RANDU-style multiplicative generator needs a positive, odd seed.
	    ibm=abs(seed)
	    ibm=ior(ibm,1)
	    seeded=.true.
	    write (*,'(A,I0,A)') 'RNG seed: ',ibm,' (CLUSTER_SEED)'
	  else
	    write (*,'(A,A,A)') 'WARNING: CLUSTER_SEED=''',trim(cseed),
     $        ''' is not a valid integer; using default seeding'
	  endif
	else if (envstat.eq.-1) then
	  write (*,'(A)') 'WARNING: CLUSTER_SEED value too long; '//
     $      'using default seeding'
	endif

	if (.not.seeded) then
	  call fdate(idate)
	  read (idate,'(11X,I2)') jh
	  read (idate,'(14X,I2)') jm
	  read (idate,'(17X,I2)') js
	  ibm=65539*js+16807*jm+(11038+10437)*jh
	  if (mod(ibm,2).eq.0) ibm=ibm+1
	  ibm=mod(ibm,100000)

	  if (flag.eq.1) ibm=65539

	  if (flag.eq.1) then
	    write (*,'(A,I0,A)') 'RNG seed: ',ibm,' (fixed)'
	  else
	    write (*,'(A,I0,A)') 'RNG seed: ',ibm,' (wall-clock)'
	  endif
	endif

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
