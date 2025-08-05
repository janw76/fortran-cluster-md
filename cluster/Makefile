#-------------------------------------------------------------------------------
# MakeFile zu Cluster und Hilfsprogrammen
#-------------------------------------------------------------------------------
# file: Makefile
# date: 03-feb-2002
#-------------------------------------------------------------------------------

CLUOBJ= main.o \
	init-cube1.o kps.o scalev.o stoddard.o\
	neighbor.o resize.o stopdead.o rates.o chg-natom.o\
	for-a-lj-n.o for-w-i9.o \
	int-verlet-dmm.o \
	read-3d.o read-rst.o \
	wrt-acc.o wrt-ave.o wrt-csi.o wrt-rst.o wrt-vel.o wrt-xyz.o wrt-fpt.o \
	chk-nucleation.o


PBCOBJ=	main.o \
	init-cube1.o kps.o scalev.o stoddard-pbc-n.o chg-natom.o\
	neighbor-pbc.o resize.o  stopdead-m.o rates.o\
	for-a-lj-n-pbc.o for-w-dummy.o \
	int-verlet-dmm.o \
	read-3d.o read-rst.o \
	wrt-acc.o wrt-ave.o wrt-csi.o wrt-rst.o wrt-vel.o wrt-xyz.o \
	chk-nucleation.o

MAVOBJ=	main.o \
	init-cube1.o kps.o scalev.o stoddard-pbc-n.o \
	neighbor-pbc-mav.o resize.o stopdead.o rates.o chg-natom.o\
	for-a-lj-n-pbc.o for-w-dummy.o \
	int-verlet-mav.o \
	read-3d.o read-rst.o \
	wrt-acc.o wrt-ave.o wrt-csi.o wrt-rst.o wrt-vel.o wrt-xyz.o \
	chk-nucleation.o

NUC=		main.o \
        init-cube1.o kps.o scalev-m.o stoddard-pbc-n-m.o \
        neighbor-pbc-mav.o resize.o stopdead-m.o rates.o chg-natom.o\
        for-a-lj-n-pbc-m.o for-w-dummy.o \
        int-verlet-mav.o \
        read-3d.o read-rst.o \
        wrt-acc.o wrt-ave.o wrt-csi.o wrt-rst.o wrt-vel.o wrt-xyz-m-vmd.o \
	chk-nucleation.o

SILENT=		main.o \
        init-cube1.o kps.o scalev-m.o stoddard-pbc-n-m.o \
        neighbor-pbc-mav-silent.o resize.o stopdead-m.o rates.o chg-natom.o\
        for-a-lj-n-pbc-m.o for-w-dummy.o \
        int-verlet-mav.o \
        read-3d.o read-rst.o \
        wrt-acc.o wrt-ave.o wrt-csi.o wrt-rst.o wrt-vel.o wrt-xyz-m-vmd.o wrt-fpt.o \
	chk-nucleation.o

FRENKEL=		main.o \
        init-cube1.o kps.o scalev-m.o stoddard-pbc-n-m-fre.o \
        neighbor-pbc-mav-silent.o resize.o stopdead-m.o rates.o chg-natom.o\
        for-a-lj-n-pbc-m.o for-w-dummy.o \
        int-verlet-mav.o \
        read-3d.o read-rst.o \
        wrt-acc.o wrt-ave.o wrt-csi-fre.o wrt-rst.o wrt-vel.o wrt-xyz-m-vmd.o wrt-fpt.o\
	chk-nucleation.o

DIFFUSION=		main.o \
        init-center.o kps.o scalev-m.o stoddard-pbc-n-m-fre.o \
        neighbor-pbc-mav-silent.o resize.o stopdead-m.o rates.o chg-natom.o\
        for-a-lj-n-pbc-m.o for-w-dummy.o \
        int-verlet-mav.o \
        read-3d.o read-rst.o \
        wrt-acc.o wrt-ave.o wrt-csi-fre.o wrt-rst.o wrt-vel.o wrt-xyz-m-vmd.o \
	chk-nucleation.o    
    
MOVIEFRENKEL=		main.o \
        init-cube1.o kps.o scalev-m.o stoddard-pbc-n-m-fre-mov.o \
        neighbor-pbc-mav-silent.o resize.o stopdead-m.o rates.o chg-natom.o\
        for-a-lj-n-pbc-m.o for-w-dummy.o \
        int-verlet-mav.o \
        read-3d.o read-rst.o \
        wrt-acc.o wrt-ave.o wrt-csi-fre.o wrt-rst.o wrt-vel.o wrt-xyz-m-vmd.o \
	chk-nucleation.o

NOSE=		main.o \
        init-cube1.o kps.o scalev-nh-m.o stoddard-pbc-n-m.o \
        neighbor-pbc-mav-silent.o resize.o stopdead-m.o rates.o chg-natom.o\
        for-a-lj-n-pbc-m.o for-w-dummy.o \
        int-verlet-nh-mav.o \
        read-3d.o read-rst.o \
        wrt-acc.o wrt-ave.o wrt-csi.o wrt-rst.o wrt-vel.o wrt-xyz-m-vmd.o \
	chk-nucleation.o

MULOBJ= main.o \
        init-cube1.o kps.o scalev-m.o stoddard-pbc-n-m.o \
        neighbor-pbc-mav.o resize.o stopdead-m.o rates.o chg-natom.o\
        for-a-lj-n-pbc-m.o for-w-dummy.o \
        int-verlet-mav.o \
        read-3d.o read-rst.o \
        wrt-acc.o wrt-ave.o wrt-csi.o wrt-rst.o wrt-vel.o wrt-xyz-m.o \
	chk-nucleation.o

TSTOBJ= main.o \
        init-cube1.o kps.o scalev-m.o stoddard-pbc-n-m.o \
        neighbor-pbc-mav.o resize.o stopdead-m.o rates.o chg-natom.o\
        for-a-lj-n-pbc-m.o for-w-dummy.o \
        int-verlet-mav.o \
        read-3d.o read-rst.o \
        wrt-acc.o wrt-ave.o wrt-csi.o wrt-rst.o wrt-vel-m-ats1.o wrt-xyz-m-ats1.o

RSTOBJ=	main.o \
        init-cube1.o kps.o scalev-m.o stoddard-pbc-n-m.o \
        neighbor-pbc-mav-silent.o resize.o stopdead-m.o rates.o chg-natom.o\
        for-a-lj-n-pbc-m.o for-w-dummy.o \
        int-verlet-mav.o \
        read-3d.o read-rst.o \
        wrt-acc.o wrt-ave.o wrt-csi.o wrt-rst-multi.o wrt-vel.o wrt-xyz-m.o \
	chk-nucleation.o

THINOUTOBJ=	thinout.o wrt-xyz.o read-3d.o read-xyz.o

SELATS_OBJ=	selats.o wrt-xyz-m1.o read-3d.o read-xyz-m.o

SELCLU_OBJ=	selclu.o stoddard.o wrt-xyz1.o read-3d.o read-xyz.o

SELCLUPOBJ=	selclu.o stoddard-pbc.o wrt-xyz3.o read-3d.o read-xyz.o

SIZER__OBJ=	sizer.o read-3d.o

MESHER_OBJ=	mesher.o read-3d.o

VELTESTOBJ=	veltest.o read-3d.o

PRESSW_OBJ=	pressw.o read-3d.o for-w-i9.o wrt-ave.o read-xyz.o

MARKER_OBJ=	marker.o wrt-xyz-m.o read-3d.o read-3d.o

PRESS__OBJ=	press.o read-3d.o

DENSITYOBJ=	density.o stoddard.o read-3d.o read-xyz.o

DENSITPOBJ=	densitp.o stoddard-pbc.o read-3d.o read-xyz.o

SYSDENSOBJ=	sysdens.o read-3d.o read-xyz.o

COLGAS_OBJ=	colgas.o stoddard.o wrt-xyz-m.o read-3d.o read-xyz.o

COLGASPOBJ=	colgas.o stoddard-pbc.o wrt-xyz-m.o read-3d.o read-xyz.o

RDF____OBJ=	rdf.o read-3d.o read-xyz.o

TIMEFILOBJ=	timefile.o read-3d.o

MAKELOGOBJ=	makelog.o read-3d.o

BARS___OBJ= 	bars.o read-3d.o

REVERSEOBJ=	reverse.o read-3d.o read-xyz.o wrt-xyz.o

RATEXYZOBJ=	rate-xyz.o read-3d.o read-xyz.o stoddard-pbc.o rates.o

STOTESTOBJ=	stotest.o read-3d.o read-xyz.o stoddard-pbc-n.o \
		wrt-csi.o neighbor-pbc.o

RATSTATOBJ=	ratstat.o

RATSTTTOBJ=	ratstat-t.o

BIGCL__OBJ=	bigcl.o stoddard.o wrt-xyz-m.o read-3d.o read-xyz.o

BIGCLP_OBJ=	bigcl.o stoddard-pbc.o wrt-xyz-m.o read-3d.o read-xyz.o

FOCUSCLP_O=	focusclp.o stoddard-pbc.o read-3d.o read-xyz.o wrt-xyz.o

PASSAGE=	passage-time.o

MPT_OBJ=	mpt.o

TEMPCL=		tempavg.o

TEMPCL2=	tempavg2.o

MTEMP=		meantemp.o

THIST=		t-hist.o

DIFFUS=		diffusivity.o

PROBAB=	probability-distrib.o


# SGI compiler flags (IRIX 6.5)
# FLAGS= -O3 -fullwarn -woff 2274
# LINUX RedHat Options
# FLAGS= -O3
# CRAY compiler options ---SET FF=fort77!!!---
# FLAGS= -i32 -O3 -Xm -M 1110
# sunfire Solaris 8 flags
# FLAGS= -xtarget=ultra3 -xarch=v9b
# Linux / G77 / Intel
# FLAGS= -O3 -malign-double
# Intel Linux Fortran Compiler
FLAGS= -O3 -warn -traceback
# FF=	f77
FF=	ifort

cluster:	$(CLUOBJ)
		$(FF) $(FLAGS) $(CLUOBJ) -o cluster

pbc:		$(PBCOBJ)
		$(FF) $(FLAGS) $(PBCOBJ) -o cluster

mav:		$(MAVOBJ)
		$(FF) $(FLAGS) $(MAVOBJ) -o cluster

mul:		$(MULOBJ)
		$(FF) $(FLAGS) $(MULOBJ) -o cluster-m

vmd:		$(VMD)
		$(FF) $(FLAGS) $(VMD) -o clu-vmd

nuc:		$(NUC)
		$(FF) $(FLAGS) $(NUC) -o clu-nuc

silent:		$(SILENT)
		$(FF) $(FLAGS) $(SILENT) -o clu-silent

frenkel:	$(FRENKEL)
		$(FF) $(FLAGS) $(FRENKEL) -o clu-fre
        
diffusion:	$(DIFFUSION)
		$(FF) $(FLAGS) $(DIFFUSION) -o diffusion
        
moviefrenkel:	$(MOVIEFRENKEL)
		$(FF) $(FLAGS) $(MOVIEFRENKEL) -o clu-fre-mov

nose:		$(NOSE)
		$(FF) $(FLAGS) $(NOSE) -o clu-nose

rst:		$(RSTOBJ)
		$(FF) $(FLAGS) $(RSTOBJ) -o clu-rst

test:		$(TSTOBJ)
		$(FF) $(FLAGS) $(TSTOBJ) -o clutest

thinyout:	bdir $(THINOUTOBJ)
		$(FF) $(FLAGS) $(THINOUTOBJ) -o $(HOME)/bin/thinout

selats:		bdir $(SELATS_OBJ)
		$(FF) $(FLAGS) $(SELATS_OBJ) -o $(HOME)/bin/selats

selclu:		bdir $(SELCLU_OBJ)
		$(FF) $(FLAGS) $(SELCLU_OBJ) -o $(HOME)/bin/selclu

selclup:	bdir $(SELCLUPOBJ)
		$(FF) $(FLAGS) $(SELCLUPOBJ) -o $(HOME)/bin/selclu

sizer:		bdir $(SIZER__OBJ)
		$(FF) $(FLAGS) $(SIZER__OBJ) -o $(HOME)/bin/sizer

mesher:		bdir $(MESHER_OBJ)
		$(FF) $(FLAGS) $(MESHER_OBJ) -o $(HOME)/bin/mesher

veltest:	bdir $(VELTESTOBJ)
		$(FF) $(FLAGS) $(VELTESTOBJ) -o $(HOME)/bin/veltest

pressw:		bdir $(PRESSW_OBJ)
		$(FF) $(FLAGS) $(PRESSW_OBJ) -o $(HOME)/bin/pressw

marker:		bdir $(MARKER_OBJ)
		$(FF) $(FLAGS) $(MARKER_OBJ) -o $(HOME)/bin/marker

press:		bdir $(PRESS__OBJ)
		$(FF) $(FLAGS) $(PRESS__OBJ) -o $(HOME)/bin/press

density:	bdir $(DENSITYOBJ)
		$(FF) $(FLAGS) $(DENSITYOBJ) -o $(HOME)/bin/density

densityp:	bdir $(DENSITPOBJ)
		$(FF) $(FLAGS) $(DENSITPOBJ) -o $(HOME)/bin/density

sysdens:	bdir $(SYSDENSOBJ)
		$(FF) $(FLAGS) $(SYSDENSOBJ) -o $(HOME)/bin/sysdens

colgas:		bdir $(COLGAS_OBJ)
		$(FF) $(FLAGS) $(COLGAS_OBJ) -o $(HOME)/bin/colgas

colgasp:	bdir $(COLGASPOBJ)
		$(FF) $(FLAGS) $(COLGASPOBJ) -o $(HOME)/bin/colgas

rdf:		bdir $(RDF____OBJ)
		$(FF) $(FLAGS) $(RDF____OBJ) -o $(HOME)/bin/rdf

timefile:	bdir $(TIMEFILOBJ)
		$(FF) $(FLAGS) $(TIMEFILOBJ) -o $(HOME)/bin/timefile

makelog:	bdir $(MAKELOGOBJ)
		$(FF) $(FLAGS) $(MAKELOGOBJ) -o $(HOME)/bin/makelog

bars:		bdir $(BARS___OBJ)
		$(FF) $(FLAGS) $(BARS___OBJ) -o $(HOME)/bin/bars

reverse:	bdir $(REVERSEOBJ)
		$(FF) $(FLAGS) $(REVERSEOBJ) -o $(HOME)/bin/reverse

rate-xyz:	bdir $(RATEXYZOBJ)
		$(FF) $(FLAGS) $(RATEXYZOBJ) -o $(HOME)/bin/rate-xyz

stotest:	bdir $(STOTESTOBJ)
		$(FF) $(FLAGS) $(STOTESTOBJ) -o $(HOME)/bin/stotest

ratstat:	bdir $(RATSTATOBJ)
		$(FF) $(FLAGS) $(RATSTATOBJ) -o $(HOME)/bin/ratstat

ratstatt:	bdir $(RATSTTTOBJ)
		$(FF) $(FLAGS) $(RATSTTTOBJ) -o $(HOME)/bin/ratstat-t

bigcl:		bdir $(BIGCL__OBJ)
		$(FF) $(FLAGS) $(BIGCL__OBJ) -o $(HOME)/bin/bigcl

bigclp:		bdir $(BIGCLP_OBJ)
		$(FF) $(FLAGS) $(BIGCL__OBJ) -o $(HOME)/bin/bigcl

focusclp:	bdir $(FOCUSCLP_O)
		$(FF) $(FLAGS) $(FOCUSCLP_O) -o $(HOME)/bin/focusclp

passagetime:	bdir $(PASSAGE)
		$(FF) $(FLAGS) $(PASSAGE) -o $(HOME)/bin/passage-time

mpt:		bdir $(MPT_OBJ)
		$(FF) $(FLAGS) $(MPT_OBJ) -o $(HOME)/bin/mpt

tempavg:	bdir $(TEMPCL)
		$(FF) $(FLAGS) $(TEMPCL) -o $(HOME)/bin/tempavg

meantemp:	bdir $(MTEMP)
		$(FF) $(FLAGS) $(MTEMP) -o $(HOME)/bin/meantemp

t-hist:		bdir $(THIST)
		$(FF) $(FLAGS) $(THIST) -o $(HOME)/bin/t-hist

tempavg2:	bdir $(TEMPCL2)
		$(FF) $(FLAGS) $(TEMPCL2) -o $(HOME)/bin/tempavg2

diffusivity:	bdir $(DIFFUS)
		$(FF) $(FLAGS) $(DIFFUS) -o $(HOME)/bin/diffusivity
		
probability:	bdir $(PROBAB)
		$(FF) $(FLAGS) $(PROBAB) -o $(HOME)/bin/probability-distrib

clean:
	touch remove.o
	rm *.o

archive:
	touch archive.tar.gz
	rm archive.tar.gz
	tar -cvf archive.tar *.f *.inc Makefile sample.3d README \
                             helpers/*
	gzip archive.tar

bdir:	
	@if test ! -d $(HOME)/bin ; then \
	  echo 'Directory '$(HOME)/bin' missing ... creating it!' ; \
	  mkdir $(HOME)/bin ; \
	  echo "Don't forget to add it to your PATH-environment variable" ; \
	  sleep 10 ; \
	fi

.f.o:	
	$(FF) $(FLAGS) -c $<


#------------------------------------------------------------------------------
# Dependencies for main program modules
#------------------------------------------------------------------------------

main.o:		main.f \
		const.inc atomc.inc filep.inc timep.inc energ.inc boxpp.inc \
		atomp.inc

init-cube1.o:	init-cube1.f \
		const.inc atomp.inc atomc.inc boxpp.inc energ.inc
        
init-center.o:	init-center.f \
		const.inc atomp.inc atomc.inc boxpp.inc energ.inc

init-2cubes.o:	init-2cubes.f \
		const.inc atomp.inc atomc.inc boxpp.inc energ.inc

init-2atom.o:	init-2atom.f \
		const.inc atomp.inc atomc.inc boxpp.inc energ.inc

init-ball1.o:	init-cube1.f \
		const.inc atomp.inc atomc.inc boxpp.inc energ.inc

kps.o:		kps.f

scalev.o:	scalev.f \
		const.inc atomp.inc atomc.inc energ.inc

scalev-m.o:	scalev-m.f \
		const.inc atomp.inc atomc.inc energ.inc timep.inc

scalev-nh-m.o:	scalev-nh-m.f \
		const.inc atomp.inc atomc.inc energ.inc timep.inc

scalev-clu.o:	scalev-clu.f \
		const.inc atomp.inc atomc.inc energ.inc cvtpp.inc

resize.o:	resize.f \
		const.inc atomc.inc atomp.inc boxpp.inc

resize2.o:	resize2.f \
		const.inc atomc.inc atomp.inc boxpp.inc

stopdead.o:	stopdead.f \
		const.inc atomc.inc atomp.inc energ.inc

rates.o:	rates.f \
		const.inc atomp.inc atomc.inc cvtpp.inc

chg-natom.o:	chg-natom.f \
		const.inc atomp.inc atomc.inc cvtpp.inc

stoddard.o:	stoddard.f \
		const.inc atomp.inc atomc.inc cvtpp.inc

stoddard-pbc.o:	stoddard-pbc.f \
		const.inc atomp.inc atomc.inc cvtpp.inc boxpp.inc

stoddard-pbc-n.o:	stoddard-pbc-n.f \
		const.inc atomp.inc atomc.inc cvtpp.inc boxpp.inc

stoddard-pbcsw-n.o:	stoddard-pbcsw-n.f \
		const.inc atomp.inc atomc.inc cvtpp.inc boxpp.inc

stoddard-pbc-n-m.o:     stoddard-pbc-n-m.f \
                const.inc atomp.inc atomc.inc cvtpp.inc boxpp.inc

stoddard-pbc-n-m-fre-mov.o:     stoddard-pbc-n-m-fre-mov.f \
                const.inc atomp.inc atomc.inc cvtpp.inc boxpp.inc

stoddard-dummy.o:	stoddard-dummy.f

neighbor-dummy.o:	neighbor-dummy.f \
		const.inc dmmpp.inc

neighbor.o:	neighbor.f \
		const.inc atomp.inc atomc.inc dmmpp.inc

neighbor-pbc.o:	neighbor-pbc.f \
		const.inc atomp.inc atomc.inc boxpp.inc dmmpp.inc

neighbor-pbc-mav.o:	neighbor-pbc-mav.f \
		const.inc atomp.inc atomc.inc boxpp.inc dmmpp.inc

neighbor-pbc-mav-silent.o:	neighbor-pbc-mav-silent.f \
		const.inc atomp.inc atomc.inc boxpp.inc dmmpp.inc

neighbor-pbcsw-mav.o:	neighbor-pbcsw-mav.f \
		const.inc atomp.inc atomc.inc boxpp.inc dmmpp.inc

for-a-dummy.o:	for-a-dummy.f

for-a-lj1.o:	for-a-lj1.f \
		const.inc atomp.inc atomc.inc energ.inc

for-a-lj2.o:	for-a-lj2.f \
		const.inc atomp.inc atomc.inc energ.inc

for-a-lj3.o:	for-a-lj3.f \
		const.inc atomp.inc atomc.inc energ.inc

for-a-lj4.o:	for-a-lj4.f \
		const.inc atomp.inc atomc.inc energ.inc

for-a-lj-n.o:	for-a-lj-n.f \
		const.inc atomp.inc atomc.inc energ.inc

for-a-lj-n-pbc.o:	for-a-lj-n-pbc.f \
		const.inc atomp.inc atomc.inc energ.inc boxpp.inc

for-a-lj-n-pbc-m.o:	for-a-lj-n-pbc-m.f \
		const.inc atomp.inc atomc.inc energ.inc boxpp.inc

for-a-lj-n-pbcsw.o:	for-a-lj-n-pbcsw.f \
		const.inc atomp.inc atomc.inc energ.inc boxpp.inc

for-w-dummy.o:	for-w-dummy.f

for-w-i12.o:	for-w-i12.f \
		const.inc atomp.inc atomc.inc boxpp.inc energ.inc

for-w-i9.o:	for-w-i9.f \
		const.inc atomp.inc atomc.inc boxpp.inc energ.inc

for-w-i9g.o:	for-w-i9g.f \
		const.inc atomp.inc atomc.inc boxpp.inc energ.inc

for-w-iX.o:	for-w-iX.f \
		const.inc atomp.inc atomc.inc boxpp.inc energ.inc

for-w-i9-dbl.o:	for-w-i9-dbl.f \
		const.inc atomp.inc atomc.inc boxpp.inc energ.inc

for-w-lj.o:	for-w-lj.f \
		const.inc atomp.inc atomc.inc boxpp.inc energ.inc

int-euler1.o:	int-euler1.f \
		const.inc atomp.inc atomc.inc timep.inc energ.inc

int-euler2.o:	int-euler2.f \
		const.inc atomp.inc atomc.inc timep.inc energ.inc

int-gear3.o:	int-gear3.f \
		const.inc atomp.inc atomc.inc timep.inc energ.inc

int-gear4.o:	int-gear4.f \
		const.inc atomp.inc atomc.inc timep.inc energ.inc

int-gear5.o:	int-gear5.f \
		const.inc atomp.inc atomc.inc timep.inc energ.inc

int-verlet.o:	int-verlet.f \
		const.inc atomp.inc atomc.inc timep.inc energ.inc

int-verlet-dmm.o:	int-verlet-dmm.f \
		const.inc atomp.inc atomc.inc timep.inc energ.inc dmmpp.inc

int-verlet-mav.o:	int-verlet-mav.f \
		const.inc atomp.inc atomc.inc timep.inc energ.inc dmmpp.inc

int-verlet-nh-mav.o:	int-verlet-nh-mav.f \
		const.inc atomp.inc atomc.inc timep.inc energ.inc dmmpp.inc

read-3d.o:	read-3d.f \
		const.inc filep.inc atomp.inc timep.inc energ.inc boxpp.inc \
		dmmpp.inc cvtpp.inc

read-rst.o:	read-rst.f \
		const.inc atomp.inc atomc.inc filep.inc energ.inc

read-xyz.o:	read-xyz.f \
		const.inc atomp.inc atomc.inc boxpp.inc timep.inc

read-xyz-m.o:	read-xyz-m.f \
		const.inc atomp.inc atomc.inc boxpp.inc timep.inc


read-vel.o:	read-vel.f \
		const.inc atomp.inc atomc.inc 

wrt-acc.o:	wrt-acc.f \
		const.inc atomp.inc atomc.inc

wrt-ave.o:	wrt-ave.f \
		const.inc atomp.inc boxpp.inc filep.inc energ.inc

wrt-csi.o:	wrt-csi.f \
		const.inc atomp.inc atomc.inc cvtpp.inc

wrt-csi-fre.o:	wrt-csi-fre.f \
		const.inc atomp.inc atomc.inc cvtpp.inc

wrt-rst.o:	wrt-rst.f \
		const.inc atomp.inc energ.inc filep.inc

wrt-rst-multi.o:	wrt-rst-multi.f \
		const.inc atomp.inc energ.inc filep.inc

wrt-vel.o:	wrt-vel.f \
		const.inc atomp.inc atomc.inc 

wrt-vel-m-ats1.o:	wrt-vel-m-ats1.f \
		const.inc atomp.inc atomc.inc 

wrt-xyz.o:	wrt-xyz.f \
		const.inc atomp.inc atomc.inc boxpp.inc timep.inc

wrt-xyz-m.o:	wrt-xyz-m.f \
		const.inc atomp.inc atomc.inc boxpp.inc timep.inc

wrt-xyz-m-ats1.o:	wrt-xyz-m-ats1.f \
		const.inc atomp.inc atomc.inc boxpp.inc timep.inc

wrt-xyz-m1.o:	wrt-xyz-m1.f \
		const.inc atomp.inc atomc.inc boxpp.inc timep.inc

wrt-xyz1.o:	wrt-xyz1.f \
		const.inc atomp.inc atomc.inc boxpp.inc cvtpp.inc

wrt-xyz2.o:	wrt-xyz2.f \
		const.inc atomp.inc atomc.inc boxpp.inc cvtpp.inc

wrt-xyz3.o:	wrt-xyz3.f \
		const.inc atomp.inc atomc.inc boxpp.inc cvtpp.inc
		
wrt-fpt.o:	wrt-fpt.f \
		const.inc atomp.inc atomc.inc cvtpp.inc timep.inc

chk-nucleation.o:	chk-nucleation.f \
		const.inc atomp.inc atomc.inc cvtpp.inc timep.inc

#------------------------------------------------------------------------------
# Dependencies for helper programs
#------------------------------------------------------------------------------

thinout.o:	helpers/thinout.f \
		const.inc atomp.inc atomc.inc timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/thinout.f

selats.o:	helpers/selats.f \
		const.inc atomp.inc atomc.inc timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/selats.f

selclu.o:	helpers/selclu.f \
		const.inc atomp.inc atomc.inc timep.inc filep.inc cvtpp.inc
		$(FF) $(FLAGS) -c helpers/selclu.f

sizer.o: 	helpers/sizer.f \
		const.inc atomp.inc atomc.inc timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/sizer.f

mesher.o:	helpers/mesher.f \
		const.inc atomp.inc atomc.inc timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/mesher.f

veltest.o: 	helpers/veltest.f \
		timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/veltest.f

pressw.o: 	helpers/pressw.f \
		const.inc atomp.inc atomc.inc timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/pressw.f

marker.o:	helpers/marker.f \
		const.inc atomp.inc atomc.inc timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/marker.f

press.o: 	helpers/press.f \
		timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/press.f

density.o: 	helpers/density.f \
		const.inc atomp.inc atomc.inc boxpp.inc timep.inc filep.inc \
		cvtpp.inc
		$(FF) $(FLAGS) -c helpers/density.f

densitp.o: 	helpers/densitp.f \
		const.inc atomp.inc atomc.inc boxpp.inc timep.inc filep.inc \
		cvtpp.inc
		$(FF) $(FLAGS) -c helpers/densitp.f

sysdens.o: 	helpers/sysdens.f \
		const.inc atomp.inc atomc.inc boxpp.inc timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/sysdens.f

colgas.o: 	helpers/colgas.f \
		const.inc atomp.inc timep.inc filep.inc cvtpp.inc
		$(FF) $(FLAGS) -c helpers/colgas.f

rdf.o:	 	helpers/rdf.f \
		const.inc atomp.inc atomc.inc timep.inc boxpp.inc filep.inc
		$(FF) $(FLAGS) -c helpers/rdf.f

timefile.o:	helpers/timefile.f \
		const.inc timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/timefile.f

makelog.o:	helpers/makelog.f \
		const.inc timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/makelog.f

bars.o:		helpers/bars.f \
		const.inc atomp.inc timep.inc filep.inc cvtpp.inc
		$(FF) $(FLAGS) -c helpers/bars.f

reverse.o:	helpers/reverse.f \
		const.inc atomp.inc timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/reverse.f

rate-xyz.o:	helpers/rate-xyz.f \
		const.inc atomp.inc timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/rate-xyz.f

stotest.o:	helpers/stotest.f \
		const.inc atomp.inc timep.inc filep.inc
		$(FF) $(FLAGS) -c helpers/stotest.f

ratstat.o:	helpers/ratstat.f 
		$(FF) $(FLAGS) -c helpers/ratstat.f

ratstat-t.o:	helpers/ratstat-t.f 
		$(FF) $(FLAGS) -c helpers/ratstat-t.f

bigcl.o: 	helpers/bigcl.f \
		const.inc atomp.inc timep.inc filep.inc cvtpp.inc
		$(FF) $(FLAGS) -c helpers/bigcl.f

focusclp.o: 	helpers/focusclp.f \
		const.inc atomp.inc atomc.inc boxpp.inc timep.inc filep.inc \
		cvtpp.inc
		$(FF) $(FLAGS) -c helpers/focusclp.f

meantemp.o:     helpers/meantemp.f
		$(FF) $(FLAGS) -c helpers/meantemp.f

passage-time.o:	helpers/passage-time.f
		$(FF) $(FLAGS) -c helpers/passage-time.f

mpt.o:		helpers/mpt.f
		$(FF) $(FLAGS) -c helpers/mpt.f

tempavg.o:	helpers/tempavg.f
		$(FF) $(FLAGS) -c helpers/tempavg.f

tempavg2.o:	helpers/tempavg2.f
		$(FF) $(FLAGS) -c helpers/tempavg2.f

t-hist.o:	helpers/t-hist.f
		$(FF) $(FLAGS) -c helpers/t-hist.f

diffusivity.o:	helpers/diffusivity.f
		$(FF) $(FLAGS) -c helpers/diffusivity.f

probability-distrib.o:	helpers/probability-distrib.f
		$(FF) $(FLAGS) -c helpers/probability-distrib.f