*===============================================================================
*	CLUSTER STATISTICS
*===============================================================================
* file: cvtpp.inc
* date: 19-jan-2001
*-------------------------------------------------------------------------------

	integer cl(mxa),clf(mxa),cs(mxa),csf(0:mxa),cc(mxa),ccf(mxa),csi(0:mxa)
	integer csif(0:mxa),cnum,cnumf,csover,clwait,clbrkt
	integer cgs(mxa),cgsa(mxa),cfpt(mxa),cfptf(mxa),cgsf(0:mxa)
	integer cgsfa(0:mxa)
	double precision cek(mxa),ceks(0:mxa)
	double precision cekf(0:mxa), ceksf(0:mxa)
	integer ccn,cen,ctrat,cstep,cbiggest,cbiggestf,clookfor,clookforf
	double precision mcts, omcts
	character*1  clnuctrl
	common /CVTPP/ cek,cekf,ceks,ceksf,mcts,omcts,cl,clf,cs,csf,cc,ccf,csi,
     $		       csif,cgs,cgsf,cgsa,cgsfa,cnum,cnumf,cfpt,cfptf,ccn,cen,
     $		       ctrat,cstep,csover,cbiggest,cbiggestf,clookfor,clookforf,
     $		       clwait,clbrkt,clnuctrl
*-------------------------------------------------------------------------------
* 	Common Block CVTPP (cvtpp.inc)
* cl()		: label of atom (cluster number, 0: unassigned)
* cs()		: size of cluster N
* cc()		: center-checker (0: not yet center, 1: was center)
* csi()		: size distribution of clusters
* cek()		: kinetic energy of cluster N
* ceks()	: kinetic energy summed over all clusters of one size
* cnum		: number of cluster/last cluster
* clf()		: label of atom (cluster number, 0: unassigned) for frenkel
* csf()		: size of cluster N for frenkel
* ccf()		: center-checker (0: not yet center, 1: was center) for frenkel
* csif()	: size distribution of frenkel clusters 
* cekf()	: kinetic energy of frenkel cluster N
* ceksf()	: kinetic energy summed over all frenkel clusters of one size
* cnumf		: number of cluster/last cluster for frenkel
* cgs		: number of carrier gas atoms in cluster
* cgsa		: number of carrier gas atoms in all clusters of size i
* cgsf		: number of carrier gas atoms in frenkel cluster
* cgsfa		: number of carrier gas atoms in all frenkel clusters of size i
* ccn		: number of atoms condensing on the cluster during ctrat
* cen		: number of atoms evaporating from the cluster during ctrat
* mcts		: temperature of biggest cluster (averaged)
* omcts		: temperature of old biggest cluster (averaged)
* ctrat		: time span for evap/cond-Statistics
* cstep		: time step of the last stoddard-call
* csover	: minimum cluster size before abort the simulation
* clwait	: wait time before abort when cluster is found
* clbrkt	: first time step a cluster>csover is found
* clnuctrl	: check for cluster in system to abort or not?
*		  X: yes, - = no
* cfpt()	: first passage time of given cluster size (Stillinger)
* cfptf()	: first passage time of given cluster size (tW/F)
*-------------------------------------------------------------------------------
