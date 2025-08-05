*===============================================================================
*	CLUSTER STATISTICS
*===============================================================================
* file: cvtpp.inc
* date: 19-jan-2001
*-------------------------------------------------------------------------------

	integer cl(mxa),cs(mxa),cc(mxa),csi(0:mxa),cnum
	double precision cek(mxa),ceks(0:mxa)
	integer ccn,cen,ctrat,cstep
	double precision mcts,omcts
	common /CVTPP/ cek,ceks,cl,cs,cc,csi,cnum,
     $		       ccn,cen,ctrat,cstep,mcts,omcts

*-------------------------------------------------------------------------------
* 	Common Block CVTPP (cvtpp.inc)
* cl()		: label of atom (cluster number, 0: unassigned)
* cs()		: size of cluster N
* cc()		: center-checker (0: not yet center, 1: was center)
* csi()		: size distribution of clusters
* cek()		: kinetic energy of cluster N
* ceks()	: kinetic energy summed over all clusters of one size
* cnum		: number of cluster/last cluster
* ccn		: number of atoms condensing on the cluster during ctrat
* cen		: number of atoms evaporatin from the cluster during ctrat
* mcts		: temperature of biggest cluster (averaged)
* omcts		: temperature of old biggest cluster (averaged)
* ctrat		: time span for evap/cond-Statistics
* cstep		: time step of the last stoddard-call
*-------------------------------------------------------------------------------
