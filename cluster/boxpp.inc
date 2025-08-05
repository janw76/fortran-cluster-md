*===============================================================================
*	BOX PARAMETERS
*===============================================================================
* file: boxpp.inc
*-------------------------------------------------------------------------------

	double precision boxx,boxy,boxz
	double precision beng,bpot
	double precision boxxn,boxyn,boxzn
	double precision boxxs,boxys,boxzs
	character*1 bcor,bshr,bperx,bpery,bperz
	integer boxt
	common /BOXPP/ boxx,boxy,boxz,boxxn,boxyn,boxzn,
     $	               boxxs,boxys,boxzs,bpot,beng,boxt,
     $		       bcor,bshr,bperx,bpery,bperz

*-------------------------------------------------------------------------------
*	Common Block BOXPP (boxpp.inc)
* boxx		: Length of the box in X-direction
* boxy		: Length of the box in Y-direction
* boxz		: Length of the box in Z-direction
* beng		: Energy parameter for repulsive wall potential (in K)
* bpot		: Exponent for potential (unused)
* boxxn		: New length of box in X-direction
* boxyn		: New length of box in Y-direction
* boxzn		: New length of box in Z-direction
* boxxs		: maximum step length in X-direction for box resize
* boxys		: maximum step length in Y-direction for box resize
* boxzs		: maximum step length in Z-direction for box resize
* bcor		: show box corners as H-atoms in first step of xyz-file
* boxt		: Time steps between box resizes
*		  -: do not show box corners
*		  X: show box corners
*		  A: show box corners in every time step (not implemented yet)
* bshr		: resize control flag
*		  X: resize until new size is reached
*		  -: do not resize
* bperx		: Periodicity of box in X-direction, X: Yes, -: No
* bpery		: Periodicity of box in Y-direction, X: Yes, -: No
* bperz		: Periodicity of box in Z-direction, X: Yes, -: No
*-------------------------------------------------------------------------------
