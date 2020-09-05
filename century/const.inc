c*****************************************************************
c	$RCSfile: const.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: const.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2004/08/30 16:52:23  conklind
cImplement /save_vars/.
c
cRevision 1.1  2003/11/20 21:44:37  conklind
cInitial revision
c
c******************************************************************

c This source file contains constants, distinguishable as variable
c names in all caps, and also contains a labeled common /save_vars/
c All the persistent local variables in the original Century source
c files have been gathered up into the /save_vars/ common so that
c they can be reinitialized to zero between cells in gridded runs
c by the subroutine cen_zero().  "persistent" local variables are
c those declared with the keyword "save", meaning that their values
c persist from one invocation of the function to the next - like 
c local "static" variables in C.
c /save_vars/ added by Conklin  8/30/04

c******************************************************************
c Constant definitions begin here.

	integer CRPSYS, FORSYS, SAVSYS,
     +          MAXIEL, N, P, S,
     +          MAXLYR, MONTHS, SRFC, SOIL,
     +          CPARTS, ABOVE, BELOW,
     +          FPARTS, LEAF, FROOT, FBRCH, LWOOD, CROOT,
     +          NEWFOR, OLDFOR,
     +          IMIN, IMAX, IVAL,
     +          INTCPT, SLOPE,
     +          UNLABL, LABELD, ISOS,
     +          ESTOR, ESOIL, ENFIX, EFERT
	real PEEDEE, PI, PI2

c****   CRPSYS, FORSYS, SAVSYS is an enumeration for the system
	parameter (CRPSYS = 1)
	parameter (FORSYS = 2)
	parameter (SAVSYS = 3)

c****   MAXIEL is the total # of elements available
c       N = Nitrogen, P = Phosphorus, S = Sulphur
	parameter (MAXIEL = 3)
	parameter (N = 1)
	parameter (P = 2)
        parameter (S = 3)

c****   MAXLYR is the maximum number of layers
	parameter (MAXLYR = 10)

c****   MONTHS is the number of months
	parameter (MONTHS = 12)

c****   SRFC, SOIL are enumerations for surface, soil
	parameter (SRFC = 1)
	parameter (SOIL = 2)

c****   CPARTS is the number of parts in the grassland/crop system:
c       ABOVE is an enumeration for aboveground
c       BELOW is an enumeration for belowground
	parameter (CPARTS = 2)
	parameter (ABOVE = 1)
	parameter (BELOW = 2)

c****   FPARTS is the number of parts in the forest system:
c       LEAF  =  leaf forest part
c       FROOT =  fine root forest part
c       FBRCH =  fine branch forest part
c       LWOOD =  large wood forest part
c       CROOT =  coarse root forest part
	parameter (FPARTS = 5)
	parameter (LEAF = 1)
	parameter (FROOT = 2)
	parameter (FBRCH = 3)
	parameter (LWOOD = 4)
	parameter (CROOT = 5)

c****   NEWFOR, OLDFOR are the new and old forests array pointers
	parameter (NEWFOR = 1)
	parameter (OLDFOR = 2)

c****   IMIN, IMAX, IVAL are enumerations for indexing arrays
	parameter (IMIN = 1)
	parameter (IMAX = 2)
	parameter (IVAL = 3)

c****   INTCPT, SLOPE are the intercept and slope array pointers
	parameter (INTCPT = 1)
	parameter (SLOPE = 2)

c****   UNLABL, LABELD are the unlabeled, labeled array pointers
	parameter (UNLABL = 1)
	parameter (LABELD = 2)
c****   ISOS is the total number of isotopes (unlabeld, labeled)
	parameter (ISOS = 2)

c****   ESTOR, ESOIL, ENFIX, EFERT are enumerations used only in
c	restrp.f, growth.f, trees.f
	parameter (ESTOR = 1)
	parameter (ESOIL = 2)
	parameter (ENFIX = 3)
	parameter (EFERT = 4)

c****   Constant values
	parameter (PEEDEE=0.0112372)
	parameter (PI=3.141592653589793)
	parameter (PI2=6.283185307179586)
	
	
c******************************************************************
c Labeled common /save_vars/ 
c Note that each variable name includes a suffix identifying the
c procedure which uses the variable.  For example,
c "prevsys_cen_step" is the prevsys variable local to the 
c cen_step() procedure, while "prevsys_readblk" is a different
c prevsys variable local to the readblk() procedure.

      common /save_vars/ prevsys_cen_step, 
     - itab_randu(64), init_randu, ictr_randu,
     - r_ran1(97), ix1_ran1, ix2_ran1, ix3_ran1,
     - prevsys_readblk,
     - plntd_schedl, curcult_schedl, curfert_schedl, 
     - curfire_schedl, curgraz_schedl, curharv_schedl,
     - curirri_schedl, curomad_schedl, curtrm_schedl, 
     - savedfert_schedl,
     - save_varsEnd
      
      integer prevsys_cen_step,
     - itab_randu, ictr_randu,
     - ix1_ran1, ix2_ran1, ix3_ran1,
     - prevsys_readblk,
     - plntd_schedl
      
      logical init_randu
       
      real r_ran1,
     - savedfert_schedl
      
      character*5 curcult_schedl, curfert_schedl, curfire_schedl,
     - curgraz_schedl, curharv_schedl, curirri_schedl, 
     - curomad_schedl, curtrm_schedl
      
      real save_varsStart, save_varsEnd
      equivalence (save_varsStart, prevsys_cen_step)
      save /save_vars/
    
      
