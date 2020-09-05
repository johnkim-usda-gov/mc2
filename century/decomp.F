
c*****************************************************************
c	$RCSfile: decomp.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: decomp.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:17:08  bachelet
c.header
c
c******************************************************************
C...DECOMP.F

c...Decomposition Submodel (rewritten by vek 04/91)

	subroutine decomp(dtm,decsys)

	real      dtm
	integer   decsys

	include 'const.inc'

 

c                   ********** LITTER **********

c...Decompose structural and metabolic components for surface and soil.
c        print *, 'calling litdec() now'
        call litdec(dtm)


c                   *********** WOOD ***********

c...If the system is a forest or savanna...
c     Decompose dead fine branches, large wood, and coarse roots.
c     Dead fine roots are in the soil structural compartment.

        if (decsys .eq. FORSYS) then
c          print *, 'calling woodec() now'
          call woodec(dtm)
        endif

c                 ***** SOIL ORGANIC MATTER *****

c...Decompose som1 (surface and soil), som2, and som3.

c        print *, 'calling somdec() now'
        call somdec(dtm)

        return
        end
