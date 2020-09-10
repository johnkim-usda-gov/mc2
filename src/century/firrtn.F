
c*****************************************************************
c	$RCSfile: firrtn.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: firrtn.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:25:46  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

c...FIRRTN
      subroutine firrtn()

c...Elemental return from a fire event.

c...Called from:	frem

      include 'const.inc'
      include 'dovars.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'zztim.inc'

c...Local Variables
      integer	iel
      real	egain(MAXIEL)

c...LITTER BURN

c...Only burn litter if a forest system.  Litter is burned in
c   grem.f for the savanna system.

      do 10 iel = 1, nelem
	egain(iel) = 0.0
 10   continue

      if (dofire(FORSYS)) then
	call litburn(egain)
      endif

c...Return from TREE compartments
c   Carbon return is usually 0.  It is ignored since it 
c   would be returned as charcoal.  N, P, and S returns
c   go to the top layer of minerl.  EGAIN will contain
c   the total returns for N, P, and S across pools.

      do 20 iel = 1, nelem
	 egain(iel) = egain(iel) +
     $                remf(1) * retf(1,iel+1) * rleave(iel) +
     $                remf(2) * retf(2,iel+1) * fbrche(iel) +
     $                remf(3) * retf(3,iel+1) * rlwode(iel) +
     $                remf(4) * retf(2,iel+1) * wood1e(iel) +
     $                remf(5) * retf(3,iel+1) * wood2e(iel) +
     $                remf(3) * retf(3,iel+1) * forstg(iel)
       	 call flow(esrsnk(iel),minerl(1,iel),time,egain(iel))
 20   continue

      return
      end
