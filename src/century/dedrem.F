
c*****************************************************************
c	$RCSfile: dedrem.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: dedrem.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:17:37  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

c...DEDREM
      subroutine dedrem(accum)

c...Removal of dead wood due to cutting or fire in a forest.

c...Called from:	frem

      include 'const.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'zztim.inc'

c...Argument Declarations
      real      accum(ISOS)

c...Local Variables
      integer	iel
      real	closs, eloss(MAXIEL)

c...Remove dead FINE BRANCHES

      closs = remf(4) * wood1c
      tcrem = tcrem + closs
      call csched(closs,wd1cis(LABELD),wood1c,
     $             wd1cis(UNLABL),csrsnk(UNLABL),
     $             wd1cis(LABELD),csrsnk(LABELD),
     $             1.0,accum)

      do 10 iel = 1, nelem
c... Bachelet 11-5-97
	if (wood1c .ne. 0.0) then
		eloss(iel) = closs * (wood1e(iel) / wood1c)
	else
		eloss(iel) = 0.0
	endif	
        terem(iel) = terem(iel) + eloss(iel)
        call flow(wood1e(iel),esrsnk(iel),time,eloss(iel))
 10   continue

c...Remove dead LARGE WOOD

      closs = remf(5) * wood2c
      tcrem = tcrem + closs
      call csched(closs,wd2cis(LABELD),wood2c,
     $             wd2cis(UNLABL),csrsnk(UNLABL),
     $             wd2cis(LABELD),csrsnk(LABELD),
     $             1.0,accum)

      do 20 iel = 1, nelem
c... Bachelet 11-5-97
	if(wood2c .ne. 0.0) then
		eloss(iel) = closs * (wood2e(iel) / wood2c)
	else
		eloss(iel) = 0.0
	endif
        terem(iel) = terem(iel) + eloss(iel)
        call flow(wood2e(iel),esrsnk(iel),time,eloss(iel))
 20   continue

      return
      end
