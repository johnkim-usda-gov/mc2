      
c*****************************************************************
c	$RCSfile: litburn.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: litburn.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:53:27  bachelet
cheader
c
c******************************************************************
	subroutine litburn(litrme)

c...Simulate removal of litter by fire for the month.

C...Called from:      grem.f
c                     frem.f
 
C NOTES:
c...fdfrem(2) - is fdrem(2) in grem.f (fraction of litter layer
c               removed by fire)
c...litrme(MAXIEL) is an accumulator of elements to be returned to
c       system after a fire.  All carbon is assumed lost.


      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

      real     litrme(MAXIEL)

c...Local variables
      integer  iel
      real     stremc,streme(MAXIEL),metrmc,metrme(MAXIEL) 



c...Remove carbon fraction of structural component of litter

       if (strucc(SRFC) .gt. 0.) then
          stremc = fdfrem(2) * strucc(SRFC)
          sdrema = sdrema + stremc
          call csched(stremc,strcis(SRFC,LABELD),strucc(SRFC),
     $                strcis(SRFC,UNLABL),csrsnk(UNLABL),
     $                strcis(SRFC,LABELD),csrsnk(LABELD),
     $                1.0,sdrmai)
         
c...Remove element fraction of structural component of litter

          do 40 iel = 1, nelem
             streme(iel) = stremc*struce(SRFC,iel)/strucc(SRFC)
             litrme(iel) = litrme(iel)+streme(iel)
c...Corrected SDRMAE from SDRMAI 8-30-90  -rm
             sdrmae(iel) = sdrmae(iel)+streme(iel)

             call flow(struce(1,iel),esrsnk(iel),time,streme(iel))
40        continue

       endif


c...Remove carbon fraction of metabolic component of litter

       if (metabc(SRFC) .gt. 0.) then
          metrmc = fdfrem(2) * metabc(SRFC)
          sdrema = sdrema + metrmc
          call csched(metrmc,metcis(SRFC,LABELD),metabc(SRFC),
     $                metcis(SRFC,UNLABL),csrsnk(UNLABL),
     $                metcis(SRFC,LABELD),csrsnk(LABELD),
     $                1.0,sdrmai)
     

c...Remove element fraction of metabolic component of litter

          do 50 iel = 1, nelem
             metrme(iel) = metrmc*metabe(1,iel)/metabc(1)
             litrme(iel) = litrme(iel) + metrme(iel)
c...Corrected SDRMAE from SDRMAI 8-30-90  -rm
             sdrmae(iel) = sdrmae(iel) + metrme(iel)

             call flow(metabe(SRFC,iel),esrsnk(iel),time,metrme(iel))
50        continue

       endif
     

       return
       end
