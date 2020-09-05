
c*****************************************************************
c	$RCSfile: cmplig.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: cmplig.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:01:45  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


        subroutine cmplig(cursys,fligni,wdlig,pltlig)

        include 'const.inc'
        include 'param.inc'
        include 'wth.inc'

        include 'zztim.inc'

        integer  cursys
        real     fligni(2,CPARTS), wdlig(FPARTS), pltlig(CPARTS)

c...Compute plant lignin; returns the fraction of residue which will
c     lignin.

c...Local variables
        integer  m
        real     arain

c...Cursys tells whether a crop, forest or savanna system is being simulated.
        if (cursys .eq. CRPSYS .or. cursys .eq. SAVSYS) then
c...crop or savanna system: lignin contend depends on annual rainfall
          arain = 0.
          do 10 m = 1, MONTHS
            arain = arain + prcurr(m)
10        continue
          if (arain .eq. 0.) then
            do 20 m = 1, MONTHS
               arain = arain + precip(m)
20            continue
          endif
        endif




        if (cursys .eq. CRPSYS) then
c...Crop/grass
          pltlig(ABOVE)=fligni(INTCPT,ABOVE)+fligni(SLOPE,ABOVE)*arain
          pltlig(BELOW)=fligni(INTCPT,BELOW)+fligni(SLOPE,BELOW)*arain
        
        else if (cursys .eq. FORSYS) then
c...Forest system (leaves,fine roots=surface,soil)
          pltlig(SRFC)=wdlig(LEAF)
          pltlig(SOIL)=wdlig(FROOT)

        else if (cursys .eq. SAVSYS) then
c...Savanna

          pltlig(SRFC) = (wdlig(LEAF)+fligni(INTCPT,SRFC) +
     $                    fligni(SLOPE,SRFC) * arain) / 2.0
          pltlig(SOIL) = (wdlig(FROOT)+fligni(INTCPT,SOIL) +
     $                    fligni(SLOPE,SOIL) * arain) / 2.0
        endif

c...Check range for pltlig; the hard-coded values should be replaced with
c     parameter values someday
        pltlig(SRFC) = max(0.02, pltlig(SRFC))
	pltlig(SRFC) = min(0.5, pltlig(SRFC))
	pltlig(SOIL) = max(0.02, pltlig(SOIL))
	pltlig(SOIL) = min(0.5, pltlig(SOIL))

c       if (pltlig(SRFC).lt.0.02 .or. pltlig(SRFC).gt.0.5) then
c         if (pltlig(SOIL).lt.0.02 .or. pltlig(SOIL).gt.0.5) then
c           call message('   Plant lignin is out of the range')
c           call message('   of 0.02 < pltlig < 0.5')
c           STOP
c         endif
c       endif

        return
        end
