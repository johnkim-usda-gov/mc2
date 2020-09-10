
c*****************************************************************
c	$RCSfile: irrigt.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: irrigt.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:49:29  bachelet
cheader
c
c******************************************************************
C...IRRIGT.F

      real function irrigt(month)

      integer   month

c...Simulate irrigation
 
      include 'param.inc'
      include 'parcp.inc'
      include 'plot1.inc'
      include 'wth.inc'

c...Local variables
      integer layr
      real toth2o


c...Check tave so that irrigation water 
c     is not added as snow
      if (tave .le. 0.0) then
        irract = 0.0

c...Add amount given by user
      else
        if (auirri .eq. 0) then
          irract = irramt

c...Add amount automatically to field capacity
        else if (auirri .eq. 1) then
          if (avh2o(1)/awhc .le. fawhc) then
            toth2o = prcurr(month)
            do 10 layr = 1, nlaypg
              toth2o = toth2o + asmos(layr)
10          continue
            irract = max(twhc - toth2o,0.0)
          else
            irract = 0.0
          endif

c...Add amount automatically to nominated amount
        else if (auirri .eq. 2) then
          if (avh2o(1)/awhc .le. fawhc) then
            irract = irraut
          else
            irract = 0.0
          endif

c...Add amount automatically to field capacity plus PET
        else if (auirri .eq. 3) then
          if (avh2o(1)/awhc .le. fawhc) then
            toth2o = prcurr(month)
            do 20 layr = 1, nlaypg
              toth2o = toth2o + asmos(layr)
20          continue
            irract = max(twhc + pet - toth2o,0.0)
          else
            irract = 0.0
          endif
        endif
      endif

      irrtot = irrtot + irract
      irrigt = irract

      return
      end
