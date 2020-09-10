
c*****************************************************************
c	$RCSfile: pevap.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: pevap.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:23:15  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

C...PEVAP.F

      real function pevap(mnum)
 
      integer   mnum

      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'site.inc'

c...Local variables
      integer   i
      real      avgmth(MONTHS), e, highest, lowest, monpet,
     $          ra, t, td, tm, tr


      sitelev = 0.0

c...Determine max and min temperatures
      do 10 i = 1, 12
        avgmth(i) = (maxtmp(i) + mintmp(i)) / 2.0
10    continue
      highest = avgmth(1)
      lowest = avgmth(1)
      do 20 i = 2, 12
        highest = max(highest, avgmth(i))
        lowest = min(lowest, avgmth(i))
20    continue
      lowest = MAX(lowest, -10.)

c...Determine average temperature range
      ra = abs(highest - lowest)

c...Temperature range calculation
      tr = maxtmp(mnum) - MAX(-10., mintmp(mnum))

      t = tr/2.+mintmp(mnum)
      tm = t+0.006*sitelev
      td = 0.0023*sitelev+0.37*t+0.53*tr+.35*ra-10.9
      e = ((700.*tm/(100.-abs(sitlat)))+15.*td)/(80.-t)
      monpet = (e*30.)/10.
      if (monpet .lt. 0.5) then
        monpet = 0.5
      endif

c...fwloss(4) is a modifier for PET loss.   vek may90
      pevap = monpet*fwloss(4)

      return
      end
