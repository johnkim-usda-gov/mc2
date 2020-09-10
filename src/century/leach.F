	
c*****************************************************************
c	$RCSfile: leach.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: leach.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:52:28  bachelet
cheader
c
c******************************************************************
	subroutine leach(amov, nelem, nlayer, minerl, minlch, frlech,
     $	                 stream, basef, stormf)

	integer   nelem, nlayer
	real      amov(10), minerl(10,3), minlch, frlech(3),
     $            stream(8), basef, stormf
 
c...This routine computes the leaching of nitrogen, phosphorus,
c     and sulfur.
c...Written 2/92 -rm

c...Called From:	SIMSOM

c...amtlea:    amount leached
c...linten:    leaching intensity
c...strm:      storm flow
c...base:      base flow

c...Local variables
        integer   iel, lyr, nxt
        real      amtlea(10), base, linten, strm


	do 20 iel = 1,nelem
	   strm = 0.0
	   base = 0.0
	   stream(iel+1) = 0.0
	   do 10 lyr = 1,nlayer
	      amtlea(lyr) = 0.0
	      nxt = lyr + 1

c...AMOV > 0. indicates a saturated water flow out of layer lyr
	      if ((amov(lyr) .gt. 0.0) .and. 
     $		  (minerl(lyr,iel) .gt. 0.0)) then
	         linten = min(1.0 - (minlch - amov(lyr)) / 
     $				minlch, 1.0)
		 linten = max(linten, 0.0)
		 amtlea(lyr) = frlech(iel) * minerl(lyr,iel) * linten

c...If you are at the bottom layer, compute storm flow.
		 if (lyr .eq. nlayer) then
                   strm = amtlea(lyr) * stormf
                 endif
		 minerl(lyr,iel) = minerl(lyr,iel) - amtlea(lyr)
		 minerl(nxt,iel) = minerl(nxt,iel) + (amtlea(lyr)
     $				 - strm)
	      endif
10	   continue

c...Compute base flow and mineral stream flows.
c...Make sure first that there is something in that layer otherwise
c...base can become extremely small and create problems DB 10-1-98
           if(minerl(nxt,iel).gt. 0.01) then
	   	base = minerl(nxt,iel) * basef
           else
		base = 0.0
           endif
c...  another trap to keep minerl >=0  DB 10-13-98
	   if(minerl(nxt,iel) .lt. base) base = minerl(nxt,iel)

	   minerl(nxt,iel) = minerl(nxt,iel) - base

c...Note: stream flow indices differ from mineral element 
c         indices by 1 (eg  stream(2) is stream flow for nitrogen).
	   stream(iel+1) = strm + base
 20	continue

	return
	end
