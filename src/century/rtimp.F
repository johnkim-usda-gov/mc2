
c*****************************************************************
c	$RCSfile: rtimp.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: rtimp.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/28 17:55:01  bachelet
cheader
c
c******************************************************************
C...RTIMP.F

	real function rtimp(riint, rictrl, rootc)

	real      riint, rictrl, rootc

c...This function calculates and returns a value between 0-1 which is the impact
c     of root biomass on available nutrients.  It is used in the calculation of
c     total plant production in RESTRP and SRESTRP.


c...Called From:       GROWTH
c                      SIMSOM
c                      TREES
	if((-rictrl * rootc * 2.5).lt. -10.0) then
		rtimp = 1.0
	else
        	rtimp = (1.0 - riint * exp(-rictrl * rootc * 2.5))
	endif

c... D. Bachelet catching underflows

        return
        end
