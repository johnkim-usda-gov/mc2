
c*****************************************************************
c	$RCSfile: ramp.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: ramp.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:29:00  bachelet
cheader
c
c******************************************************************
C...RAMP.F

	real function ramp(x, x1, y1, x2, y2)

	real      x, x1, y1, x2, y2

c...This function models a "ramp":
c              /-----
c             /
c       -----/


c...Function declarations
	real      line
	external  line


	if (x .lt. x1) then 
          ramp = y1
	else if (x .gt. x2) then
          ramp = y2
        else
          ramp = line(x, x1, y1, x2, y2)
c...where line = (y2 - y1) / (x2 - x1) * (x - x2) + y2
	endif

	return
	end
