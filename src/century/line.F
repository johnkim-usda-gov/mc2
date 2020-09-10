	
c*****************************************************************
c	$RCSfile: line.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: line.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:53:00  bachelet
cheader
c
c******************************************************************
	real function line(x, x1, y1, x2, y2)
	real x, x1, y1, x2, y2

c****	This function is the generic equation of a line from
c	two points.
c	CALLED FROM:	co2eff.f

	line = (y2 - y1) / (x2 - x1) * (x - x2) + y2

	return
	end
