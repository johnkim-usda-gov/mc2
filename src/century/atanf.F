
c*****************************************************************
c	$RCSfile: atanf.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: atanf.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 19:59:14  bachelet
cheader
c
c******************************************************************

c		Copyright 1993 Colorado State University
c			All Rights Reserved



	real function atanf(x,a,b,c,d)
	real x, a, b, c, d
c 
c  ******************** flowlib ********************* 
c 
c			(run-time sub-set of modaid, exclusive of modctl) 
c 
c			release 1.0  (first formal release of modaid) 
c 
c			james m. vevea
c			natural resource ecology lab
c			colorado state university 
c			fort collins, colorado	80523 
c 
c      this routine is functionally equivalent to the routine of the
c		same name, described in the publication: 
c 
c		some graphs and their functional forms 
c		technical report no. 153 
c		william parton and georgo innis	(1972) 
c		natural resource ecology lab.
c		colorado state university
c		fort collins, colorado  80523
c 
	include 'pi.inc'
	atanf = b+ (c / PI) * atan(PI * d * (x - a))
	return 
	end
