c*****************************************************************
c	$RCSfile: forrem.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: forrem.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:45:33  conklind
cInitial revision
c
c******************************************************************
c23456789012345678901234567890123456789012345678901234567890123456789012
c**** Forest removal variables
c----     evntyp = type of event
c                  0 for a cutting event (this includes storms)
c                  1 for a fire event
c
c----     fd(ipool) = fraction of fine and coarse roots that die
c
c----     remf(ipool) = removal fraction
c                       ipool indicates which state variable
c
c----     retf(ipool,iel) = return fraction
c                           ipool indicates which state variable
c                           iel indicates C, N, P, or S
c
 	common /forrem/ evntyp,fd(2),mf_fd(2),remf(5),mf_remf(5),
     -    retf(3,4), forremEnd
	integer evntyp
	real fd, mf_fd, remf, mf_remf, retf

        real forremStart, forremEnd
        equivalence (forremStart, evntyp)
        save /forrem/
