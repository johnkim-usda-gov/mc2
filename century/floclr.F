
c*****************************************************************
c	$RCSfile: floclr.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: floclr.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:26:42  bachelet
cheader
c
c******************************************************************

c		Copyright 1993 Colorado State University
c			All Rights Reserved



	subroutine floclr
c 
c	floclr clears the flow stack; this can be done at any time, and
c	    has the effect of IGNORING any outstanding flows; the idea is
c	    that it should be possible to restart a simulation from a known,
c	    clean state, with no outstanding flows.
c 
	include 'flow.inc'
	
	nflows = 0 

	return
	end
