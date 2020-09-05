c*****************************************************************
c	$RCSfile: flow.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: flow.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2004/08/30 16:45:49  conklind
cAdd zzzfloStart and zzzfloEnd.
c
cRevision 1.1  2003/11/20 21:45:24  conklind
cInitial revision
c
c******************************************************************
c
c	These comments relate to the entire 'flow' system:
c
c	The variables to be moved via flow MUST remain active at
c	    ALL times from the call to flow to the call to flowup;
c	    More specifically, they should normally be in common blocks
c	    which are defined in a routine (probably your main program)
c	    which remains active during the entire process!  If this isn't
c	    done, I can't tell, and you will get very strange answers!
c
	integer LENFST
	parameter (LENFST=250)
	real WLEN
	parameter (WLEN=4.0)
	
	save /zzzflo/

	common /zzzflo/ nflows,flostk(4,LENFST),dummy(1),zzzfloEnd

        real nflows, flostk, dummy, zzzfloStart, zzzfloEnd
	integer iflost(4,LENFST)
	equivalence (flostk(1,1),iflost(1,1))
	equivalence (zzzfloStart, nflows)
