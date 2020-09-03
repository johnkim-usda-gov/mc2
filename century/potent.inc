c*****************************************************************
c	$RCSfile: potent.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: potent.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:48:35  conklind
cInitial revision
c
c******************************************************************

c               Copyright 1993 Colorado State University
c                       All Rights Reserved


	common/potent/agp,tgprod,pcropc,pforc,eup(5,3),
     -    potentEnd
     
	real agp,tgprod,pcropc,pforc,eup

      real potentStart, potentEnd
      equivalence (potentStart, agp)
      save /potent/
