c*****************************************************************
c	$RCSfile: pass.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: pass.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:48:13  conklind
cInitial revision
c
c******************************************************************
	common/pass/pass1,
     -    passEnd
     
	logical pass1

      real passStart, passEnd
      equivalence (passStart, pass1)
      save /pass/
