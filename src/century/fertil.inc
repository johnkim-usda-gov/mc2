c*****************************************************************
c	$RCSfile: fertil.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: fertil.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:45:18  conklind
cInitial revision
c
c******************************************************************
      common/fertil/aufert,feramt(3),
     -  fertilEnd 
      real aufert,feramt

      real fertilStart, fertilEnd
      equivalence (fertilStart, aufert)
      save /fertil/
