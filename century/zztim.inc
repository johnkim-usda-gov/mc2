c*****************************************************************
c	$RCSfile: zztim.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: zztim.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:49:45  conklind
cInitial revision
c
c******************************************************************
c**** time is in a separate common block for compatability with TIME0
c
      common /zztim/ time,
     -  zztimEnd
     
      real time

      real zztimStart, zztimEnd
      equivalence (zztimStart, time)
      save /zztim/

