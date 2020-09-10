c*****************************************************************
c	$RCSfile: wth.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: wth.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/19 23:38:32  conklind
cInitial revision
c
c******************************************************************
      common/wth/prcurr(12),prcnxt(12),tmn2m(12),tmx2m(12),
     -wthEnd
     
      real prcurr, prcnxt, tmn2m, tmx2m
           
       real wthStart, wthEnd
       equivalence (wthStart, prcurr(1))
       save /wth/
