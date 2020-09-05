c*****************************************************************
c	$RCSfile: site.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: site.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/19 23:34:57  conklind
cInitial revision
c
c******************************************************************
       common/site/sitlat,sitlng,sitpot,sand,silt,clay,sitelev,
     - siteEnd
     
       real sitlat,sitlng,sitpot,sand,silt,clay,sitelev
           
       real siteStart, siteEnd
       equivalence (siteStart, sitlat)
       save /site/
