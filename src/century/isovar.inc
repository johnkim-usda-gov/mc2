c*****************************************************************
c	$RCSfile: isovar.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: isovar.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:47:24  conklind
cInitial revision
c
c******************************************************************
c**** c14frc tells what fraction of the carbon in new plant tissue is labeled.
      common/isovar/cisofr,cisotf,
     -  isovarEnd      
      real cisofr, cisotf

      real isovarStart, isovarEnd
      equivalence (isovarStart, cisofr)
      save /isovar/
