c*****************************************************************
c	$RCSfile: pheno.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: pheno.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:48:26  conklind
cInitial revision
c
c******************************************************************

c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      common/pheno/dayhrs, hrsinc,
     -  phenoEnd
c
c     dayhrs - # of hours in a day based on julian date
c     hrsinc - indicates if the length of days in increasing
c              or decreasing

      real dayhrs
      logical hrsinc

      real phenoStart, phenoEnd
      equivalence (phenoStart, dayhrs)
      save /pheno/


