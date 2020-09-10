      
c*****************************************************************
c	$RCSfile: weathr.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: weathr.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:54:12  bachelet
cheader
c
c******************************************************************
	subroutine weathr(precip,prcstd,prcskw,mintmp,maxtmp)

      real      precip(12), prcstd(12), prcskw(12),
     $          mintmp(12), maxtmp(12)

c...Determines the current year values for precipitation and
c     temperature and next year values for precipitation for predicting
c     production potential.

      include 'chrvar.inc'
      include 'const.inc'
      include 'wth.inc'
      
c...Function declarations
      real      anorm, sknorm
      external  anorm, sknorm

c...Local variables
      integer   m
      character label*4
	  
c...Current year precipitation, regardless of wthr type
      do 10 m = 1, MONTHS
        prcurr(m) = prcnxt(m)
10    continue

        do 20 m = 1, MONTHS
          prcurr(m) = precip(m)
          prcnxt(m) = precip(m)
20      continue
        return


      end
