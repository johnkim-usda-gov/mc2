
c*****************************************************************
c	$RCSfile: frespr.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: frespr.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:29:58  bachelet
cheader
c
c******************************************************************
C...FRESPR.F

       real function frespr(temp,nitrog)

       real      temp, nitrog


c...This function calculates the maintenance respiration per month
c     for the forest materials (fine root, coarse root, fine branch,
c     large wood, and leaf) using an equation derived by Mike Ryan.
c     Maintenance respiration and tissue N content are strongly
c     correlated.

c...Ref:  A simple method for estimating gross carbon budget
c           for vegetation in forest ecosystems.  (Appendix A)
c           Michael G. Ryan
c           Tree Physiology, February 1, 1991
c      
c      temp   :  tave - for leaf and stem, stemp - for root
c      nitrog :  rleave(1), rlwode(1), fbrche(1), croote(1), or froote(1)
c

       frespr =  ((.0106 / 4.0) * (12. / 14.) * 24 * 30 *
     $           exp((alog(2.0) / 10.) * temp) * nitrog) * 1.1

       return
       end
