      
c*****************************************************************
c	$RCSfile: prcgrw.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: prcgrw.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:26:23  bachelet
cheader
c
c******************************************************************
	subroutine prcgrw(imnth)

      integer imnth

c...Compute growing season precipiation

      include 'const.inc'
      include 'parcp.inc'
      include 'wth.inc'

c...Local variables
      integer m


      grwprc = 0.

      do 10 m = imnth, MONTHS
        grwprc = grwprc + prcurr(m)
10    continue

      do 20 m = 1, imnth - 1
        grwprc = grwprc + prcnxt(m)
20    continue

      return
      end
