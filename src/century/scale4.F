
c*****************************************************************
c	$RCSfile: scale4.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: scale4.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:42:52  bachelet
cheader
c
c******************************************************************
c----------------------------------------------------------------------
      subroutine scale4 (f, ic, val1, val2, val3, val4, result)

      real f, result, val, val1, val2, val3, val4
      integer ic
      dimension val(4)

      val(1) = val1
      val(2) = val2
      val(3) = val3
      val(4) = val4

      result = (f * val(ic+1)) + ((1. - f) * val(ic))

      return
      end
