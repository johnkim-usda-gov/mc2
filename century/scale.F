
c*****************************************************************
c	$RCSfile: scale.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: scale.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:42:14  bachelet
cheader
c
c******************************************************************
c----------------------------------------------------------------------
      subroutine scale (f, val1, val2, result)

      real f, val1, val2, result


      result = (f * val1) + ((1. - f) * val2)

      return
      end
