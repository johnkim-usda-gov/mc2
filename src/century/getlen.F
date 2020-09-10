      
c*****************************************************************
c	$RCSfile: getlen.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: getlen.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:31:05  bachelet
cheader
c
c******************************************************************
	integer function getlen(name)

      character*(*) name
      integer rlen, max

      max = len(name)
      do 10 rlen = 1, max
         if (name(rlen:rlen) .le. ' ')  goto 20
 10   continue
      rlen = 0
 
 20   if (rlen .eq. 0) then
         getlen = max
      else
         getlen = rlen - 1
      endif

      return
      end
