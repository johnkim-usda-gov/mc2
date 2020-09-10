      
c*****************************************************************
c	$RCSfile: message.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: message.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:54:51  bachelet
cheader
c
c******************************************************************
	subroutine message(string)

      character*(*) string

c...Write out the string on standard output
c     Done in a separate subroutine so that specific OS methods of
c     writing to standard output need only to make changes here.

      write(*,*) string

      return
      end                                                                       
