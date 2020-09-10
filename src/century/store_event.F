      
c*****************************************************************
c	$RCSfile: store_event.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: store_event.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:47:51  bachelet
cheader
c
c******************************************************************
	subroutine store_event()

c-----C. Daly 10-25-96 - stores sceduled events for use during
c-----one-step operation.

      include 'schvar.inc'
      include 'chrvar.inc'
      include 'timvar.inc'
      include 't0par.inc'

      include 'mapssvar.inc'

c-----If return interval is greater than repeat years, skip event.
      if(timary(ttlind, 1) .gt. rptyrs) then
        write(*,*) ' Return interval of ',timary(ttlind, 1),
     .         ' for ',cmdary(ttlind), 'is greater than ',
     .         ' simulation years of ',rptyrs, '. No scheduling done.'
        return
      endif
 
          timstor(ttlind, 1) = timary(ttlind, 1)
          timstor(ttlind, 2) = timary(ttlind, 2)
          cmdstor(ttlind) = cmdary(ttlind)
          typstor(ttlind) = typary(ttlind)
          fltstor(ttlind, 1) = fltary(ttlind, 1)
          ttlind_stor = ttlind

          if (diags) write(*,'(a,i5,a,i5,i3,a,a4,3a,g12.5)')  
     $    ' Storing event #', ttlind,'   ->', timary(ttlind, 1),
     $    timary(ttlind, 2),' "',cmdary(ttlind),'" "',
     $    typary(ttlind),'"',fltary(ttlind, 1)

        return
        end
