      
c*****************************************************************
c	$RCSfile: update_sched.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: update_sched.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:52:45  bachelet
cheader
c
c******************************************************************
	subroutine update_sched(crtyr)

c-----C. Daly - update_sched moves this year's scheduled events
c-----from storage arrays into the working arrays.

      include 'chrvar.inc'
      include 'schvar.inc'
      include 'timvar.inc'
      include 't0par.inc'

      include 'mapssvar.inc'

      integer i, crtyr

      ttlind = 0

      do 40 i = 1, ttlind_stor

c-------Check to see if event is to occur this year.  timstor carries
c-------the nth repeat year each event is to occur.  Events that
c-------had return intervals greater than years_to_run have
c-------timstor values of zero, so skip them to avoid mod errors.
  
        if (timstor(i, 1) .le. 0.) go to 40

        if(mod(crtyr, timstor(i, 1)) .lt. dt * 0.5) then

c---------If event occurs this year, put it in.
          ttlind = ttlind + 1
          cmdary(ttlind) = cmdstor(i)
          typary(ttlind) = typstor(i)
          fltary(ttlind, 1) = fltstor(i,1)
          timary(ttlind, 1) = crtyr
          timary(ttlind, 2) = timstor(i,2)

cx          write(*,100) ttlind, timary(ttlind, 1),timary(ttlind, 2),
cx     .                 cmdary(ttlind),typary(ttlind),fltary(ttlind, 1)

        endif

40    continue

      return

100   format(' update_sched # ',i4,'  -> ',i5,i3,' "',a4,'" "',a5,
     .       '"',g12.5)
      end
