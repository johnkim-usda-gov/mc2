      
c*****************************************************************
c	$RCSfile: repeat_event.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: repeat_event.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 19:43:51  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:36:05  bachelet
cheader
c
c******************************************************************
	subroutine repeat_event()

c-----C. Daly 2-7-95 - For one_step, use yrnum as a return interval.
c-----So, a yrnum of 4 means to repeat the event every 4th year, 
c-----starting with year 4; yrnum=1 means repeat event every year 
c-----starting with yr 1.

c-----C. Daly 5-3-96
c-----Broke out of readblk into subroutine.

      include 'schvar.inc'
      include 'chrvar.inc'
      include 'timvar.inc'
      include 't0par.inc'

      integer year_new, ttlind_orig

c-----If return interval is greater than repeat years, skip event.
      if(timary(ttlind, 1) .gt. rptyrs) then
        write(*,*) ' Return interval of ',timary(ttlind, 1),
     .         ' for ',cmdary(ttlind), 'is greater than ',
     .         ' simulation years of ',rptyrs, '. No scheduling done.'
        return
      endif

      ttlind_orig = ttlind
 
      do 90 year_new = timary(ttlind, 1) + 1, rptyrs

        if(mod(year_new, timary(ttlind_orig, 1)) .lt. dt * 0.5) then
          ttlind = ttlind + 1
          timary(ttlind, 1) = year_new
          timary(ttlind, 2) = timary(ttlind_orig, 2)
          cmdary(ttlind) = cmdary(ttlind_orig)
          typary(ttlind) = typary(ttlind_orig)
          fltary(ttlind, 1) = fltary(ttlind_orig, 1)

          write(*,'(a,i5,a,i5,i3,a,a4,3a,g12.5)') ' Read event #', 
     $    ttlind,'   ->', timary(ttlind, 1),timary(ttlind, 2),' "',
     $    cmdary(ttlind),'" "',typary(ttlind),'"',fltary(ttlind, 1)
        endif

90      continue

        return
        end
