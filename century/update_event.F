      
c*****************************************************************
c	$RCSfile: update_event.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: update_event.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:52:08  bachelet
cheader
c
c******************************************************************
	subroutine update_event(cmd_string,typ_string,fltnum,crtyr,
     .               mon,tevent)

c-----C. Daly - overwrite existing events or add new ones
c-----C. Daly 5-3-96.  Update_event subroutine broken out from 
c-----cen_step.

      include 'chrvar.inc'
      include 'schvar.inc'
      include 'timvar.inc'

      character cmd_string*4
      character typ_string*5
      integer flag
      integer tevent
      integer crtyr
      integer fltnum
      integer mon

c-----Search for matching event name, then flag event number.
      flag = 0
      do 50 tevent = 1, ttlind
        if (cmdary(tevent) .eq. cmd_string) flag = tevent
50    continue
      if ( flag .ne. 0 ) then
        tevent = flag
      else
        ttlind = ttlind + 1
        tevent = ttlind
      endif

c-----Overwrite event with new information, or add to end of list.
      cmdary(tevent) = cmd_string
      typary(tevent) = typ_string
      fltary(tevent, 1) = fltnum
      timary(tevent, 1) = crtyr
      timary(tevent, 2) = mon

c      write(*,100) tevent, timary(tevent, 1),timary(tevent, 2),
c     .  cmdary(tevent),typary(tevent),fltary(tevent, 1)
      return

100   format(' update_event # ',i4,'  -> ',i5,i3,' "',a4,'" "',a5,
     .       '"',g12.5)
      end
