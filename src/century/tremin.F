
c*****************************************************************
c	$RCSfile: tremin.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: tremin.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:51:34  bachelet
cheader
c
c******************************************************************
C...TREMIN.F

      subroutine tremin(tomatch,curtrm)

      character*5 tomatch
      character*5 curtrm

c...Read in the new tree removal type

      include 'forrem.inc'
      include 'const.inc'

      include 'mapssvar.inc'

c...Local variables
      integer   i, TREMLNS
      real      temp
      character fromdat*5, name*20, string*80
      character*132 filename

c...Number of lines to read for each tree removal type
      parameter (TREMLNS = 20)
      filename = ""
      filename(1:pathlen) = filpath
      filename(pathlen+1:pathlen+8) = 'trem.100'
      open(unit=11, file=filename,status='OLD')

cx      open(unit=11, file='/data/liveoak1/gclim_tools/models/century/model/'//'trem.100',status='OLD')
20    continue
      read(11, 100, end=200) fromdat
      if (tomatch .ne. fromdat) then 
        do 25 i = 1, TREMLNS 
          read(11, *) temp, name
25      continue
        goto 20
      else
        read(11, *) temp, name
        evntyp = int(temp)
        call ckdata('schedl','evntyp',name)
        read(11, *) remf(LEAF), name 
        call ckdata('schedl','remf',name)
        read(11, *) remf(FROOT), name 
        call ckdata('schedl','remf',name)
        read(11, *) remf(FBRCH), name 
        call ckdata('schedl','remf',name)
        read(11, *) remf(LWOOD), name 
        call ckdata('schedl','remf',name)
        read(11, *) remf(CROOT), name 
        call ckdata('schedl','remf',name)
        read(11, *) fd(1), name 
        call ckdata('schedl','fd',name)
        read(11, *) fd(2), name 
        call ckdata('schedl','fd',name)
        read(11, *) retf(1,1), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(1,2), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(1,3), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(1,4), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(2,1), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(2,2), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(2,3), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(2,4), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(3,1), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(3,2), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(3,3), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(3,4), name 
        call ckdata('schedl','retf',name)
        close(11)
        curtrm = tomatch 
      endif


      return

100   format(a5)

200   continue
      call message('   Error reading in values from the trem.100 file.')
      string = '   Looking for type: ' // tomatch
      call message(string)
      STOP

      end
