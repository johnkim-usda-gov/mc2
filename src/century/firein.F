
c*****************************************************************
c	$RCSfile: firein.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: firein.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:25:18  bachelet
cheader
c
c******************************************************************
C...FIREIN.F

      subroutine firein(tomatch,curfire)

      character*5 tomatch
      character*5 curfire

c...Read in the new fire type
c
      include 'parcp.inc'
      include 'const.inc'

      include 'mapssvar.inc'

c...Local variables
      integer   i, FIRELNS
      real      temp
      character fromdat*5, name*20, string*80
      character*132 filename


c...Number of lines to read for each fire type
      parameter (FIRELNS = 10)
      filename = ""
      filename(1:pathlen) = filpath
      filename(pathlen+1:pathlen+8) = 'fire.100'
      open(unit=11, file=filename,status='OLD')

cx      open(unit=11, file='/data/liveoak1/gclim_tools/models/century/model/'//'fire.100',status='OLD')
20    continue
      read(11, 100, end=200) fromdat
      if (tomatch .ne. fromdat) then 
        do 25 i = 1, FIRELNS 
          read(11, *) temp, name
25      continue
        goto 20
      else
        read(11, *) flfrem, name 
        call ckdata('schedl','flfrem',name)
        read(11, *) fdfrem(1), name 
        call ckdata('schedl','fdfrem',name)
        read(11, *) fdfrem(2), name 
        call ckdata('schedl','fdfrem',name)
        read(11, *) ffcret, name 
        call ckdata('schedl','ffcret',name)
        read(11, *) fret(N), name 
        call ckdata('schedl','fret',name)
        read(11, *) fret(P), name 
        call ckdata('schedl','fret',name)
        read(11, *) fret(S), name 
        call ckdata('schedl','fret',name)
        read(11, *) frtsh, name 
        call ckdata('schedl','frtsh',name)
        read(11, *) fnue(ABOVE), name 
        call ckdata('schedl','fnue',name)
        read(11, *) fnue(BELOW), name 
        call ckdata('schedl','fnue',name)
        close(11)
        curfire = tomatch 
      endif


      return

100   format(a5)

200   continue
      call message('   Error reading in values from the fire.100 file.')
      string = '   Looking for type: ' // tomatch
      call message(string)
      STOP

      end
