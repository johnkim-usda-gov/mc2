
c*****************************************************************
c	$RCSfile: grazin.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: grazin.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:32:06  bachelet
cheader
c
c******************************************************************
C...GRAZIN.F

      subroutine grazin(tomatch,curgraz)

      character*5 tomatch
      character*5 curgraz

c...Read in the new graze type

      include 'parcp.inc'
      include 'const.inc'

      include 'mapssvar.inc'

c...Local variables
      integer   i, GRAZLNS
      real      temp
      character fromdat*5, name*20, string*80
      character*132 filename


c...Number of lines to read for each graze type
      parameter (GRAZLNS = 11)
      filename = ""
      filename(1:pathlen) = filpath
      filename(pathlen+1:pathlen+8) = 'graz.100'
      open(unit=11, file=filename,status='OLD')

cx      open(unit=11, file='/data/liveoak1/gclim_tools/models/century/model/'//'graz.100',status='OLD')
20    continue
      read(11, 100, end=200) fromdat
      if (tomatch .ne. fromdat) then 
        do 25 i = 1, GRAZLNS 
          read(11, *) temp, name
25      continue
        goto 20
      else
        read(11, *) flgrem, name 
        call ckdata('schedl','flgrem',name)
        read(11, *) fdgrem, name 
        call ckdata('schedl','fdgrem',name)
        read(11, *) gfcret, name 
        call ckdata('schedl','gfcret',name)
        read(11, *) gret(N), name 
        call ckdata('schedl','gret',name)
        read(11, *) gret(P), name 
        call ckdata('schedl','gret',name)
        read(11, *) gret(S), name 
        call ckdata('schedl','gret',name)
        read(11, *) temp, name
        grzeff = int(temp)
        call ckdata('schedl','grzeff',name)
        read(11, *) fecf(N), name 
        call ckdata('schedl','fecf',name)
        read(11, *) fecf(P), name 
        call ckdata('schedl','fecf',name)
        read(11, *) fecf(S), name 
        call ckdata('schedl','fecf',name)
        read(11, *) feclig, name 
        call ckdata('schedl','feclig',name)
        close(11)
        curgraz = tomatch 
      endif


      return

100   format(a5)

200   continue
      call message('   Error reading in values from the graz.100 file.')
      string = '   Looking for type: ' // tomatch
      call message(string)
      STOP

      end
