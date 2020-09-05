c*****************************************************************
c	$RCSfile: sitein.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: sitein.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.1  2003/12/04 00:19:05  conklind
cInitial revision
c
c******************************************************************
C...SITEIN.F

      subroutine sitein(ivopt)

      logical   ivopt

c...Read the parameter file in the form needed by time0.

      include 'const.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot3.inc'
c...schvar.inc is included here because it has the code_flags
      include 'schvar.inc'
      include 'site.inc'
      include 'wth.inc'
      include 'mapssvar.inc'

c...Local variables
      integer   i, j
      character*6 name
      real      temp, dum(MONTHS)

c-----C. Daly - Commenting out reads for tmin, tmax, precip, sand,
c-----silt, clay.
c-----C. Daly 5-21-97.  Also commenting out nlayer, nlaypg. 
c-----C. Daly 5-21-97 Bringing in which_soils, so reads will be made for
c-----sand, silt, clay, nlayer, and nlaypg if no gridded data available. 

      read(7,*)
      read(7,*)

      do 10 i = 1, MONTHS
        read(7,*) dum(i), name
cx        read(7,*) precip(i), name
cx        call ckdata('sitein','precip',name)
10    continue

      do 20 i = 1, MONTHS
        read(7,*) prcstd(i), name
        call ckdata('sitein','prcstd',name)
20    continue

      do 30 i = 1, MONTHS
        read(7,*) prcskw(i), name
        call ckdata('sitein','prcskw',name)
30    continue

      do 40 i = 1, MONTHS
        read(7,*) dum(i), name
cx        read(7,*) tmn2m(i), name
cx        call ckdata('sitein','tmn2m',name)
40    continue

      do 50 i = 1, MONTHS
        read(7,*) dum(i), name
cx        read(7,*) tmx2m(i), name
cx        call ckdata('sitein','tmx2m',name)
50    continue

      read(7,*)
      read(7,*) temp, name
      ivauto = int(temp)
      call ckdata('sitein','ivauto',name)
      read(7,*) temp, name
      nelem = int(temp)
      call ckdata('sitein','nelem',name)
      if (nelem .gt. 1 .or. nelem .lt. 1) then
        print *, "*** NELEM input value in vvegTypeNN.100 file is ", nelem 
        print *, " but MC2 expects nelem to be 1."
        print *, "MC2 will set nelem to 1 and continue."
        nelem = 1
      endif
 
      read(7,*) temp, name
cx      read(7,*) sitlat, name
cx      call ckdata('sitein','sitlat',name)
c...Make sure latitude is a positive number
      sitlat = abs(sitlat)
      read(7,*) sitlng, name
      call ckdata('sitein','sitlng',name)

      if (which_soils .lt. 1 .or. which_soils .gt. 2) then
        read(7,*) sand, name
        call ckdata('sitein','sand',name)
      else
        read(7,*) temp, name
      endif
      if (which_soils .lt. 1 .or. which_soils .gt. 2) then
        read(7,*) silt, name
        call ckdata('sitein','silt',name)
      else
        read(7,*) temp, name
      endif
      if (which_soils .lt. 1 .or. which_soils .gt. 2) then
        read(7,*) clay, name
        call ckdata('sitein','clay',name)
      else
        read(7,*) temp, name
      endif

      read(7,*) temp, name
cx      read(7,*) bulkd, name
cx      call ckdata('sitein','bulkd',name)

      if (code_flags(10) 
     $     .or. which_soils .lt. 1 .or. which_soils .gt. 2) then
        read(7,*) temp, name
        nlayer = int(temp)
        call ckdata('sitein','nlayer',name)
        if (nlayer .gt. 10) nlayer = 10
        if (code_flags(1) .and. nlayer .gt. 9) nlayer = 9
      else
        read(7,*) temp, name
      endif

      read(7,*) temp, name
      nlaypg = int(temp)
      call ckdata('sitein','nlaypg',name)

c-----C. Daly 5-21-97.  If nlaypg is greater than nlayer, set to nlayer
      if (nlaypg .gt. nlayer) nlaypg = nlayer

      read(7,*) drain, name
      call ckdata('sitein','drain',name)
      read(7,*) basef, name
      call ckdata('sitein','basef', name)
      read(7,*) stormf, name
      call ckdata('sitein','stormf', name)
      read(7,*) temp, name
      swflag = int(temp)
      call ckdata('sitein','swflag', name)

      do 60 i = 1, MAXLYR
        read(7,*) awilt(i), name
        call ckdata('sitein','awilt', name)
60    continue

      do 70 i = 1, MAXLYR
        read(7,*) afiel(i), name
        call ckdata('sitein','afiel', name)
70    continue

      read(7,*) ph, name
      call ckdata('sitein','ph',name)
      read(7,*) pslsrb, name
      call ckdata('sitein','pslsrb',name)
      read(7,*) sorpmx, name
      call ckdata('sitein','sorpmx',name)
      read(7,*)
      read(7,*) epnfa(INTCPT), name
      call ckdata('sitein','epnfa',name)
      read(7,*) epnfa(SLOPE), name
      call ckdata('sitein','epnfa',name)
      read(7,*) epnfs(INTCPT), name
      call ckdata('sitein','epnfs',name)
      read(7,*) epnfs(SLOPE), name
      call ckdata('sitein','epnfs',name)
      read(7,*) satmos(INTCPT), name
      call ckdata('sitein','satmos',name)
      read(7,*) satmos(SLOPE), name
      call ckdata('sitein','satmos',name)
      read(7,*) sirri, name
      call ckdata('sitein','sirri',name)

c...If extending, do not read in initial conditions
      if (ivopt) then
        return
      endif

      read(7,*)
      read(7,*) som1ci(SRFC,UNLABL), name
      call ckdata('sitein','som1ci',name)
      read(7,*) som1ci(SRFC,LABELD), name
      call ckdata('sitein','som1ci',name)
      read(7,*) som1ci(SOIL,UNLABL), name
      call ckdata('sitein','som1ci',name)
      read(7,*) som1ci(SOIL,LABELD), name
      call ckdata('sitein','som1ci',name)

      read(7,*) som2ci(UNLABL), name
      call ckdata('sitein','som2ci',name)
      read(7,*) som2ci(LABELD), name
      call ckdata('sitein','som2ci',name)

      read(7,*) som3ci(UNLABL), name
      call ckdata('sitein','som3ci',name)
      read(7,*) som3ci(LABELD), name
      call ckdata('sitein','som3ci',name)

      do 90 i = SRFC, SOIL
        do 80 j = 1, MAXIEL
          read(7,*) rces1(i,j), name
          call ckdata('sitein','rces1',name)
80      continue
90    continue

      do 100 i = 1, MAXIEL
        read(7,*) rces2(i), name
        call ckdata('sitein','rces2',name)
100   continue

      do 110 i = 1, MAXIEL
        read(7,*) rces3(i), name
        call ckdata('sitein','rces3',name)
110   continue

      read(7,*) clittr(SRFC,UNLABL), name
      call ckdata('sitein','clittr',name)
      read(7,*) clittr(SRFC,LABELD), name
      call ckdata('sitein','clittr',name)
      read(7,*) clittr(SOIL,UNLABL), name
      call ckdata('sitein','clittr',name)
      read(7,*) clittr(SOIL,LABELD), name
      call ckdata('sitein','clittr',name)

      do 130 i = SRFC, SOIL
        do 120 j = 1, MAXIEL
          read(7,*) rcelit(i,j), name
          call ckdata('sitein','rcelit',name)
120     continue
130   continue

      read(7,*) aglcis(UNLABL), name
      call ckdata('sitein','aglcis',name)
      read(7,*) aglcis(LABELD), name
      call ckdata('sitein','aglcis',name)

      do 140 i = 1, MAXIEL
        read(7,*) aglive(i), name
        call ckdata('sitein','aglive',name)
140   continue

      read(7,*) bglcis(UNLABL), name
      call ckdata('sitein','bglcis',name)
      read(7,*) bglcis(LABELD), name
      call ckdata('sitein','bglcis',name)

      do 150 i = 1, MAXIEL
        read(7,*) bglive(i), name
        call ckdata('sitein','bglive',name)
150   continue

      read(7,*) stdcis(UNLABL), name
      call ckdata('sitein','stdcis',name)
      read(7,*) stdcis(LABELD), name
      call ckdata('sitein','stdcis',name)

      do 160 i = 1, MAXIEL
        read(7,*) stdede(i), name
        call ckdata('sitein','stdede',name)
160   continue

      read(7,*)
      read(7,*) rlvcis(UNLABL), name
      call ckdata('sitein','rlvcis',name)
      read(7,*) rlvcis(LABELD), name
      call ckdata('sitein','rlvcis',name)

      do 170 i = 1, MAXIEL
        read(7,*) rleave(i), name
        call ckdata('sitein','rleave',name)
170   continue

      read(7,*) fbrcis(UNLABL), name
      call ckdata('sitein','fbrcis',name)
      read(7,*) fbrcis(LABELD), name
        call ckdata('sitein','fbrcis',name)

      do 180 i = 1, MAXIEL
        read(7,*) fbrche(i), name
      call ckdata('sitein','fbrche',name)
180   continue

      read(7,*) rlwcis(UNLABL), name
      call ckdata('sitein','rlwcis',name)
      read(7,*) rlwcis(LABELD), name
      call ckdata('sitein','rlwcis',name)

      do 190 i = 1, MAXIEL
        read(7,*) rlwode(i), name
        call ckdata('sitein','rlwode',name)
190   continue

      read(7,*) frtcis(UNLABL), name
      call ckdata('sitein','frtcis',name)
      read(7,*) frtcis(LABELD), name
      call ckdata('sitein','frtcis',name)

      do 200 i = 1, MAXIEL
        read(7,*) froote(i), name
        call ckdata('sitein','froote',name)
200   continue

      read(7,*) crtcis(UNLABL), name
      call ckdata('sitein','crtcis',name)
      read(7,*) crtcis(LABELD), name
      call ckdata('sitein','crtcis',name)

      do 210 i = 1, MAXIEL
        read(7,*) croote(i), name
        call ckdata('sitein','croote',name)
210   continue

      read(7,*) wd1cis(UNLABL), name
      call ckdata('sitein','wd1cis',name)
      read(7,*) wd1cis(LABELD), name
      call ckdata('sitein','wd1cis',name)
      read(7,*) wd2cis(UNLABL), name
      call ckdata('sitein','wd2cis',name)
      read(7,*) wd2cis(LABELD), name
      call ckdata('sitein','wd2cis',name)
      read(7,*) wd3cis(UNLABL), name
      call ckdata('sitein','wd3cis',name)
      read(7,*) wd3cis(LABELD), name
      call ckdata('sitein','wd3cis',name)

      read(7,*) w1lig, name
      call ckdata('sitein','w1lig',name)
      read(7,*) w2lig, name
      call ckdata('sitein','w2lig',name)
      read(7,*) w3lig, name
      call ckdata('sitein','w3lig',name)
      read(7,*)

       do 230 i = 1, MAXIEL
        do 220 j = 1, MAXLYR
          read(7,*) minerl(j,i), name
          call ckdata('sitein','minerl',name)
220     continue
230   continue

      do 240 i = 1, MAXIEL
        read(7,*) parent(i), name
        call ckdata('sitein','parent',name)
240   continue

      do 250 i = 1, MAXIEL
        read(7,*) secndy(i), name
        call ckdata('sitein','secndy',name)
250   continue

      read(7,*) occlud, name
      call ckdata('sitein','occlud',name)
      read(7,*)

      do 260 i = 1, MAXLYR
        read(7,*) rwcf(i), name
        call ckdata('sitein','rwcf',name)
260   continue

      read(7,*) snlq, name
      call ckdata('sitein','snlq',name)
      read(7,*) snow, name
      call ckdata('sitein','snow', name)

      close(unit=7)

      return
      end
