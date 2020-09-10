
c*****************************************************************
c	$RCSfile: cropin.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: cropin.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:03:09  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c-----C. Daly 8-1-97.  Modified function to read in C3/C4
c-----values of several parameters from crop.100 file.
c-----These values are used in cropmix to determine values for
c-----each cell based on c3c4_ratio.   These parameters
c-----are read in as default, C3, C4.
 
      subroutine cropin(tomatch)

      character*5 tomatch

c...Read in the new crop type

      include 'chrvar.inc'
      include 'comput.inc'
      include 'const.inc'
      include 'isovar.inc'
      include 'ligvar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'seq.inc'

      include 'mapssvar.inc'
      include 'supergrass.inc'
      include 'supertree.inc'

c...Local variables
      integer   i, j, layer
      real      del13c, temp
      character name*6, string*80

      call reset_crop_params()
      
c     skip over header
      call get_crop1(temp, name)
      
c-----Read in as default, C3, C4
      call get_crop3( prdx(1), prdx1(1), prdx2(1), name)
      call ckdata('cropin','prdx',name)

c-----Read in as default, C3, C4
      do 30 i = 1, 4
        call get_crop3( ppdf(i,1), ppdf1(i,1), ppdf2(i,1), name)
        call ckdata('cropin','ppdf',name)
30    continue

      call get_crop1( temp, name)
      bioflg = int(temp)
      call ckdata('cropin','bioflg',name)
      call get_crop1( biok5, name)
      call ckdata('cropin','biok5',name)
      call get_crop1( pltmrf, name)
      call ckdata('cropin','pltmrf',name)
      call get_crop1( fulcan, name)
      call ckdata('cropin','fulcan',name)

      do 40 i = 1, 3
        call get_crop1( frtc(i), name)
        call ckdata('cropin','frtc',name)
40    continue
      if (frtc(3) .eq. 0) frtc(3) = 1

      call get_crop1( biomax, name)
      call ckdata('cropin','biomax',name)

c-----Read in as default, C3, C4
      do 60 i = 1, 2
        do 50 j = 1, MAXIEL
          call get_crop3( pramn(j,i), pramn1(j,i), pramn2(j,i), name)
          call ckdata('cropin','pramn',name)
50      continue
60    continue

c-----Read in as default, C3, C4
      do 80 i = 1, 2
        do 70 j = 1, MAXIEL
          call get_crop3( pramx(j,i), pramx1(j,i), pramx2(j,i), name)
          call ckdata('cropin','pramx',name)
70      continue
80    continue

c-----Read in as default, C3, C4
      do 100 i = 1, 2
        do 90 j = 1, MAXIEL
          call get_crop3( prbmn(j,i), prbmn1(j,i), prbmn2(j,i), name)
          call ckdata('cropin','prbmn',name)
90      continue
100   continue

c-----Read in as default, C3, C4
      do 120 i = 1, 2
        do 110 j = 1, MAXIEL
          call get_crop3( prbmx(j,i), prbmx1(j,i), prbmx2(j,i), name)
          call ckdata('cropin','prbmx',name)
110     continue
120   continue

      do 140 i = ABOVE, BELOW
        do 130 j = INTCPT, SLOPE
          call get_crop1( fligni(j,i), name)
          call ckdata('cropin','fligni',name)
130     continue
140   continue

      call get_crop1( himax, name)
      call ckdata('cropin','himax',name)
      call get_crop1( hiwsf, name)
      call ckdata('cropin','hiwsf',name)
      call get_crop1( temp, name)
      himon(1) = int(temp)
      call ckdata('cropin','himon',name)
      call get_crop1(temp, name)
      himon(2) = int(temp)
      call ckdata('cropin','himon',name)

      do 150 i = 1, MAXIEL
        call get_crop1( efrgrn(i), name)
        call ckdata('cropin','efrgrn',name)
150   continue

      call get_crop1( vlossp, name)
      call ckdata('cropin','vlossp',name)

      do 160 i = 1, 4
        call get_crop1( fsdeth(i), name)
        call ckdata('cropin','fsdeth',name)
160   continue

      call get_crop1( fallrt, name)
      call ckdata('cropin','fallrt',name)
      call get_crop1( rdr, name)
      call ckdata('cropin','rdr',name)
      call get_crop1( rtdtmp, name)
      call ckdata('cropin','rtdtmp',name)

      do 170 i = 1, MAXIEL
        call get_crop1( crprtf(i), name)
        call ckdata('cropin','crprtf',name)
170   continue

      call get_crop1( snfxmx(CRPSYS), name)
      call ckdata('cropin','snfxmx',name)
      call get_crop1( del13c, name)
      call ckdata('cropin','del13c',name)
      call get_crop1( co2ipr(CRPSYS), name)
      call ckdata('cropin','co2ipr',name)
      call get_crop1( co2itr(CRPSYS), name)
      call ckdata('cropin','co2itr',name)

      do 190 i = IMIN, IMAX
        do 180 j = 1, MAXIEL
          call get_crop1( co2ice(CRPSYS,i,j), name)
          call ckdata('cropin','co2ice',name)
180     continue
190   continue

      call get_crop1( co2irs(CRPSYS), name)
      call ckdata('cropin','co2irs',name)

c     check for end of data
      call get_crop1(temp, name)
      call ckdata('cropin','end',name)

c...Hold on to the current crop just read in
      curcrp = tomatch

c...Determine the 'numerical value' of the curcrp, 
c     for use as an output variable
      crpval = 0
      do 200 i = 1, 5
        if (curcrp(i:i) .ne. ' ') then
	  if (curcrp(i:i) .ge. '0' .and. curcrp(i:i) .le. '9') then
            crpval = crpval +
     -               ((ichar(curcrp(i:i)) - ichar('0')) / 10.0)
          else
            crpval = crpval + (ichar(curcrp(i:i)) - ichar('A')) + 1
          endif
        endif
200   continue

c...Calculate available water holding capacity 11/91 lh
      awhc = 0.0
      do 210 layer = 1, nlaypg
        awhc = awhc + (afiel(layer) - awilt(layer)) * adep(layer)
210   continue

c...Recalculate lignin
      call cmplig(cursys,fligni,wdlig,pltlig)

c...Calculate cisofr as 13C if 13C labeling
      if (labtyp .eq. 2) then
         stop 'Unexpected code for MC1'
        cisofr = del13c * PEEDEE * 1.0e-03 + PEEDEE
        cisofr = 1 / (1/cisofr + 1)
      endif

c-----Caculate c3/c4 mixed values of some parameters.
      call cropmix


      return

      end
