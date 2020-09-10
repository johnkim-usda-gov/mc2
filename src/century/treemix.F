
c*****************************************************************
c	$RCSfile: treemix.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: treemix.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:49:32  bachelet
cheader
c
c******************************************************************

      subroutine treemix      

c-----C. Daly 11-22-96
c-----Index "mix" describes mix between subtropical
c-----evergreen broadleaf, temperate deciduous, and boreal 
c-----evergreen needleleaf for eastern US. This will be done
c-----for deciduous tree types (CWT0 only. Index
c-----ranges from 100 (all evergreen broadleaf) to 0 (all
c-----deciduous broadleaf) to -100 (all evergreen needleleaf).
c-----CWT is the only temperate deciduous tree type in Century.
c-----CSU is providing us with some mixed type tree.100 sections,
c-----so this is how things will look:
c
c-----If tree = CWT only.
c-----Keep same schedule and 100 file, but change tree.100 section 
c-----as follows:

c-----C. Daly
c-----Now allowing WPINE into mix calcs.  Changed 75WP threshold from
c----- -62.5 to -50.0.

c-----C. Daly 7-30-97
c-----New mix index calculated by Jim's routine MixIndex.
c-----Now, index varies from 0 (pure deciduous) to 100 (pure evergreen).
c-----Using "supertree" tree type, so instead of changing tree type
c-----based on this index, characteristics of the tree are changed.
c-----Here is present list of parameters changed:

c-----prdx(3) - max production
c-----ppdf(1-4) - temperature curve parameters
c-----cerfor(1,1,1), (1,2,1), (1,3,1), (1,4,1), (1,5,1) - min C/N ratios
c-----cerfor(2,1,1), (2,2,1), (2,3,1), (2,4,1), (2,5,1) - max C/N ratios
c-----fcfrac(1,2),(2,2),(3,2),(4,2),(5,2) - C allocation
c-----leafdr(1-12) - monthly leaf death rate
c-----btolai - SLA
c-----maxlai - max potential LAI
c-----wooddr(1-5) - monthly death rates
c-----sitpot - N given to grasses

c-----Note that decid is always set to 0, unless drought-deciduous, 
c-----then set to 2.

c-----C. Daly - 10-29-97 Added klai to mixed index scaling.

c-----C. Daly - 1-23-98.  Expanded routine to handle 4 columns of
c-----parameters from tree.100, not just 2.  Columns in tree.100,
c-----with respective values of a reformulated mx_index, are:
c----- col 1  -200  Deciduous needleleaf (larch type)
c----- col 2  -100  Evergreen needleleaf (evergreen conifers)
c----- col 3     0  Deciduous broadleaf  (temperate deciduous trees)
c----- col 4   100  Evergreen broadleaf  (tropical evergreen trees)

c-----C. Daly 2-12-98.  Changing routine again to make use of two
c-----mix indexes, one for ppt (ppt_index) and one for temp (tmp_index).
c-----Now, there will be two gradients of tree types: one for dry
c-----areas and one for wet areas.  The wet gradient is DN-EN-DB-EB,
c-----while the dry gradient replaces DB with EN: DN-EN-EN-EB.  If cell 
c-----is in a transition between dry and wet, parameters from both dry 
c-----and wet are included in a weighted average. 
 
c-----Values of the two indexes are as follows:

c-----ppt_index   0=wet 100=dry 1-99 transition
c-----tmp_index   -200=DN -100=EN 100=DB

      include 'chrvar.inc'
      include 'const.inc'
      include 'isovar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'site.inc'
      include 'plot4.inc'

      include 'mapssvar.inc'
      include 'supertree.inc'
      include 'vemap2.inc'

c-----------------------------------------------------------------------
      integer switch
c-----If treemix is to be switched off, set to 0; set to 1 to switch on.
c-----If switched off, default values of supertree will be used.

      data switch /1/  

c-----C. Daly 1-23-98.  Set mx_index thresholds to link to tree.100 cols.
      real thresh(4)
      data  thresh / -200., -100., 0., 100. / 

c-----C. Daly 7/16/98.  Setting columns of a needleleaf value that will
c-----enable us to get a mx_index from the method below.  See needleleaf
c-----code.  Mx_index of 1 indicates 100% needleleaf.

      real fac, chkfrc, ppti
      integer i, ic, j, k
c-----------------------------------------------------------------------

      if (switch.eq.1) then

              
        if (diags) write(*,*) 'Based on tmp_index of ',tmp_index,
     .      ' and ppt_index of ',ppt_index,  
     .      '  some tree.100 params changed. '


c-------C. Daly 1-23-98. Determine the first of the two columns parameters
c-------will be interpolated from.  Then scale mx_index to 0-1.
c-------C. Daly 2-17-98.  Now using tmp_index for this.
  
        if (tmp_index .lt. thresh(2) ) then
          ic = 1
        else
          if (tmp_index .ge. thresh(2) .and.
     .       tmp_index .lt. thresh(3)) then
            ic = 2
          else
            ic = 3
          endif
        endif
cx        write(*,*) 'ic = ', ic

c-------Get mix index factor from original -200 - 100 range.
c-------Now it will range from 0-1, with ic indicating which
c-------tree.100 columns to interpolate between.
c-------C. Daly 2-17-98.  Now using tmp_index for this.

        fac = (tmp_index-thresh(ic))/abs(thresh(ic+1)-thresh(ic))
        if (fac .gt. 1.0) fac = 1.0
        if (fac .lt. 0.0) fac = 0.0

c-------C. Daly 2-17-98.  If ppt_index=0, just do full column interpolation
c-------as before.  If ppt_index=100, replace column 3 value (DB) with
c-------column 2 value (EN).  If 0<ppt_index<100, do both and perform
c-------a weighted average.

c-------WET CASE----------------

c-------Max production
        do 30 i = 2, 3
          call scale4 (fac, ic, prdx1(i), prdx2(i), 
     .                prdx3(i), prdx4(i), prdxw(i))
cx          write(*,*) 'prdxw(',i,') = ',prdxw(i)
30      continue

c-------Temp curves
        do 40 i = 1, 4
          call scale4 (fac, ic, ppdf1(i,2), ppdf2(i,2), 
     .                 ppdf3(i,2), ppdf4(i,2), ppdfw(i,2))
cx          write(*,*) 'ppdfw(',i,', 2) = ',ppdfw(i,2)
40      continue

c-------C/N ratios
        do 70 i = IMIN, IVAL
          do 60 j = 1, FPARTS
            do 50 k = 1, MAXIEL
              call scale4 (fac, ic, cerfor1(i,j,k),
     .        cerfor2(i,j,k), cerfor3(i,j,k), cerfor4(i,j,k),
     .        cerforw(i,j,k))
50          continue
60        continue
70      continue

c-------Carbon allocation
        do 90 j = NEWFOR, OLDFOR
	  chkfrc = 0.0
          do 80 i = 1, FPARTS
            call scale4 (fac, ic, fcfrac1(i, j),fcfrac2(i, j),
     .                  fcfrac3(i, j), fcfrac4(i, j), fcfracw(i, j))
	    chkfrc = chkfrc + fcfracw(i, j)
cx            write(*,*) 'fcfracw(',i,',',j,') = ',fcfracw(i,j)
80	  continue
          if (abs(1.0 - chkfrc) .gt. 0.000001) then
	    call message('Total of all mixed FCFRAC != 1')
            STOP
          endif
90      continue

c-------Monthly leaf death rate
        do 100 i = 1, MONTHS
          call scale4(fac, ic, leafdr1(i), leafdr2(i), 
     .                leafdr3(i), leafdr4(i), leafdrw(i))
cx          write(*,*) 'leafdrw(',i,') = ',leafdrw(i)
100     continue

c-------SLA
        call scale4(fac, ic, btolai1, btolai2, 
     .              btolai3, btolai4, btolaiw)
cx        write(*,*) 'btolaiw = ', btolaiw

c-------klai - 
        call scale4(fac, ic, klai1, klai2, 
     .              klai3, klai4, klaiw)
cx        write(*,*) 'klaiw = ', klaiw

c-------maxlai - max potential LAI
        call scale4(fac, ic, maxlai1, maxlai2, 
     .              maxlai3, maxlai4, maxlaiw)
cx        write(*,*) 'maxlaiw = ', maxlaiw

c-------wooddr(1-5) - monthly death rates
        do 140 i = 1, FPARTS
          call scale4( fac, ic, wooddr1(i), wooddr2(i), 
     .                 wooddr3(i), wooddr4(i), wooddrw(i))
cx          write(*,*) 'wooddrw(',i,') = ',wooddrw(i)
140     continue

c-------sitpot - N given to grasses
        call scale4(fac, ic, sitpot1, sitpot2, 
     .              sitpot3, sitpot4, sitpotw)
cx        write(*,*) 'sitpotw = ', sitpotw

c-------Needleleaf index
        call scale4(fac, ic, 100., 100., 0., 0., needlew)
c          write(*,*) 'needlew = ', needlew

c-------Evergreen index
        call scale4(fac, ic, 0., 100., 0., 100., evergw)



c-------DRY CASE----------------

c-------Max production
        do 230 i = 2, 3
          call scale4 (fac, ic, prdx1(i), prdx2(i), 
     .                prdx2(i), prdx4(i), prdxd(i))
cx          write(*,*) 'prdxd(',i,') = ',prdxd(i)
230      continue

c-------Temp curves
        do 240 i = 1, 4
          call scale4 (fac, ic, ppdf1(i,2), ppdf2(i,2), 
     .                 ppdf2(i,2), ppdf4(i,2), ppdfd(i,2))
cx          write(*,*) 'ppdfd(',i,', 2) = ',ppdfd(i,2)
240      continue

c-------C/N ratios
        do 270 i = IMIN, IVAL
          do 260 j = 1, FPARTS
            do 250 k = 1, MAXIEL
              call scale4 (fac, ic, cerfor1(i,j,k), cerfor2(i,j,k),
     .        cerfor2(i,j,k), cerfor4(i,j,k), cerford(i,j,k))
250          continue
260        continue
270      continue

c-------Carbon allocation
        do 290 j = NEWFOR, OLDFOR
	  chkfrc = 0.0
          do 280 i = 1, FPARTS
            call scale4 (fac, ic, fcfrac1(i, j),fcfrac2(i, j),
     .                  fcfrac2(i, j), fcfrac4(i, j), fcfracd(i, j))
	    chkfrc = chkfrc + fcfracd(i, j)
cx            write(*,*) 'fcfracd(',i,',',j,') = ',fcfracd(i,j)
280	  continue
          if (abs(1.0 - chkfrc) .gt. 0.000001) then
	    call message('Total of all mixed FCFRAC != 1')
            STOP
          endif
290      continue

c-------Monthly leaf death rate
        do 300 i = 1, MONTHS
          call scale4(fac, ic, leafdr1(i), leafdr2(i), 
     .                leafdr2(i), leafdr4(i), leafdrd(i))
cx          write(*,*) 'leafdrd(',i,') = ',leafdrd(i)
300     continue

c-------SLA
        call scale4(fac, ic, btolai1, btolai2, 
     .              btolai2, btolai4, btolaid)
cx        write(*,*) 'btolaid = ', btolaid

c-------klai - 
        call scale4(fac, ic, klai1, klai2, 
     .              klai2, klai4, klaid)
cx        write(*,*) 'klaid = ', klaid

c-------maxlai - max potential LAI
        call scale4(fac, ic, maxlai1, maxlai2, 
     .              maxlai2, maxlai4, maxlaid)
cx        write(*,*) 'maxlaid = ', maxlaid

c-------wooddr(1-5) - monthly death rates
        do 340 i = 1, FPARTS
          call scale4( fac, ic, wooddr1(i), wooddr2(i), 
     .                 wooddr2(i), wooddr4(i), wooddrd(i))
cx          write(*,*) 'wooddrd(',i,') = ',wooddrd(i)
340     continue

c-------sitpot - N given to grasses
        call scale4(fac, ic, sitpot1, sitpot2, 
     .              sitpot2, sitpot4, sitpotd)
cx        write(*,*) 'sitpotd = ', sitpotd

c-------Needleleaf index
        call scale4(fac, ic, 100., 100., 100., 0., needled)
cx          write(*,*) 'needled = ', needled

c-------Evergreen index
        call scale4(fac, ic, 0., 100., 100., 100., evergd)



c-------WEIGHTED AVERAGE OF WET AND DRY GRADIENTS

        ppti = ppt_index * 0.01

c-------Max production
        do 430 i = 2, 3
          prdx(i) = ppti*prdxd(i) + (1.-ppti)*prdxw(i)
cx          write(*,*) 'prdx(',i,') = ',prdx(i)
430      continue

c-------Temp curves
        do 440 i = 1, 4
          ppdf(i,2) = ppti*ppdfd(i,2) + (1.-ppti)*ppdfw(i,2)
cx          write(*,*) 'ppdf(',i,', 2) = ',ppdf(i,2)
440      continue

c-------C/N ratios
        do 470 i = IMIN, IVAL
          do 460 j = 1, FPARTS
            do 450 k = 1, MAXIEL
              cerfor(i,j,k) = ppti*cerford(i,j,k) + 
     .       (1.-ppti)*cerforw(i,j,k)
450          continue
460        continue
470      continue

c-------Carbon allocation
        do 490 j = NEWFOR, OLDFOR
	  chkfrc = 0.0
          do 480 i = 1, FPARTS
	    chkfrc = chkfrc + fcfrac(i, j)
            fcfrac(i,j) = ppti*fcfracd(i,j) + (1.-ppti)*fcfracw(i,j)
cx            write(*,*) 'fcfrac(',i,',',j,') = ',fcfrac(i,j)
480	  continue
          if (abs(1.0 - chkfrc) .gt. 0.000001) then
	    call message('Total of all mixed FCFRAC != 1')
            STOP
          endif
490      continue

c-------Monthly leaf death rate
        do 500 i = 1, MONTHS
          leafdr(i) = ppti*leafdrd(i) + (1.-ppti)*leafdrw(i)
cx          write(*,*) 'leafdr(',i,') = ',leafdr(i)
500     continue

c-------SLA
        btolai = ppti*btolaid + (1.-ppti)*btolaiw
cx        write(*,*) 'btolai = ', btolai

c-------klai - 
        klai = ppti*klaid + (1.-ppti)*klaiw
cx        write(*,*) 'klai = ', klai

c-------maxlai - max potential LAI
        maxlai = ppti*maxlaid + (1.-ppti)*maxlaiw
cx        write(*,*) 'maxlai = ', maxlai

c-------wooddr(1-5) - monthly death rates
        do 540 i = 1, FPARTS
          wooddr(i) = ppti*wooddrd(i) + (1.-ppti)*wooddrw(i)
cx          write(*,*) 'wooddr(',i,') = ',wooddr(i)
540     continue

c-------sitpot - N given to grasses
        sitpot = ppti*sitpotd + (1.-ppti)*sitpotw
cx        write(*,*) 'sitpot = ', sitpot

c-------Needleleaf index
        needle_index = ppti*needled + (1.-ppti)*needlew
        mx_index = needle_index
c        write(*,*) 'mx_index = ', mx_index

c-------Evergreen index
        everg_index = ppti*evergd + (1.-ppti)*evergw


      endif
      
      return
      end
