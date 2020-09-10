
c*****************************************************************
c	$RCSfile: treein.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: treein.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:49:10  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

c-----C. Daly 12-4-96.  Added call to mix(), which determines
c-----if tree type should change to a mixed type.

c-----C. Daly 7-31-97.  Modified function to read in evergreen/
c-----deciduous values of several parameters from tree.100 file.
c-----These values are used in treemix to determine values for
c-----each cell based on the mx_index value.   These parameters
c-----are read in as default, evergreen, deciduous. 

c-----C. Daly - 1-23-98.  Expanded routine to handle 4 columns of
c-----parameters from tree.100, not just 2.  Columns in tree.100,
c-----with respective values of a reformulated mx_index, are:
c----- col 1  -200  Deciduous needleleaf (larch type)
c----- col 2  -100  Evergreen needleleaf (evergreen conifers)
c----- col 3     0  Deciduous broadleaf  (temperate deciduous trees)
c----- col 4   100  Evergreen broadleaf  (tropical evergreen trees)


C...TREEIN.F

      subroutine treein(tomatch)

      character tomatch*5

c...Read in the new forest type
 
      include 'chrvar.inc'
      include 'const.inc'
      include 'isovar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'site.inc'

      include 'mapssvar.inc'
      include 'supertree.inc'

c...Local variables
      integer   i, j, k
      real      del13c, temp, chkfrc, chkfrc1, chkfrc2
      character fromdat*5, name*20, string*80
      
      name = '                    '

      call reset_tree_params()

c     skip over header
      call get_tree1(temp, name)
      
      call get_tree1(temp, name)
      decid = int(temp)
      call ckdata('treein','decid',name)

c-----Read in as default, DN, EN, DB, EB
      do 30 i = 2, 3
        call get_tree5( prdx(i), prdx1(i), prdx2(i), 
     .              prdx3(i), prdx4(i), name)        
        call ckdata('treein','prdx',name)
30    continue

c-----Read in as default, DN, EN, DB, EB
      do 40 i = 1, 4
        call get_tree5( ppdf(i,2), ppdf1(i,2), ppdf2(i,2), 
     .              ppdf3(i,2), ppdf4(i,2), name)
        call ckdata('treein','ppdf',name)
40    continue

c-----Read in as default, DN, EN, DB, EB
      do 70 i = IMIN, IVAL
        do 60 j = 1, FPARTS
          do 50 k = 1, MAXIEL
            call get_tree5(cerfor(i,j,k),cerfor1(i,j,k),cerfor2(i,j,k),
     .                 cerfor3(i,j,k), cerfor4(i,j,k), name)
            call ckdata('treein','cerfor',name)
50        continue
60      continue
70    continue

      call get_tree1( decw1, name)
      call ckdata('treein','decw1',name)
      call get_tree1( decw2, name)
      call ckdata('treein','decw2',name)
      call get_tree1( decw3, name)
      call ckdata('treein','decw3',name)

c-----Read in as default, DN, EN, DB, EB
      do 90 j = NEWFOR, OLDFOR
c...Added check to make sure fractions add to 1
	chkfrc = 0.0
	chkfrc1 = 0.0
	chkfrc2 = 0.0
        do 80 i = 1, FPARTS
          call get_tree5( fcfrac(i, j),fcfrac1(i, j),fcfrac2(i, j), 
     .                fcfrac3(i, j), fcfrac4(i, j), name)
	  chkfrc = chkfrc + fcfrac(i, j)
	  chkfrc1 = chkfrc1 + fcfrac1(i, j)
	  chkfrc2 = chkfrc2 + fcfrac2(i, j)
          call ckdata('treein','fcfrac',name)
80	continue

        if (abs(1.0 - chkfrc) .gt. 0.000001) then
	  call message('Total of all FCFRAC != 1')
          STOP
        endif
        if (abs(1.0 - chkfrc1) .gt. 0.000001) then
	  call message('Total of all FCFRAC1 != 1')
          STOP
        endif
        if (abs(1.0 - chkfrc2) .gt. 0.000001) then
	  call message('Total of all FCFRAC2 != 1')
          STOP
        endif
90    continue

c-----Read in as default, DN, EN, DB, EB
      do 100 i = 1, MONTHS
        call get_tree5( leafdr(i), leafdr1(i), leafdr2(i), 
     .              leafdr3(i), leafdr4(i), name)
        call ckdata('treein','leafdr',name)
100   continue

c-----Read in as default, DN, EN, DB, EB
      call get_tree5( btolai, btolai1, btolai2,  
     .           btolai3, btolai4, name)
c-----Read in as default, DN, EN, DB, EB
      call ckdata('treein','btolai',name)
      call get_tree5( klai,  klai1, klai2,  
     .            klai3, klai4, name)
      call ckdata('treein','klai',name)
      call get_tree1( laitop, name)
      call ckdata('treein','laitop',name)
c-----Read in as default, DN, EN, DB, EB
      call get_tree5( maxlai, maxlai1, maxlai2,  
     .            maxlai3, maxlai4, name)
      call ckdata('treein','maxlai',name)
      call get_tree1( maxldr, name)
      call ckdata('treein','maxldr',name)

      do 120 i = 1, MAXIEL
        call get_tree1( forrtf(i), name)
        call ckdata('treein','forrtf',name)
120   continue

      call get_tree1( sapk, name)
      call ckdata('treein','sapk',name)
      call get_tree1( swold, name)
      call ckdata('treein','swold',name)

      do 130 i = 1, FPARTS
        call get_tree1( wdlig(i), name)
        call ckdata('treein','wdlig',name)
130   continue

c-----Read in as default, DN, EN, DB, EB
      do 140 i = 1, FPARTS
        call get_tree5( wooddr(i), wooddr1(i), wooddr2(i), 
     .              wooddr3(i), wooddr4(i), name)
        call ckdata('treein','wooddr',name)
140   continue

      call get_tree1( snfxmx(FORSYS), name)
      call ckdata('treein','snfxmx',name)
      call get_tree1( del13c, name)
      call ckdata('treein','del13c',name)

      call get_tree1( co2ipr(FORSYS), name)
      call ckdata('cropin','co2ipr',name)
      call get_tree1( co2itr(FORSYS), name)
      call ckdata('cropin','co2itr',name)

      do 160 i = IMIN, IMAX
        do 150 j = 1, MAXIEL
          call get_tree1( co2ice(FORSYS,i,j), name)
          call ckdata('cropin','co2ice',name)
150     continue
160   continue

      call get_tree1( co2irs(FORSYS), name)
      call ckdata('cropin','co2irs',name)

      call get_tree1( basfc2, name)
      call ckdata('treein','basfc2',name)
      call get_tree1( basfct, name)
      call ckdata('treein','basfct',name)

c-----Read in as default, DN, EN, DB, EB
      call get_tree5( sitpot, sitpot1, sitpot2,  
     .            sitpot3, sitpot4, name)
      call ckdata('treein','sitpot',name)
 
c     check for end of data
      call get_tree1(temp, name)
      call ckdata('treein','end',name)

c...Hold on to the current tree just read in
      curtre = tomatch

c...Calculate cisotf as 13C if 13C labeling
      if (labtyp .eq. 2) then
         stop 'Unexpected code for MC1'
        cisotf = del13c * PEEDEE * 1.0e-03 + PEEDEE
        cisotf = 1 / (1/cisotf + 1)
      endif

c-----Caculate ever/decid mixed values of some parameters.
c      write (*,*) 'calling treemix in treein.F.'
      call treemix
c      write (*,*) 'called treemix in treein.F.'


      return

      end
