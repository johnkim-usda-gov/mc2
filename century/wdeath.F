
c*****************************************************************
c	$RCSfile: wdeath.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: wdeath.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:53:45  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

c...WDEATH

      subroutine wdeath (tave, wfunc)

      real tave, wfunc

c...Death of leaves, fine branches, large wood, fine roots, and coarse roots.

      include 'const.inc'
      include 'dovars.inc'
      include 'isovar.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'pheno.inc'
      include 'plot3.inc'
      include 'timvar.inc'
      include 'zztim.inc'
      include 'forrem.inc'

      include 'mapssvar.inc'
      include 'vemap2.inc'

c...Local variables
      integer iel
      real accum(ISOS), ctodie, etodie, fr14, recres(MAXIEL), tostore

      real ctodie_decid, ctodie_ever, ctodie_ever1, ctodie_ever2

c______LENIHAN 5/23/97

      real lfdr, wddrfb, wddrcw

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0


      lfdr   = leafdr(month)
      wddrfb = wooddr(FBRCH)
      wddrcw = wooddr(LWOOD)
      
      if ( did_burn ) then

        if (mf_dr(1) .gt. 0.) lfdr   = lfdr + mf_dr(1)
        if (mf_dr(2) .gt. 0.) wddrfb = wddrfb + mf_dr(2)
        if (mf_dr(3) .gt. 0.) wddrcw = wddrcw + mf_dr(3)
          
        if (lfdr .gt. 0.99) lfdr = 0.99
        if (wddrfb .gt. 0.99) wddrfb = 0.99
        if (wddrcw .gt. 0.99) wddrcw = 0.99

      endif

c_________________________________

c...Death of leaves
c   NOTE:  WOODDR(1) is the death rate in fall for deciduous forests
c          LEAFDR(MTH) is the monthly death rate for leaves in every
c                      case except for fall in deciduous forests.

      if (rleavc .gt. 0.0) then

c-------C. Daly 6-23-97.  Deciduousness is now a continuous function,
c-------determined by wooddr(1).  So, calculate leaf death in both
c-------evergreen and deciduous modes, and take a proportion of each.

c-------C. Daly 2-12-98.  Calculating drought-deciduousness for tropical
c-------and subtropical zones only.  This leaf drop rate is compared
c-------to the monthly turnover rate and the max taken as ctodie for
c-------evergreen.  Evergreen is used here because we want this to
c-------apply to tropical lifeforms, which are termed evergreen
c-------broadleaf by the mix_index.

c.......Evergreen trees

c       Use leaf death rate multiplier from EACHYR
        ctodie_ever1 = rleavc * lfdr * ldrmlt

c.......Drought deciduousness for evergreen broadleaf trees, for tropical
c-------and subtropical zones only. BOREAL=0,TEMPERATE=1,SUBTROPICAL=2,
c-------TROPICAL=3.
c...    changed this to only include tropical zone 4-27-98 DB
c...        if (zone .eq. 2 .or. zone .eq. 3) then 

        if (zone .eq. 3) then
	  ctodie_ever2 = rleavc * (1. - wfunc)
        else
          ctodie_ever2 = 0.
        endif

        if (rleavc .gt. 0.) then
          ddecid = ctodie_ever2 / rleavc
        else
          ddecid = 0.
        endif

c-------Take max of monthly death rate or drought-decid.
        ctodie_ever = MAX(ctodie_ever1, ctodie_ever2)

c.......Deciduous trees

c.......Cold deciduous trees JL 11-18-98
	  if (tave .lt. 13.0) then
	    ctodie_decid = rleavc
          else
	    ctodie_decid = rleavc * lfdr
          endif

c-------Apportion death rates from decid and evergreen components.
        ctodie = (wooddr(LEAF) * ctodie_decid) + 
     .           (1.-wooddr(LEAF))*ctodie_ever

c-------LENIHAN 9/15/98: Assign ctodie to leaf litterfall var 
c-------for VEMAP2 output

        lleaf = ctodie

c...Compute E/C ratios
        do 10 iel = 1, nelem
	 recres(iel) = rleave(iel) / rleavc

c...Compute flow to retranslocation storage
	  tostore = recres(iel) * ctodie * forrtf(iel)
c	  print *, 'wdeath: tostore = ', tostore
	  call flow(rleave(iel), forstg(iel), time, tostore)

c...Decrease E/C by the amount that is retranslocated
	  recres(iel) = recres(iel) * (1 - forrtf(iel))
 10     continue

	fr14 = rlvcis(LABELD) /rleavc
	call partit(ctodie, recres, 1, rlvcis, rleave, wdlig(LEAF),
     $   fr14)
      endif

c...Death of fine roots
c... Bachelet 11-5-97
      	if (frootc .ne. 0.0) then
		ctodie = frootc * wooddr(FROOT)
	else
		ctodie = 0.0
	endif

c-------LENIHAN 9/15/98: Assign ctodie to fine root litterfall var 
c-------for VEMAP2 output

        lfroot = ctodie

	do 20 iel = 1, nelem
	  recres(iel) = froote(iel) / frootc
 20     continue
	if (frootc .ne. 0.0) then
		fr14 = frtcis(LABELD) / frootc
	else
                fr14 = 0.0 
        endif 
	call partit(ctodie, recres, 2, frtcis, froote, wdlig(FROOT),
     -    fr14)

c...Fine Branches, Large Wood, and Coarse Roots go to the dead wood
c   compartments: WOOD1, WOOD2, WOOD3

c...Death of fine branches
      if (fbrchc .gt. 0.0) then
	ctodie = fbrchc * wddrfb

c-------LENIHAN 9/15/98: Assign ctodie to fine branch litterfall var 
c-------for VEMAP2 output

        lfwood = ctodie

	call csched(ctodie, cisotf, 1.0,
     $              fbrcis(UNLABL), wd1cis(UNLABL),
     $              fbrcis(LABELD), wd1cis(LABELD),
     $              1.0, accum)

	do 30 iel = 1, nelem
c... Bachelet 11-5-97
	  if(fbrchc .ne. 0.0) then
	  	etodie = ctodie * (fbrche(iel) / fbrchc)
	  else
		etodie = 0.0
	  endif
	  call flow(fbrche(iel), wood1e(iel), time, etodie)
 30     continue
      endif

c...Death of large wood
      if (rlwodc .gt. 0.0) then
	ctodie = rlwodc * wddrcw

c-------LENIHAN 9/15/98: Assign ctodie to large wood litterfall var 
c-------for VEMAP2 output

        llwood = ctodie

	call csched(ctodie, cisotf, 1.0,
     $              rlwcis(UNLABL), wd2cis(UNLABL),
     $              rlwcis(LABELD), wd2cis(LABELD),
     $              1.0, accum)

	do 40 iel = 1, nelem
c... Bachelet 11-5-97
	  if(rlwodc .ne. 0.0) then
	  	etodie = ctodie * (rlwode(iel) / rlwodc)
	  else
	  	etodie =0.0
	  endif
	  call flow(rlwode(iel), wood2e(iel), time, etodie)
 40     continue
      endif

c...Death of coarse roots
c... Bachelet 11-5-97
      if (crootc .gt. 0.0) then
	ctodie = crootc * wooddr(CROOT)

c-------LENIHAN 9/15/98: Assign ctodie to coarse root litterfall var 
c-------for VEMAP2 output

        lcroot = ctodie

	call csched(ctodie, cisotf, 1.0,
     $              crtcis(UNLABL), wd3cis(UNLABL),
     $              crtcis(LABELD), wd3cis(LABELD),
     $              1.0, accum)

	do 50 iel = 1, nelem
	  if(crootc .ne. 0.0) then
		etodie = ctodie * (croote(iel) / crootc)
	  else
	  	etodie = 0.0
	  endif
	  call flow(croote(iel), wood3e(iel), time, etodie)
 50     continue
      endif

      return
      end
