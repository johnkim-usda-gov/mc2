
c*****************************************************************
c	$RCSfile: potcrp.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: potcrp.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.3  2005/05/25 19:24:35  lenihan
cchanges for output of production scalars
c
cRevision 1.2  2004/09/24 23:23:38  bachelet
cheader
c
c******************************************************************
c              Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine potcrp (month,cancvr)
      integer month
      real cancvr

c...Compute monthly production potential based upon montly precip
c     and restrict potential production based upon method specified
c     by grzeff.

      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'plot4.inc'
      include 'potent.inc'
c...schvar.inc is included here because it has the code_flags
      include 'schvar.inc'
      include 'seq.inc'
      include 'wth.inc'
      include 'supertree.inc'
      include 'vemap2.inc'

c...Function declarations
      real     gpdf, pevap, pprdwc, pprdwc_old
      external gpdf, pevap, pprdwc, pprdwc_old

c...Local variables
      integer  iel,mode
      real     agprod, aisc, bgp, bgprod, bio, bioc, biof,
     $         bioprd, bop, ctemp, fracrc, potprd, pptprd, ratlc, rtsh, 
     $         sdlng, shdmod, subcan, temp1, temp2, temp3, tmns, tmxs,
     $         k, percent_cover,c3potprd,c4potprd
      real pprdwcA, pprdwcB, pprdwcC 

c...Compute shading modifier for savanna
      if (cursys .eq. SAVSYS) then
       if (cancvr .eq. 0) then
         aisc = 0.
c-------C. Daly 4-24-96.  Catch underflows.
       else if ( -.0035 * (rleavc*2.5)/cancvr .lt. -10.) then
         aisc = 0.0
       else
         aisc = 5 * exp(-.0035 * (rleavc*2.5)/cancvr)
       endif
       subcan = aisc/(aisc + 1.)
       shdmod = (1.0-cancvr) + (cancvr*subcan)
      else
       shdmod = 1.0
      endif

c.....LENIHAN 9/18/97 Use Beer's Law to calculate shading modifer
c.....LENIHAN 7/16/98 Calculate sliding extinction coeff k 

      call lacalc1(lait, rleavc)
      call lacalc2(lait, rlwodc, maxlai, klai)

      k = ((mx_index*0.01)*.70) + ((1.0-(mx_index*0.01))*.50)

c.....DB 10/8/98 canopy distribution reduction factor
c...  Ron's idea

      shdmod = exp(-k * lait)

      if (shdmod .lt. .001) shdmod = .001 

c... DB 10/10/98 Ron wanted LAI of 3.75, savannas to start favor grasses
c... problem is you lose shrublands. Jim had chosen 2 as threshold.

      if(lait .ge. 3.75) then
	percent_cover = 1.0
      else
        percent_cover = (1./3.75)*lait
      endif

      shdmod = shdmod*percent_cover + 1.0*(1.0-percent_cover)

c...Estimate the root/shoot ratio as a function of annual precip
c...added a trap for desert cells where precip is low

      if (frtc(1) .eq. 0.0 .and. grwprc .gt. 10.0) then
        rtsh = (bgppa + grwprc*bgppb)/(agppa + grwprc*agppb)
      else
        if (msplt .gt. frtc(3)) then
          fracrc = frtc(2)
        else
          fracrc = frtc(1) - msplt*(frtc(1)-frtc(2))/frtc(3)
        endif
	rtsh = fracrc/(1 - fracrc)
      endif


c...Estimate plant production:
      if (stemp .gt. 0.0) then

c...Change root/shoot ratio if burning occurs
	if (dofire(CRPSYS).or.dofire(SAVSYS)) then
          rtsh = rtsh + frtsh
        endif

c...Change root/shoot ratio by effect of CO2
	rtsh = rtsh * co2crs(CRPSYS)


c...Calculate temperature effect on growth

c...Removal of litter affects on soil temperature as it 
c     drives plant production (keith paustian)
        bio = aglivc * 2.5
        bio = min(bio,pmxbio)

c...Maximum temperature
        tmxs=maxtmp(month)+(25.4/(1.+18.*exp(-.20*maxtmp(month))))
     $             *(exp(pmxtmp*bio)-.13)

c...Minimum temperature
        tmns=mintmp(month)+pmntmp*bio-1.78

c...Average surface temperature
        ctemp=(tmxs+tmns)/2.

c...  Possible change in calculating the temperature effect on production
c...  DB 1-15-99

        potprd = gpdf(stemp, ppdf(1,1), ppdf(2,1), ppdf(3,1), 
     $                ppdf(4,1))
                                                                       
	c3potprd = gpdf(stemp, ppdf1(1,1), ppdf1(2,1), ppdf1(3,1),
     $                ppdf1(4,1))
	c4potprd = gpdf(stemp, ppdf2(1,1), ppdf2(2,1), ppdf2(3,1), 
     $                ppdf2(4,1))
      potprd = c3c4_index * 0.01 * c3potprd + (1-0.01*c3c4_index) *
     $   c4potprd

c...LENIHAN 4/28/05 
 
        gprd_tmp = potprd


c...Value for potential plant production is now calculated from the
c     equation of a line whose intercept changes depending on water
c     content based on soil type.  The function PPRDWC contains the
c     equation for the line.

        pet = pevap(month)
        if (pet .ge. .01) then

c...    testing sensitivity of the model to water availability
c...    between grasses and trees  DB/2-11-98
c...    more testing separating monthly rainfall as a function of lifeform canopy
c...    DB 11-20-98

c...	11-20-98 modifying role of water

	pptprd = (avh2o(1) + prcurr(month) + irract)/pet

c          pptprd = (avh2o(2) + prcurr(month) + irract)/pet

        else
         pptprd = .01
        endif
      if (pptprd .gt. 1.0) pptprd = 1.0
c...LENIHAN 4/28/05
 
        mi_grass = pptprd

c...LENIHAN 5/30/97

        mf_lmc(2) = pptprd
        if ( mf_lmc(2) .gt. 1.0 ) mf_lmc(2) = 1.0
        if ( pet .le. .01 ) mf_lmc(2) = 1.0

        call get_pprdwc3(pprdwcA, pprdwcB, pprdwcC, 0) 
        pptprd = pprdwc(wc, pptprd, pprdwcA, pprdwcB, pprdwcC)

c...LENIHAN 4/28/05

        gprd_ppt = pptprd

c...Calculate biof
	if (bioflg .eq. 1) then

c...Calculate maximum potential effect of standing dead on plant growth
c     (the effect of physical obstruction of litter and standing dead)
	  bioc = stdedc + .1*strucc(SRFC)
	  if (bioc .le. 0.) then
            bioc = .01
          endif

          if (bioc .gt. pmxbio) then
            bioc = pmxbio
          endif
	  bioprd = 1. - (bioc/(biok5+bioc))

c...Calculate the effect of the ratio of live biomass to dead biomass
c     on the reduction of potential growth rate.  the intercept of this 
c     equation ( highest negative effect of dead plant biomass ) is equal
c     to bioprd when the ratio is zero.
	  temp1 = (1. - bioprd)
	  temp2 = temp1*0.75
	  temp3 = temp1*0.25
	  ratlc = aglivc/bioc 
	  if (ratlc .le. 1.0) then
            biof = bioprd+(temp2*ratlc)
          endif
	  if (ratlc .gt. 1.0 .and. ratlc .le. 2.0) then
       	    biof = (bioprd+temp2)+temp3*(ratlc-1.)
          endif
	  if (ratlc .gt. 2.0) then
            biof = 1.0
          endif
	else
	  biof = 1.0
	endif

c...Restriction on seedling growth
c     sdlng is the fraction that prdx is reduced
        if (aglivc .gt. fulcan) then
          seedl = 0
        endif

        if (seedl .eq. 1) then
 	  sdlng = min(1.0, pltmrf + aglivc*(1-pltmrf) /fulcan)
        else
 	  sdlng = 1.0
        endif

c...Calculate potential production (biomass)

	agprod = prdx(1) * potprd * pptprd * biof * shdmod *
     +           sdlng * co2cpr(CRPSYS)

	bgprod = agprod * rtsh
	agp = agprod
	bgp = bgprod
	tgprod = agp + bgp
      else
        tgprod = 0.0
        agp = 0.0
        pcropc = 0.0
 
c-------C. Daly 12-19-97.  More initialized vars when stemp < 0.0
        bgp = 0.0
        potprd = 0.0
        pptprd = 0.0
        biof = 0.0
        sdlng = 0.0
       
        goto 10
      endif

c...Determine if grazing occurs
      if (dograz) then

c...Restrict production due to grazing:
c	grzeff = 0	grazing has no direct effect on production
c	grzeff = 1	linear impact on agp
c	grzeff = 2	quadratic impact on agp and root/shoot ratio
c	grzeff = 3	quadratic impact on root/shoot ratio
c	grzeff = 4	linear impact on root/shoot ratio
c	grzeff = 5	quadratic impact on agp and 
c                       linear impact on root/shoot ratio
c	grzeff = 6	linear impact on agb and root/shoot ratio


        if (grzeff .ne. 0) then
	  if (agp .le. 0.02) then
            agp = 0.02
          endif
	  if (grzeff .eq. 1) then
            agp = (1 - (2.21*flgrem)) * agp
            if (agp .lt. 0.02) then
              agp = 0.02
            endif
            bgp = rtsh * agp
          elseif (grzeff .eq. 2) then
            agp = (1 + 2.6*flgrem - (5.83*(flgrem**2)))*agp
            if (agp .lt. 0.02) then
              agp = 0.02
            endif
            bop = rtsh + 3.05*flgrem  - 11.78*(flgrem**2)
            if (bop .le. 0.01) then
              bop = 0.01
            endif
            bgp = agp * bop
          else if (grzeff .eq. 3) then
            bop = rtsh + 3.05*flgrem  - 11.78*(flgrem**2)
            if (bop .le. 0.01) then
              bop = 0.01
            endif
            bgp = agp * bop
          else if (grzeff .eq. 4) then
            bop = 1 - (flgrem * gremb)
            bgp = agp * bop
          else if (grzeff .eq. 5) then
            agp = (1 + 2.6*flgrem - (5.83*(flgrem**2)))*agp
            if (agp .lt. 0.02) then
              agp = 0.02
            endif
            bop = 1 - (flgrem * gremb)
            bgp = agp * bop
          else if (grzeff .eq. 6) then
            agp = (1 - (2.21*flgrem)) * agp
            if (agp .lt. 0.02) then
              agp = 0.02
            endif
            bop = 1 - (flgrem * gremb)
            bgp = agp * bop
          endif
	  tgprod = agp + bgp
          rtsh = bgp/agp
	endif
      endif

c...Compute the minimum and maximum C/N, C/P, and C/S ratios allowed
c     in plants.
c     Changed grwprc to 2.5*aglivc in calculation of cercrp

10    do 20 iel=1,nelem
        cercrp(IMIN,ABOVE,iel) =
     +          min(pramn(iel,1)+(pramn(iel,2)-pramn(iel,1)) *
     +          2.5 * aglivc / biomax,pramn(iel,2))
        cercrp(IMAX,ABOVE,iel) = 
     +          min(pramx(iel,1)+(pramx(iel,2)-pramx(iel,1)) *
     +          2.5 * aglivc / biomax,pramx(iel,2))

	cercrp(IMIN,BELOW,iel) = prbmn(iel,1)+prbmn(iel,2)*grwprc
        cercrp(IMAX,BELOW,iel) = prbmx(iel,1)+prbmx(iel,2)*grwprc
20    continue

c...If burning occurs, modify C/N ratio of shoots & roots.
      if (dofire(CRPSYS) .or. dofire(SAVSYS)) then
        cercrp(IMIN,ABOVE,N) = cercrp(IMIN,ABOVE,N) + 0.
        cercrp(IMAX,ABOVE,N) = cercrp(IMAX,ABOVE,N) + fnue(1)
        cercrp(IMIN,BELOW,N) = cercrp(IMIN,BELOW,N) + 0.
        cercrp(IMAX,BELOW,N) = cercrp(IMAX,BELOW,N) + fnue(2)
      endif

c...Added effect of co2
      do 30 iel = 1, nelem
  	cercrp(IMIN,ABOVE,iel) = cercrp(IMIN,ABOVE,iel) * 
     +                           co2cce(CRPSYS,IMIN,iel)
	cercrp(IMAX,ABOVE,iel) = cercrp(IMAX,ABOVE,iel) *
     +                           co2cce(CRPSYS,IMAX,iel)
30    continue

c...Update accumulators & compute potential C production
      ptagc = ptagc + agp/2.5
      ptbgc = ptbgc + bgp/2.5
      pcropc = tgprod / 2.5
c
      return
      end
