
c*****************************************************************
c	$RCSfile: simsom.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: simsom.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:44:13  bachelet
cheader
c
c******************************************************************
c
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine simsom()

c...Simulate flow of carbon, nitrogen, phosphorous, and sulfur.
c     This routine is executed each time step.  It calls the decomposition
c     submodel and a producer submodel.  It also includes a bunch of
c     N fixation stuff that needs to be rewritten and put in its own routine.
c
c     Added new local variable FSOL.  Added calls to new function FSFUNC
c     to calculate the amount of mineral P that is in solution.  Added
c     call to new subroutine PSCHEM, which calculates and schedules the
c     Phosphorus and Sulfur flows during decomposition.  Previously
c     this was calculated in the DECOMP routine.  -rm 6/91
      
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'isovar.inc'
      include 'ligvar.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'site.inc'
      include 't0par.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'zztim.inc'
      include 'vemap2.inc'

c...Function declarations
      real      fsfunc, rtimp, daylen
      external  fsfunc, rtimp, daylen

c...Local variables
      integer   iel, kts, lyr, iso
      real      biof, cancvr, cmn, 
     $          frlech(MAXIEL), fsol, fwdfx, fxbiom,
     $          satm, sirr, stot, tbiom, texeff, tmplen,
     $          wdbmas, wdfnp, wdfxm, wfunc


c...Added below for savanna model (rm)

      if (cursys .eq. SAVSYS ) then

        wdbmas = (fbrchc + rlwodc) * 2.0
        trbasl = wdbmas / basfct

c-------C. Daly 12-13-96.  Catch underflows. 
c-------Also made '1 - exp' '1.0 - exp'

        if (-0.064 * trbasl .lt. -10. ) then
          cancvr = 1.0
        else
          cancvr = 1.0 - exp(-0.064 * trbasl)
        endif

        if (trbasl .eq. 0.0) then
          trbasl = 0.1
        endif
      endif

c...Set aminrl for use in routines called from decomp
      do 10 iel = 1, nelem
        if (iel .eq. P) then
          fsol = fsfunc(minerl(1,P), pslsrb, sorpmx)
        else
          fsol = 1.0
        endif
        aminrl(iel) = minerl(1,iel) * fsol
10    continue

c...Determine decomposition factor and initialize accumulators
c     Added CANCVR to argument list.  -rm 6/91
      call cycle(month, cancvr, wfunc)
              
c...N Fixation
c     Does not take into account the effect of irrigation 
      if (nsnfix .eq. 1 .and. nelem .ge. P) then

c...Compute mineral N:P ratio for N-Fixation (use suface layer only)
c       rnpml1 is used to control soil N-fixation using a regression
c       equation based on Kansas data. This ratio is used only if nelem = 2.
c       rnpml1 is flagged if either minerl(1,1) or minerl(1,2) is zero.
c
        rnpml1 = minerl(1,N)/minerl(1,P)
     $           *fsfunc(minerl(1,P),pslsrb,sorpmx)

c...Wet-dry fixation of nitrogen -- monthly flow
c...Atmospheric fixation is split between monthly dry fall and
c     wdfnp is the N:P ratio control function for non-symbiotic
c     soil fixation.
c...Both wdfnp and fxbiom are relative functions
c     which vary from 0 to 1.
c...wdfnp computed as a negative natural log of rnpml1
c...symnfx is the symbiotic N fixation by legumes derived from Cole and
c     Heil (1981) using data from Walker et al. 1959. 
c... Bachelet NaN prevention

	if (rnpml1 .eq. 0.0 .or. fxnpb .eq. 0.0) then
	  wdfnp = 1.
	else
          wdfnp = min(1., ((-alog(rnpml1))/fxnpb)-.2)
	endif

c...The definitions of fxmca and fxmcb originally refered to water,
c     not biomass. (Alister 11/91)
	tbiom = aglivc+stdedc+strucc(SRFC)
	biof  = fxmca + tbiom * fxmcb
	fxbiom = 1 - biof
        fxbiom = min(1.,fxbiom)
        if (wdfnp .lt. 0 .or. fxbiom .lt. 0 .or. stemp .lt. 7.5) then
          fwdfx = 0.0
        else
          fwdfx = wdfnp * fxbiom
        endif

c...Compute N-fixation for the month

c...Wet fall depending on the monthly precipitation (wdfxma)
	if(prcann.eq.0.0) then 
                wdfxma = 0.
        else 
                wdfxma = wdfxa*   prcurr(month)/prcann
        endif
        wdfxms = fxmxs * fwdfx
        wdfxm  = wdfxma + wdfxms
		if (wdfxm .ne. wdfxm) then
                  print *,
     $        '*** simsom: wdfxm is a NaN. wdfxma, wdfxms = ',
     $        wdfxma, wdfxms
	      stop
        endif

c...Compute annual N-fixation accumulators for atmosphere and soils
        wdfxas = wdfxas + wdfxms
        wdfxaa = wdfxaa + wdfxma
        call flow(esrsnk(N),minerl(1,N),time,wdfxm)
        nfixac = nfixac+wdfxm

c...Monthly N-fixation based on annual parameters

      else
C...USE PRCURR/PRCANN INSTEAD OF DT
	if(prcann .eq. 0.0) then
		wdfxms = 0.
		wdfxma = 0.
	else
        	wdfxms = wdfxs*prcurr(month)/prcann
        	wdfxma = wdfxa*prcurr(month)/prcann
	endif
        wdfxas = wdfxas + wdfxms
        wdfxaa = wdfxaa + wdfxma
        wdfxm = wdfxma + wdfxms
        call flow(esrsnk(N),minerl(1,N),time,wdfxm)
        nfixac = nfixac + wdfxm
      endif

c...Monthly atmospheric S deposition
c commented out because minerl is now dimensioned (11,2) instead of (10,3)
c       if (nelem .eq. S) then
c	 if(prcann .eq. 0.0) then
c		satm = 0.0
c	 else
c         	satm = satmt * prcurr(month) / prcann
c	 endif
c         satmac = satmac + satm
c         if (doirri) then
c           sirr = sirri * irract * 0.01
c         else
c           sirr = 0
c         endif
c         sirrac = sirrac + sirr
c         stot = satm + sirr
c         call flow(esrsnk(S),minerl(SRFC,S),time,stot)
c       endif
         
c-----------------------------------------------------------

c...Decomposition Submodel

c...Determine  whether cultivation occurs this month.
c...Removed this call from DECOMP and put it here -rm 2/91
c...Removed it from here and put it in CYCLE.  -vk 3/91

c...Initialize stream variables for organic leaching (they are
c     monthly output).  -rm 3/92
      do 20 kts = 5, 8
        stream(kts) = 0.0
20    continue
      strm5l = 0.0
      strm5u = 0.0

c...Initialize monthly co2 accumlators (10/92)
      do 25 iso = 1, 2
        st1c2(iso) = 0.0
        st2c2(iso) = 0.0
        mt1c2(iso) = 0.0
        mt2c2(iso) = 0.0
        s11c2(iso) = 0.0
        s21c2(iso) = 0.0
        s2c2(iso)  = 0.0
        s3c2(iso)  = 0.0
25    continue

c...Call decomp routines ntspm times per month.
c...Removed the P and S chemistry from decomp and created the
c     subroutine pschem.  -rm  6/91
      do 40 kts = 1, ntspm
c        print *, 'calling decomp() now'
        call decomp(decodt,decsys)
c        print *, 'back from decomp() now'
	if (nelem .ge. P) then
          call pschem(decodt)
        endif

c...Update decomposition and nitrogen fixation flows.
        call flowup(time)
c	  print *, 'calling sumcar from simsom.F l239'
        call sumcar

c...aminrl contains the average amount of N, P, and S 
c     available in the top layer for the time period covered by 
c     dt/ntspm.  minerl contains the current value of mineral N,
c     P, and S by layer.
        do 30 iel = 1, nelem
	  if (iel .eq. P) then
            fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
          else
            fsol = 1.0
          endif
	  aminrl(iel) = (aminrl(iel) + minerl(SRFC,iel)*fsol)/2.0
30      continue
40    continue

c...Annual co2 accumulators (10/92)
      amt1c2 = amt1c2 + mt1c2(UNLABL) + mt1c2(LABELD)
      amt2c2 = amt2c2 + mt2c2(UNLABL) + mt2c2(LABELD)
      as11c2 = as11c2 + s11c2(UNLABL) + s11c2(LABELD)
      as2c2 = as2c2 + s2c2(UNLABL) + s2c2(LABELD)
      as3c2 = as3c2 + s3c2(UNLABL) + s3c2(LABELD)
      ast1c2 = ast1c2 + st1c2(UNLABL) + st1c2(LABELD)
      as21c2 = as21c2 + s21c2(UNLABL) + s21c2(LABELD)
      ast2c2 = ast2c2 + st2c2(UNLABL) + st2c2(LABELD)

c...Volatilization loss of nitrogen as a function of
c     gross mineralization
      volgm = vlossg*gromin(1)

c...desperate move by DB 3-12-99
      if(volgm .lt.0.) volgm = 0.0

c... added a trap to make sure minerl remains >0 DB 10-13-98
c      print *, "time, minerl(SRFC,N) = ", time, minerl(SRFC, N)
      if(minerl(SRFC,N) .lt. volgm) volgm = minerl(SRFC,N)
      minerl(SRFC,N) = minerl(SRFC,N) - volgm
      esrsnk(N) = esrsnk(N) + volgm

c...Set frlech to leaching fraction vek june90
c     Compute normal value for frlech.  Recompute in flood routine
c     if flooding occurs.
      texeff = fleach(1) + fleach(2) * sand
      do 50 iel = 1, nelem
        if (iel .eq. P) then
          fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
        else
          fsol = 1.0
        endif
        frlech(iel) = texeff * fleach(iel+2) * fsol
50    continue

c...Soil erosion
      if (doerod) then
        call erosn(psloss,bulkd,edepth,enrich,lhzci,lhze,nelem)
      else
        scloss = 0.0
      endif

c...Fertilization option
      if (dofert) then
        do 60 iel = 1, nelem
          esrsnk(iel) = esrsnk(iel)-feramt(iel)
          minerl(SRFC,iel) = minerl(SRFC,iel)+feramt(iel)
          fertot(iel) = fertot(iel) + feramt(iel)
60      continue
      endif

c...Available nutrients
c     tminrl is the total amount of each element available in
c     mineral form.
      do 80 iel = 1, nelem
        tminrl(iel) = 0.
        if (iel .eq. P) then
          fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
        else 
          fsol = 1.0
        endif
	       
        do 70 lyr = 1, nlayer

c...Plants can only uptake from a layer with a positive
c     value, so only the positive layers are summed here.
          if (minerl(lyr,iel) .gt. 0.)  then
     		tminrl(iel) = tminrl(iel) + minerl(lyr,iel) * fsol
          endif
70        continue
80      continue

c...Compute the fraction of labile (non-sorbed) P in the surface
c     layer available to plants

      favail(2) = max(favail(4),
     $                min(favail(4) + minerl(SRFC,N)*
     $                   (favail(5) - favail(4)) / favail(6),
     $                    favail(5)))


c...Add to fallow rain
      if (falprc .eq. 1 .and. .not. dohrvt) then
        prcfal = prcfal + prcurr(month)
      endif


c...Call the producer submodel

c...Compute daylength for use in phenology of trees
      tmplen = daylen(month, sitlat)

c...Determine if daylength is increasing or decreasing
      if (tmplen .lt. dayhrs) then
	 hrsinc = .FALSE.
      else
	 hrsinc = .TRUE.
      endif

      dayhrs = tmplen


c...Crop and Forest removal options - moved here from CROP
c     and TREES so the mineral pools are not radically changed
c     between the growth routines. - rm 7/94


c...Fire and grazing
        if (dofire(SAVSYS) .or. dograz) then
          call grem()
        endif
c...Fire or cutting events
        if (dotrem) then
          call frem()
        endif

        
c...Update state variables and accumulators and sum carbon isotopes.
      call flowup(time)
c	  print *, 'calling sumcar from simsom.F l375'
      call sumcar

      call crop(time, wfunc)
      call trees(month,cisotf, wfunc)


c...Harvest may be performed after updating flows.  Put here for
c     consistency with the Savanna model - moved calls to flowup, 
c     sumcar and harvst from CROP routine to here. -rm 7/94

      if (dohrvt) then
        call harvst(month,pltlig)
      endif




c...Leaching
      call leach(amov, nelem, nlayer, minerl, minlch, frlech, stream,
     $           basef, stormf)

c...Update state variables and accumulators and sum carbon isotopes.
      call flowup(time)
c	  print *, 'calling sumcar from simsom.F l399'
      call sumcar

c...Accumulate leached C,N,P,S
      csrsnk(UNLABL) = csrsnk(UNLABL) + strm5u
      csrsnk(LABELD) = csrsnk(LABELD) + strm5l
      stream(5) = strm5u + strm5l
      do 90 iel = 1, nelem
        esrsnk(iel) = esrsnk(iel) + stream(iel+1) + stream(iel+5)
90    continue

c...Volatilization loss as a function of the mineral n which
c     remains after uptake by plants
c...another desperate move by DB here 3-12-99
c...if (minerl(SRFC,N) .gt. 0.0) then

      if (minerl(SRFC,N) .gt. 0.1) then
        volex = vlosse*minerl(SRFC,N)*dt
        minerl(SRFC,N) = minerl(SRFC,N) - volex
        esrsnk(N) = esrsnk(N) + volex
      endif

c...Volatilization
      volgma = volgma+volgm
      volexa = volexa+volex
c      print *, "*** volgma, volgm = ", volgma, volgm

c...Production
      cproda = cproda + cprodc + cprodf

c...Net Mineralization
      do 100 iel = 1, nelem

c...Net mineralization for the mineralizing compartments
c...The structural component of litter and the wood compartments
c     are not mineralizers.  They should not be added into cmn or
c     sumnrs.
        cmn = metmnr(SRFC,iel) + metmnr(SOIL,iel) + 
     $        s1mnr(SRFC,iel) + s1mnr(SOIL,iel) +
     $        s2mnr(iel) + s3mnr(iel)
	sumnrs(iel) = sumnrs(iel) + cmn

c...soilnm is net mineralization in the soil.
        soilnm(iel) = soilnm(iel) + s1mnr(SOIL,iel) +
     $                s2mnr(iel) + s3mnr(iel) +
     $                metmnr(SOIL,iel) + strmnr(SOIL,iel) + w3mnr(iel)

c...Total net mineralization
        tnetmn(iel) = tnetmn(iel) + cmn + 
     $                strmnr(SRFC,iel) + strmnr(SOIL,iel) +
     $                w1mnr(iel) + w2mnr(iel) + w3mnr(iel)
100   continue

c...Compute output variables for printing or plotting.
      call savarp

      bio_consume = shrema + sdrema + tcrem

      return
      end
