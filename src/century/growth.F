
c*****************************************************************
c	$RCSfile: growth.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: growth.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.3  2004/12/13 21:20:04  bachelet
cgetting rid of calculation for growth (deleted it from declaration stmt)
c
cRevision 1.2  2004/09/24 20:43:32  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine growth(cisofr)

      real      cisofr

c...Simulate production for the month.

      include 'comput.inc'
      include 'const.inc'
      include 'fertil.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'zztim.inc'
      include 'dovars.inc'

c...Function declarations
      real      fsfunc, rtimp
      external  fsfunc, rtimp

c...Local variables
      integer   iel, lyr
      real      agfrac, amt, avefrc(MAXIEL), bgfrac,
     $          calcup, cfrac(CPARTS),
     $          euf(CPARTS), fsol,
     $          gnfrac, growag, growbg,
     $          rimpct, tm, uptake(4,MAXIEL)


c...Determine actual production values, restricting the C/E ratios
      if (crpgrw .eq. 1 .and. pcropc .gt. 0.0) then

c...Calculate impact of root biomass on available nutrients
        rimpct = rtimp(riint, rictrl, bglivc)

c...Calculate carbon fraction in each part
        cfrac(ABOVE) = agp / tgprod
        cfrac(BELOW) = 1 - cfrac(ABOVE)

c...Calculate savanna available fractions
        if (cursys .eq. SAVSYS) then
          if (tminrl(N) .gt. 1.5) then
            tm = 1.5
          else
            tm = tminrl(N)
          endif

c---------C. Daly 12-13-96.  Catch underflows.

          if (-.00102*tm*sitpot .lt. -10. ) then
            gnfrac = 0.0
          elseif (-1.664*exp(-.00102*tm*sitpot)*basfc2*trbasl 
     .            .lt. -10) then
            gnfrac = 0.0
          else
            gnfrac = exp(-1.664*exp(-.00102*tm*sitpot)*basfc2*trbasl)
          endif

c...Bound GNFRAC between 0 and 1
          if (gnfrac .gt. 1.0) then
            gnfrac = 1.0
          elseif (gnfrac .lt. 0.0) then
            gnfrac = 0.0
          endif

	    avefrc(N) = tminrl(N) * gnfrac
          else
  	  avefrc(N) = tminrl(N)
          endif

        avefrc(P) = tminrl(P)
        avefrc(S) = tminrl(S)


c...2003 DB checking the year

	call restrp(nelem,CPARTS,avefrc,cercrp,cfrac,pcropc,rimpct,
     $              crpstg,snfxmx(CRPSYS),cprodc,eprodc,eup,
     $              uptake,elimit,nfix,relyld,snfxac(CRPSYS),svuptk)

c...If growth occurs...
        if (cprodc .gt. 0.) then

c...Update accumulators for N, P, and S uptake
          do 20 iel = 1, nelem
            eupacc(iel) = eupacc(iel) + eprodc(iel)
c           eup(ABOVE) is the fraction of element allocated to aboveground
c           eup(BELOW) is the fraction of element allocated to belowground
            eupaga(iel) = eupaga(iel) + eup(ABOVE,iel)
            eupbga(iel) = eupbga(iel) + eup(BELOW,iel)
20        continue

c...C/N ratio for production
c          tcnpro = cprodc/eprodc(N)

c...Growth of shoots
          agfrac = agp/tgprod
          growag = cprodc * agfrac
          call csched(growag,cisofr,1.0,
     $                csrsnk(UNLABL),aglcis(UNLABL),
     $                csrsnk(LABELD),aglcis(LABELD),
     $                1.0,agcisa)

c...Growth of roots
          bgfrac = 1. - agfrac
          growbg = cprodc * bgfrac
          call csched(growbg,cisofr,1.0,
     $                csrsnk(UNLABL),bglcis(UNLABL),
     $                csrsnk(LABELD),bglcis(LABELD),
     $                1.0,bgcisa)

c...Actual uptake
          do 40 iel = 1, nelem
	    euf(ABOVE) = eup(ABOVE,iel) /eprodc(iel)
	    euf(BELOW) = eup(BELOW,iel) /eprodc(iel)

c...Take up nutrients from internal storage pool
c... adding a check on positive pools DB 10-14-98
	   if(crpstg(iel) .gt. 0.0) then
	    amt = uptake(ESTOR,iel) * euf(ABOVE)
	    call flow(crpstg(iel),aglive(iel),time,amt)
	    amt = uptake(ESTOR,iel) * euf(BELOW)
	    call flow(crpstg(iel),bglive(iel),time,amt)
	   endif

c...Take up nutrients from soil
            do 30 lyr = 1, nlayer
              if (minerl(lyr,iel) .gt. 0.) then
                fsol = 1.0
                if (iel .eq. P) then
                  fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
                endif
                calcup = uptake(ESOIL,iel) *
     $                   minerl(lyr,iel) * fsol / tminrl(iel)
	        amt = calcup * euf(ABOVE)
c... added by DB 10-2-98
		if(minerl(lyr,iel).gt. amt) then
	        call flow(minerl(lyr,iel),aglive(iel),time,amt)
		endif
	        amt = calcup * euf(BELOW)
c... added by DB 10-2-98
		if(minerl(lyr,iel).gt. amt) then
	        call flow(minerl(lyr,iel),bglive(iel),time,amt)
		endif
	      endif
30	    continue

c...Take up nutrients from nitrogen fixation
	    if (iel .eq. N .and. nfix .gt. 0.) then
	      amt = uptake(ENFIX,iel) * euf(ABOVE)
c		if(esrsnk(iel).gt. amt) then
              call flow(esrsnk(iel),aglive(iel),time,amt)
c		endif
 	      amt = uptake(ENFIX,iel) * euf(BELOW)
c		if(esrsnk(iel).gt. amt) then
	      call flow(esrsnk(iel),bglive(iel),time,amt)
c		endif
	    endif

c...Take up nutrients from fertilizer
	    if (aufert .ne. 0 .and. uptake(EFERT,iel) .gt. 0.0) then

c...Automatic fertilizer added to plant pools
              amt = uptake(EFERT,iel) * euf(ABOVE)
              call flow(esrsnk(iel),aglive(iel),time,amt)
              amt = uptake(EFERT,iel) * euf(BELOW)
              call flow(esrsnk(iel),bglive(iel),time,amt)

c...Automatic fertilizer added to mineral pool
              amt = uptake(EFERT,iel) * (1./favail(iel) - 1.)  
              fertot(iel) = fertot(iel) + uptake(EFERT,iel) + amt
              call flow(esrsnk(iel),minerl(SRFC,iel),time,amt)
            endif
40        continue
	endif

c...Else no production this month
      else
	cprodc = 0.0
        do 50 iel = 1, nelem
          eprodc(iel) = 0.0
50      continue
        nfix = 0.
      endif

      return
      end
