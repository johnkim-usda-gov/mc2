
c*****************************************************************
c	$RCSfile: cycle.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: cycle.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.3  2005/05/25 16:21:25  lenihan
cchanges for output of production scalars
c
cRevision 1.2  2004/09/24 20:16:09  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

        real function fimod(fi)
 
        real fi   

c...This funtion returns a value for an index that will modify soil T
c...The equation for the y-intercept (intcpt) is B+A*FI.  A and B
c     determine the effect of the frost index ie permafrost (non)existence
c     on soil temperature.   
		REAL :: fi1 = 0.5, fi2 = 0.67, fim1 = 1.0, fim2 = 0.5
      real slope, intercept
		slope = (fim1 - fim2)/(fi1 - fi2)
		intercept = fim1 - slope*fi1
       
       if (fi .le. fi1) then 
         fimod = fim1        
       elseif (fi .lt. fi2) then
         fimod = slope*fi + intercept	 
       else
         fimod = fim2
       endif
 
       return
       end


      subroutine cycle(month,cancvr,wfunc) 

      integer month
      real cancvr

c...Determine relative water content, available water, and
c     decomposition factor related to temperature and water.
c     Initialize mineralization accumulators.  Compute potential
c     production.

c...Changed 8-31-90:  added ELITST as a multiplier for SFCLIT.
c
c...Added savanna model (BO)

      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'wth.inc'
      include 'vemap2.inc'

c-----C. Daly 4-7-96.  Need to declare tcalc to get updated code to work.
      real tcalc
      external tcalc

c...Function declarations
      real anerob, irrigt,  pevap, fimod
      external anerob, irrigt, pevap

c...Local variables
      integer i, iel, lyr
      real aglivb, bio, co2val, pp, rprpet, sfclit, stdead, 
     -     tfunc, tmns, tmxs, wfunc, look 
	  real ttt, tt1, tt2, tt3, tt4, tt5, tt6, tt7
      character*80 string

c...Added new variable CANCVR for savanna model.  It will be passed
c     to POTGRS.  -rm 6/91

c...Add cursys .eq. 3 checks (BO)

c...Call schedl to determine scheduling options for this month
      call schedl()



c...Initialize production accumulators

c...For crops, annual accumulation starts at the month of planting.
c...For grass, annual accumulation starts at the beginning of
c     the growing season.
c...For forest, annual accumulation starts at the beginning
c     of the growing season.

      if ((dofrst .or. doplnt .or. dofone) .or.
     +		(month .eq. 1)) then
        call inprac
      endif

c...Average air temperature at 2 meters
      tave = (maxtmp(month) + mintmp(month)) / 2.0

c...Calculate RAIN for an output variable
      rain = prcurr(month)

c...If irrigating, determine actual amount of irrigation, irract
      if (doirri) then
        irract = irrigt(month)
      else
        irract = 0
      endif

c...If planting or growing month, calculate growing season precipitation
      if (doplnt .or. dofrst) then
        call prcgrw(month)
      endif

c...Initialize the mineralization accumulators for each element.
      do 30 iel = 1, MAXIEL
        do 25 lyr = 1, 2
          strmnr(lyr,iel) = 0.0
          metmnr(lyr,iel) = 0.0
          s1mnr(lyr,iel) = 0.0
25      continue
        s2mnr(iel) = 0.0
        s3mnr(iel) = 0.0
        gromin(iel) = 0.0
	w1mnr(iel) = 0.0
	w2mnr(iel) = 0.0
	w3mnr(iel) = 0.0
30    continue

c...Compute the ratio of precipitation to PET.
c   Added avh2o(3) to rprpet calculation  -rm 6/91

      pet=pevap(month)
      rprpet = (avh2o(3) + prcurr(month) + irract) / pet
c      write (*,*) 'avh2o(3), prcurr(month), irract, pet, rprpet = ',
c     .    avh2o(3), prcurr(month), irract, pet, rprpet
	  if (rprpet .ne. rprpet) then
            print *,
     $        'cycle: rprpet is a NaN. avh2o(3), prcurr(month) =',
     $        avh2o(3), prcurr(month)
	    print *, 'irract, pet, month = ', irract, pet, month
		stop
      endif

c...If the system is a microcosm, skip the rest of the routine
      if (micosm .eq. 1) then

        pet = 15.
        anerb = anerob(aneref, drain, rprpet, pet, tave, micosm)

        if (docult) then
          do 32 i = 1, 4
            cltfac(i) = clteff(i)
32        continue
        else
          do 33 i = 1, 4
            cltfac(i) = 1.0
33        continue
        endif

        return

      endif


c***********************************************************************

c...Calculate temperature...
      if (cursys .eq. FORSYS) then
c...Live biomass
        aglivb = rleavc * 2.5

c...Surface litter biomass
c...Second mod to remove effect of woodc -rm 1/91
        sfclit = (strucc(SRFC) + metabc(SRFC)) * 2.0

c...Standing dead biomass
        stdead = 0.0

      elseif (cursys .eq. SAVSYS) then
c...Live biomass
        aglivb = (rleavc + aglivc) * 2.5
        
c...Surface litter biomass
        sfclit = (strucc(SRFC) + metabc(SRFC)) * 2.0

c...Standing dead biomass
        stdead = stdedc * 2.5

      else
c...Live biomass
c...Removed addition of .25*pp to aglivb (AKM)
        aglivb = aglivc * 2.5

c...Surface litter biomass
        sfclit = (strucc(SRFC) + metabc(SRFC)) * 2.5

c...Standing dead biomass
        stdead = stdedc * 2.5
      endif

      bio = aglivb + stdead + elitst * sfclit
      bio = min(bio,pmxbio)

c-----C. Daly 12-18-97.  Constrain bio to 0.0 < bio < pmxbio to avoid
c-----math errors on exp function.
      bio = max(0.0,bio)
      
c...Maximum temperature
          tmxs=maxtmp(month)+(25.4/(1.+18.*exp(-.20*maxtmp(month))))
     1             *(exp(pmxtmp*bio)-.13)
      if (tmxs .ne. tmxs) then
        ttt = maxtmp(month)
	    tt1 = exp(pmxtmp*bio)
	    tt2 = exp(-.20*ttt)
	    tt3 = 25.4/(1. + 18.*tt2)
	    tt4 = ttt + tt3*(tt1 - 0.13)
c	    print *, 'ttt, tt1, tt2, tt3, tt4 = ', ttt, tt1, tt2, tt3, tt4
        tmxs = tt4
	    if (tmxs .ne. tmxs) then 
	      print *, 'tmxs = ', tmxs
          stop
	    endif
        print *,'*** cycle.F: NaN in tmxs has been fixed. tmxs = ', tmxs
      endif

c...Minimum temperature
      tmns=mintmp(month)+pmntmp*bio-1.78

c...Average surface temperature
c...Note: soil temperature used to calculate potential production does not
c         take into account the effect of snow (AKM)
      stemp=(tmxs+tmns)/2.
      look = fimod(frost_index)
      if (stemp .le. 0.0) then 
       stemp = stemp/look
      else
       stemp = stemp*look
      endif


c***********************************************************************


c...Shifted calculation of potential production so that it uses the current
c     month's soil temperature (AKM)

c...Determine potential production if it's during the growth season.
c...The variable, pp, is used below to recompute aglivb for use in h2olos.
      pp = 0.0

c... LENIHAN 4/28/05 initialize prd output vars
      fprd_ppt = 0.0
      fprd_tmp = 0.0
      gprd_ppt = 0.0
      gprd_tmp = 0.0
      mi_tree = 0.0
      mi_grass = 0.0

c...For a Crop System...
      if (crpgrw .eq. 1) then
        call potcrp(month,cancvr)
        pp = pcropc
      endif
	
c...For a Forest System...
      if (forgrw .eq. 1) then
        call potfor(month)
      endif
			

c***********************************************************************

c...Moisture
      if (cursys .eq. CRPSYS) then
        aglivb = aglivb + .25 * pp * 2.5
      endif

c...Determine co2 effect on transpiration, pass to h2olos
      if (cursys .eq. SAVSYS) then
        if (aglivc + rleavc .eq. 0.0) then
	  co2val = 1.0
        else
	  co2val = (co2ctr(CRPSYS)*aglivc + co2ctr(FORSYS)*rleavc) /
     +             (aglivc + rleavc)
        endif
      else
	co2val = co2ctr(cursys)
      endif

      call h2olos(month,aglivb,sfclit,stdead,co2val)

c...Average surface temperature
c     Added check for snow  -rm 8/91
      if (snow .gt. 0.0 .and. stemp .gt. 0.0) then
        stemp = 0.0
      endif

c...Effect of temperature on decomposition
      tfunc = tcalc(stemp, teff)

c...Effect of moisture on decomposition
c...Option selection for wfunc depending on idef
c        idef = 1 for linear option
c        idef = 2 for ratio option
c...Changed coefficients in both wfunc equations  -rm 6/91
c...Relative water content -- rwcf

      if (idef .eq. 1) then
        wfunc = 1./(1. + 4.0 * exp(-6.0*rwcf(1)))
      else if (idef .eq. 2) then
c-------C. Daly 4-23-96 catch underflows.
        if (-8.5 * rprpet .lt. -10.) then
          wfunc = 1.
        else
          wfunc = 1./(1. + 30.0 * exp(-8.5 * rprpet))
c          tt5 = -8.5*rprpet
c		  tt6 = exp(tt5)
c		  tt7 = 1./(1. + 30.*tt6)
c		  if (tt7.ne.wfunc) then
c		    print *, 'cycle l332: tt7, wfunc = ', tt7, wfunc
c	        print *, 'idef, rwcf(1), rprpet, tt5, tt6 = ',
c    $         idef, rwcf(1), rprpet, tt5, tt6
c			stop
c		  endif
        endif
c        write(*,*) 'rprpet, wfunc = ', rprpet, wfunc
      else
        call message(' ')
        string = '   IDEF may be either 1 or 2, but is currently: ' //
     +           char(ichar(char(idef)) + ichar('0'))
        call message(string)
        call message('   You must reset it and restart the run.')
        call message(' ')
        STOP
      endif
      if (wfunc .ge. 1.0) then
        wfunc = 1.0
      endif
      if (wfunc .ne. wfunc) then
	     print *, 'cycle l 342: wfunc is a NaN. rwcf, rprpet = ',
     $       rwcf, rprpet
	     print *, 'idef = ', idef
	     stop
      endif
c...Calculate the effect impact of anerobic conditions on decomposition
      anerb = anerob(aneref, drain, rprpet, pet, tave, micosm)

c...Combined effects of temperature and moisture
      defac=tfunc*wfunc
c      write(*,*) 'tfunc, wfunc, defac = ', tfunc, wfunc, defac
c...Bound defac to >= 0.0 12/21/92
       if (defac .lt. 0.0) defac = 0.0
       defacm(month)=defac
	   if (defac .ne. defac) then
         print *, '**** cycle: defac, tfunc, wfunc = ', 
     $       defac, tfunc, wfunc
		 stop
       endif
c...Effect of cultivation on decomposition (used in decomp routine)
c...Determine effect of cultivation on decomposition. vk 03-13-91
c...cltfac is this month's effect of cultivation on decomposition
c     of som1, som2, som3, and structural.  It is set to clteff
c     in months when cultivation occurs; otherwise it equals 1.
c...clteff is the effect of cultivation on decomposition read from
c     the cult.100 file

      if (docult) then
	do 34 i = 1, 4
          cltfac(i) = clteff(i)
34      continue
      else
	do 35 i = 1, 4
          cltfac(i) = 1.0
35      continue
      endif

      return
      end
