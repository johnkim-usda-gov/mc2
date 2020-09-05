c*****************************************************************
c	$RCSfile: eachyr.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: eachyr.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:48:15  conklind
cUse mapssvar.inc instead of daly.inc; eliminate unused variables yr, ndep_scalar, and ndep_diff.
c
cRevision 1.1  2003/11/21 00:14:01  conklind
cInitial revision
c
c******************************************************************

c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine eachyr()

c...Perform tasks that only need to be done once a year.

      include 'const.inc'
      include 'mapssvar.inc'
      include 'isovar.inc'
      include 'ligvar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'seq.inc'
      include 'wth.inc'
      include 'zztim.inc'

c...Function declarations
      real     fracis, pevap 
      external fracis, pevap 

c...Local variables
      integer  iel, ipart, mon
c     integer  itime
      real     lfncmax, lfncmin, lfncon


c...Added for savanna model (plot3.inc) BO

c...Correct for the rounding problem with time. The value for time
c     drifts downward during long runs since dt=1/12 cannot be represented
c     precisely.  At this point in the run, time should be a whole number.

c...Changed increment value from .001 to .5 to correct error in time calculation
c     occuring after year 8192. (mse 3/95).  New code if from Kendrick Killian.

      time = sign(int(abs(time)+.5),int(time))

c     if (time .lt. 0) then
c       itime = int(time - .001)
c     else
c       itime = int(time + .001)
c     endif

c     time = real(itime)

c...Reset annual accumulators to zero
      call annacc

c...Weather data
      call weathr(precip,prcstd,prcskw,mintmp,maxtmp)

c...Wet-dry fixation of N

c...Determine annual precipitation and annual PET
      prcann = 0.0
      petann = 0.0
      nfixac = 0.0
      wdfxs = 0.0
      wdfxa = 0.0
      do 10 mon = 1, MONTHS
        prcann = prcann + prcurr(mon)
        petann = petann + pevap(mon)
10    continue

c...N fixation in atmosphere
c... DB 10-16-98 Saturating the soil with N using the N deposition 
c... DB 4-19-99 limiting the amount of prcann to reflect "new" century code

      wdfxa = epnfa(INTCPT)+epnfa(SLOPE)*MIN(prcann,80.0)

c...DB MAY2003  removing the cap on N dep as a function of precip
c      wdfxa = epnfa(INTCPT)+epnfa(SLOPE)*prcann

c      wdfxa = 5.

c... reading Thornton's N dep files  4-20-99
c... assuming co2 in 1895 is 294 and in 1990 353
c... commenting it ut oct 29, 2002

c       yr = time
c       if(co2val(yr) .gt. 294.842) then
c         ndep_scalar = (cdp - hdp) / (353.636 - 294.842)
c         ndep_diff = (co2val(yr) - 294.842) * ndep_scalar
c       else
c         ndep_diff = 0.0
c       endif

c       wdfxa = hdp + ndep_diff

c... don't allow the industrial ndep levels to be less than the preindustrial levels 

c       if (wdfxa .lt. hdp) wdfxa = hdp

c.. this is a test DB 4-27-99
c.. back to version13 style: no MIN in the equation
 
c       wdfxa = epnfa(INTCPT)+epnfa(SLOPE)*prcann
c       wdfxa = epnfa(INTCPT)+epnfa(SLOPE)*MIN(prcann,80.0)

       if (wdfxa .lt. 0.) then
         wdfxa = 0.0
       endif

c... running with no N limitation we don't need deposition 
c... right: oct 29,2002
c         wdfxa = 0.
c... running same code as century aug2000, wdfxa not 0.0

c...Use annual ET unless it is the first timestep
c...No longer using the intercept in the calculation.

c      wdfxs = epnfs(INTCPT)+epnfs(SLOPE)*MIN(prcann,100.0)

c...Adding "new" century code with abiotic fixation as a function  of aet
c...DB 4/19/99 added annet to param.inc
c...Code specific to GRIDDED CENTURY
c...Code specific to POTSDAM - using ANNET
c... commenting it out when demand is met oct 29, 2002
c... DB removing cap on wdfxs MAY 2003

      if (annet .eq. 0.0) then
        wdfxs = epnfs(INTCPT)+epnfs(SLOPE)*MIN(prcann,100.0)
        wdfxs = epnfs(SLOPE)*MIN(prcann,100.0)
c        wdfxs = epnfs(INTCPT)+epnfs(SLOPE)*prcann
      else 
	 wdfxs = epnfs(2) * (annet - epnfs(1))
      endif

c... zeroing wdfxs to match Peter Thornton's inputs
c... DB 4-22-99   and oct 29,2002
c	wdfxs = 0.

c...Reset annual accumulator for evapotranspiration
c...Code specific to GRIDDED CENTURY
      annet = 0

      if (wdfxs .lt. 0.)  then
        wdfxs = 0.0
      endif
      wdfx = wdfxa+wdfxs

c...Atmospheric S deposition
      satmt = max(0.0, satmos(1) + satmos(2)*prcann)

c...Determine what fraction of the carbon in new plant tissue is labeled
      if (labtyp .eq. 0) then
        cisofr = 0.0
        cisotf = 0.0
      elseif (labtyp .eq. 1) then
c        cisofr = fracis(time,labyr)
         stop 'Unexpected code for MC1'
        cisotf = cisofr
c     elseif (labtyp .eq. 2) then
c       cropin has set cisofr
c       treein has set cisotf
      endif

c...Initialize co2 effects
      call co2eff(time)

c...Added effect of co2 for forest; done here because not calcualted
c     dynamically based on biomass like grassland/crop
c...Direct CO2 effects only C/E ratio of leaves.
      do 30 iel = 1, nelem
        ccefor(IMIN,LEAF,iel) = cerfor(IMIN,LEAF,iel) *
     $                            co2cce(FORSYS,IMIN,iel)
        ccefor(IMAX,LEAF,iel) = cerfor(IMAX,LEAF,iel) *
     $                            co2cce(FORSYS,IMAX,iel)
 30   continue

      do 50 ipart = 2, FPARTS
        do 40 iel = 1, nelem 
          ccefor(IMIN,ipart,iel) = cerfor(IMIN,ipart,iel)
          ccefor(IMAX,ipart,iel) = cerfor(IMAX,ipart,iel)
 40     continue 
 50   continue


c...Calculate leaf death rate multiplier for continuous forests 11/20/92
c...Initialize LDRMLT to 1.0
      ldrmlt = 1.0

c...Change leaf death rate multiplier if you have floating C/E ratios.
      if ((ccefor(IMIN,LEAF,N) .ne. ccefor(IMAX,LEAF,N)) .and.
     $    (ccefor(IMIN,LEAF,N) .ne. 0) .and. 
     $    (ccefor(IMAX,LEAF,N) .ne. 0) .and. 
     $    (1/ccefor(IMAX,LEAF,N) - 1/ccefor(IMIN,LEAF,N))
     $    .ne. 0) then
	  if (rleavc .gt. 0) then
	  lfncon = rleave(N) / rleavc
	  lfncmin = 1 / ccefor(IMIN,LEAF,N)
	  lfncmax = 1 / ccefor(IMAX,LEAF,N)
	  ldrmlt = 1 + (maxldr - 1) *
     $             (lfncon - lfncmin) / (lfncmax - lfncmin)
        endif
      endif

      if (cursys .ne. FORSYS) then

c...Determine what fraction of plant residue added this year
c     will be lignin.
        call cmplig(cursys,fligni,wdlig,pltlig)
      endif

      return
      end
