	
c*****************************************************************
c	$RCSfile: co2eff.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: co2eff.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:02:15  bachelet
chedaer
c
c******************************************************************
	subroutine co2eff(time)
	real time

c...Compute the effect of atmospheric CO2 concentration.

	include 'const.inc'
	include 'param.inc'
	include 'parfx.inc'
	include 'plot1.inc'
        include 'mapssvar.inc'

c...Function declarations
	real      effect, line, ramp
	external  line, ramp

c...Local variables
	integer   iel, mnmx, system, yr
	real      co2conc


c...Reset all effects to 1.0
	do 30 system = CRPSYS, FORSYS
	  co2cpr(system) = 1.0
	  co2ctr(system) = 1.0
	  co2crs(system) = 1.0
	  do 20 mnmx = IMIN, IMAX
	    do 10 iel = 1, nelem
	      co2cce(system,mnmx,iel) = 1.0
10	    continue
20	  continue
30	continue

c...If there is no co2 effect, return; 
c     otherwise, calculate a new co2 effect
	if (co2sys .lt. 0) then
	  return
	endif

c...Calculate a co2 concentration
	if (co2rmp .eq. 0) then
c...Steady
          co2conc = co2ppm(1)
	else
c...Ramping
cx	  co2conc = ramp(time,co2tm(1),co2ppm(1),co2tm(2),co2ppm(2))
          yr = time
c         co2conc = co2val(yr)
          call get_co2(yr, co2conc)

c         print *,'*** co2eff() time = ', time, '  co2conc = ', co2conc
	endif
	if (co2conc .le. 0) then
	  print *,'*** co2eff.f: co2conc .le. 0'
	  stop
	endif

c...Calculate effect on production
        do 40 system = CRPSYS, FORSYS
	  co2cpr(system) = effect(co2ipr(system),co2conc)
40      continue
c      print *,'*** co2ipr(FORSYS), co2cpr(FORSYS) = ', 
c     $    co2ipr(FORSYS), co2cpr(FORSYS)

c...Calculate effect on PET
        do 50 system = CRPSYS, FORSYS
	  co2ctr(system) = effect(co2itr(system),co2conc)
50	continue
c      print *,'*** co2itr(FORSYS), co2ctr(FORSYS) = ', 
c     $    co2itr(FORSYS), co2ctr(FORSYS)

c...Calculate effect on C/E
	do 80 iel = 1, nelem
	  do 70 mnmx = IMIN, IMAX
            do 60 system = CRPSYS, FORSYS
	      co2cce(system,mnmx,iel) =
     $                 effect(co2ice(system,mnmx,iel),co2conc)
60          continue
70        continue
80      continue

c...Calculate effect on root/shoot
c     Reference co2 concentration = 350.0 at 1.0
        do 90 system = CRPSYS, FORSYS
	  co2crs(system) = line(co2conc,350.0,1.0,700.0,co2irs(system))
90      continue

	return
	end




	real function effect(co2input, co2conc)

	real      co2input, co2conc

c...Reference co2 concentration = 350.0
	effect = 1 + (co2input-1) / (log10(2.0))*(log10(co2conc/350.0))

	return
	end
