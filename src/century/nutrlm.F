
c*****************************************************************
c	$RCSfile: nutrlm.F,v $
c	$Revision: 1.3 $
c	$Date: 2004/09/24 20:55:49 $
c	$Locker:  $
c	
c$Log: nutrlm.F,v $
cRevision 1.3  2004/09/24 20:55:49  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

C...NUTRLM.F
c
c...Nutrient limitation for plants is based on demand

      subroutine nutrlm(nelem, nparts, cprodl, eprodl, maxec,
     $                  maxeci, mineci, cfrac, eavail, nfix,
     $                  snfxmx, snfxac, elimit, eup)

      include 'const.inc'
      include 'plot4.inc'
	  include 'mapssvar.inc'

      integer nelem, nparts
      real cprodl, elimit, eprodl(MAXIEL), maxec(MAXIEL), 
     $     eavail(MAXIEL), nfix, snfxmx, snfxac,
     $     eup(nparts,MAXIEL), mineci(nparts,MAXIEL),
     $     maxeci(nparts,MAXIEL), cfrac(nparts)

c...Local Variables
c   NOTE:  Local variables cannot have adjustable array size.  ECFOR
c          is set to the largest array size which may occur.
      integer 	iel, ipart
      real 	mxnfix, cpbe(MAXIEL), demand, ecfor(5,MAXIEL)

c...Definitions of Local variables
c	mxnfix		- Maximum N fixation
c	cpbe		- Carbon production limited by element
c       ecfor		- Actual E/C ratio by part 

c...Compute production limitation

      do 10 iel = 1, nelem

c...Initialize fixation to 0
	mxnfix = 0.0
c...DEMAND based on the maximum E/C ratio.
        demand = cprodl * maxec(iel)

c...N FIXATION
      if (iel .eq. N) then
        if (unl_N_flag) then
          mxnfix = demand
        else
          mxnfix = snfxmx * cprodl
        endif
        symb = mxnfix
      endif

c...New C/E ratios by part based on E available.
c...new code from Melannie 8/8/00 in case eavail <0
c...not quite ready for show time yet

c        if(eavail(iel) .lt. 0.0) then
c          do 18 ipart = 1, nparts
c            ecfor(ipart,iel) = mineci(ipart, iel)
c 18       continue

c        else if ((eavail(iel) + mxnfix) .gt. demand) then

        if ((eavail(iel) + mxnfix) .gt. demand) then
	  do 20 ipart = 1, nparts
	    ecfor(ipart,iel) = maxeci(ipart,iel)
 20       continue

	else
	  do 30 ipart = 1, nparts
	    ecfor(ipart,iel) = mineci(ipart,iel) +
     $                    (maxeci(ipart,iel) - mineci(ipart,iel))
     $			* (eavail(iel) + mxnfix) / demand

 30       continue
	endif

c...Initialize local variables to zero
	cpbe(iel) = 0.0

c...Total potential production with nutrient limitation
	do 40 ipart = 1, nparts
	  cpbe(iel) = cpbe(iel) + cfrac(ipart) * ecfor(ipart,iel)
 40     continue
c.. updated with new code from melannie 8/16/00 in case cpbe goes <0
c... again not quite ready for showtime

c	  cpbe(iel) =  max(0.0, (eavail(iel) + mxnfix) / cpbe(iel))
          cpbe(iel) = (eavail(iel) + mxnfix) / cpbe(iel)

c...Put automatic fertilization here when necessary
c	if (aufert .gt. 0)

 10   continue

c...Compute the limiting element
c...added clause that elimit needs initializing 
c...if no element is limiting DB 1998

      do 50 iel = 1, nelem
        if (cprodl .gt. cpbe(iel)) then 
	  cprodl = cpbe(iel)
	  elimit = 1.
	else
	  elimit = 0.
        endif
 50   continue

c...Recompute EPRODL
      do 70 iel = 1, nelem
        eprodl(iel) = 0.0

c...Total potential production with nutrient limitation
	do 60 ipart = 1, nparts
	  eup(ipart,iel) = cprodl * cfrac(ipart) *
     $                       ecfor(ipart,iel)
	  eprodl(iel) = eprodl(iel) + eup(ipart,iel)
 60     continue

c...Check to make sure the total flow won't exceed what's available.
c...DEBUG
	  if (eprodl(N) - (eavail(N) + mxnfix) .gt. .000001) then
	     call message('In NUTRLM: EPRODL > FIX + AVAIL')
	     write(*,*) 'NPARTS: ', nparts, ' EPRODL: ', eprodl(iel)
	     write(*,*) 'EAVAIL: ', eavail(iel) 
	     write(*,*) 'MXNFIX: ', mxnfix
c	     stop
         endif

 70   continue

c...Compute N fixation which actually occurs and add to the
c   N fixation accumulator.
      nfix = max(eprodl(N) - eavail(N), 0.0)
      snfxac = snfxac + nfix

      return
      end

