      
c*****************************************************************
c	$RCSfile: csched.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: csched.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:14:10  bachelet
cheader
c
c******************************************************************
	subroutine csched(cflow,protop,probot,
     $                  froma1,tob1,
     $                  froma2,tob2,
     $                  frac,accum)

c...Schedule C flows for decomposition from Box A to Box B
c...written by vek 05/91

      real    cflow, protop, probot, froma1, tob1, froma2, tob2,
     $        frac,accum(2)

c...Input:
c     cflow  = total C flow from Box A to Box B
c     protop = proportion top, either the labeled pool 
c                OR cisofr or cisotf
c     probot = proportion bottom, either the total pool or 1
c     frac   = amount of fractionation
c     time   = simulation time (passed in /zztim/)

c...Transput:
c     froma1 = state variable for unlabeled C in Box A
c     tob1   = state variable for unlabeled C in Box B
c     froma2 = state variable for labeled C in Box A
c     tob2   = state variable for labeled C in Box B
c     accum(1) = unlabeled accumulator
c     accum(2) = labeled accumulator

      include 'const.inc'
      include 'zztim.inc'

c...Local variables
      real    cflow1, cflow2


c...Determine the amount of labeled to flow
c----Changed to prevent NaN.********
      if(probot .ne. 0.0) then
        cflow2 = cflow * (protop / probot) * frac
      else
        return 
      endif

c...Determine the amount of unlabeled to flow
      cflow1 = cflow - cflow2

c...Flow the amounts
      if (cflow1 .ne. cflow1) then
          print *, "cflow1 is a NaN" 
      endif
      if (cflow2 .ne. cflow2) then
          print *, "cflow2 is a NaN"
      endif    
      call flow(froma1,tob1,time,cflow1)

      call flow(froma2,tob2,time,cflow2)

c...Accumulate the amounts flowed
      accum(UNLABL) = accum(UNLABL) + cflow1
      accum(LABELD) = accum(LABELD) + cflow2

      return
      end
