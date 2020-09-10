
c*****************************************************************
c	$RCSfile: frem.F,v $
c	$Revision: 1.2 $
c	$Date: 2004/09/24 20:29:05 $
c	$Locker:  $
c	
c$Log: frem.F,v $
cRevision 1.2  2004/09/24 20:29:05  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

c...FREM
      subroutine frem()

c...Forest removal - fire or cutting (cutting includes storms)
c   Includes litter burning in forest systems.

c...Called from:	simsom

      include 'const.inc'
      include 'forrem.inc'
      include 'param.inc'

c...Local variables
      real	accum(ISOS)

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c_______LENIHAN 4/28/97

      if ( did_burn ) then
        remf(1)   = mf_remf(1)
        remf(2)   = mf_remf(2)
        remf(3)   = mf_remf(3)
        remf(4)   = mf_remf(4)
        remf(5)   = mf_remf(5)
      endif

c...Removal for both CUT and FIRE events
      if ((evntyp .eq. 0) .or. (evntyp .eq. 1)) then

c...Live Removal
	 call livrem(accum)

c...Dead Removal
	 call dedrem(accum)

c...Death of Roots
	 call killrt(accum)

      endif

c...Returns from cutting event
      if (evntyp .eq. 0) then
	 call cutrtn(accum)

c...Returns from fire event
      else if (evntyp .eq. 1) then
	 call firrtn()

      endif

      return
      end

