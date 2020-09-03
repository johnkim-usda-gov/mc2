
c*****************************************************************
c	$RCSfile: grochk.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: grochk.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:40:25  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

C...GROCHK.F

      integer function grochk(tave)

c...This function determines whether or not there is
c   potential growth and if it is the first month of
c   the growing season.
c
c   grochk = 1:  Month of greenup in deciduous forests
c   grochk = 0:  No greenup
c

      include 'pheno.inc'
      include 'timvar.inc'
      include 'zztim.inc'

      real	tave


c...If it is spring and the temperature is high enough and
c     you haven't already reapportioned carbon to leaves...
      if (hrsinc .and. (tave .gt. 13.0)) then
	   grochk = 1

      else
	 grochk = 0
      endif

      return
      end
	 
