c*****************************************************************
c	$RCSfile: schvar.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: schvar.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 18:55:13  conklind
cInitial revision
c
c******************************************************************
c23456789012345678921234567893123456789412345678951234567896123456789712
c**** Scheduling option variables
c     evtptr	array index of next event in cmdary, timary, typary
c     fltary    array holding real (float) values for each event
c     rptyrs    number of years actually included in the block
c     timary	array holding year and month for each event in the block
c     ttlind    total number of events read in a block

      common/schvar/ evtptr, fltary(2500,6), rptyrs, timary(2500,2),
     -               ttlind, code_flags(10), runoff_x_intercept, 
     -               runoff_slope, lait_low_lim,
     -	             schvarEnd

      integer evtptr, rptyrs, timary, ttlind
      real fltary
      logical code_flags
      real runoff_x_intercept, runoff_slope, lait_low_lim
           
      real schvarStart, schvarEnd
      equivalence (schvarStart, evtptr)
      save /schvar/
