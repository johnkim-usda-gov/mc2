c*****************************************************************
c	$RCSfile: dovars.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: dovars.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:45:00  conklind
cInitial revision
c
c******************************************************************
c23456789012345678921234567893123456789412345678951234567896123456789712
c**** Scheduling option variables
c     docult	time for a cultivation event
c     doerod	time for an erosion event
c     dofert	time for a fertilization event
c     dofire	time for a fire event
c                          (3) - grass/crop,forest,savannah
c     doflst	time for a forest last month growth event
c     dofone	time for a forest first month growth event
c     dofrst	time for a crop first growth event
c     dograz	time for a grazing event
c     dohrvt	time for a harvest event
c     doirri	time for an irrigation event
c     dolast	time for a crop last growth event
c     doomad	time for a organic matter addition event
c     doplnt	time for a planting event
c     dosene	time for a senescence event
c     dotrem	time for a forest removal event

      common/dovars/docult, doerod, dofert, dofire(3), doflst, dofone,
     -              dofrst, dograz, dohrvt, doirri, dolast, doomad,
     -              doplnt, dosene, dotrem,
     -              dovarsEnd
      logical docult, doerod, dofert, dofire, doflst, dofone
     - ,dofrst, dograz, dohrvt, doirri, dolast, doomad
     - ,doplnt, dosene, dotrem

      real dovarsStart, dovarsEnd
      equivalence (dovarsStart, docult)
      save /dovars/
