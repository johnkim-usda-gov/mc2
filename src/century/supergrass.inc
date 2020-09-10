c*****************************************************************
c	$RCSfile: supergrass.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: supergrass.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:49:01  conklind
cInitial revision
c
c******************************************************************
c      This is the supertree set of evergreen/deciduous variables.
c
c
c-----prdx1(1) and prdx2(1) - c3/c4 max production (in supertree.inc)
c-----ppdf1(1-4) and ppdf2(1-4) - c3/c4 temperature curve parameters
c---------------------------------(in supertree.inc)
c-----pramn1 and pramn2 - min/max c/n, biomass = 0/max
c-----pramx1 and pramx2 - min/max c/n, biomass = 0/max
c-----prbmn1 and prbmn2 - intercept c/n as function of ppt
c-----prbmx1 and prbmx2 - slope c/n as function of ppt


      common /suprg/ pramn1(3,2), pramx1(3,2), prbmn1(3,2), prbmx1(3,2),
     .               pramn2(3,2), pramx2(3,2), prbmn2(3,2), prbmx2(3,2),
     .               suprgEnd

      real  pramn1, pramx1, prbmn1, prbmx1, pramn2, pramx2, prbmn2, 
     .               prbmx2        

      real suprgStart, suprgEnd
      equivalence (suprgStart, pramn1(1,1))
      save /suprg/
