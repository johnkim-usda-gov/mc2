c*****************************************************************
c	$RCSfile: supertree.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: supertree.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 19:58:42  conklind
cInitial revision
c
c******************************************************************
c      This is the supertree set of evergreen/deciduous variables.

c------C. Daly 1-25-98.  Expanded to four classes: DN, EN, DB, EB

c------C. Daly 2-17-98.  Added wet and dry classes.

c------C. Daly 7-16-98.  Added wet and dry needleleaf index.
c
c
c-----prdx1(3) and prdx2(3) - ever/decid max production
c-----ppdf1(1-4) and ppdf2(1-4) - ever/decid temperature curve parameters
c-----cerfor1 and cerfor2(1,1,1), (1,2,1), (1,3,1), (1,4,1), (1,5,1)
c-------------------- - ever/decid min C/N ratios
c-----cerfor1 and cerfor2(2,1,1), (2,2,1), (2,3,1), (2,4,1), (2,5,1)
c-------------------- - ever/decid max C/N ratios
c-----fcfrac1 and fcfrac2(1,2),(2,2),(3,2),(4,2),(5,2) - ever/decid C alloc
c-----leafdr1 and leafdr2(1-12) - ever/decid monthly leaf death rate
c-----btolai1 and btolai2 - ever/decid SLA
c-----maxlai1 and maxlai2 - ever/decid max potential LAI
c-----wooddr1(1-5) and wooddr2(1-5)- ever/decid monthly death rates
c-----sitpot1 and sitpot2 - ever/decid N given to grasses
c-----kali1 and klai2 - 

      common /suprt/ prdx1(3),prdx2(3),prdx3(3),prdx4(3),prdxw(3),
     .          prdxd(3),ppdf1(4,2),ppdf2(4,2),ppdf3(4,2),
     .          ppdf4(4,2),ppdfw(4,2),ppdfd(4,2),
     .          cerfor1(3,5,3),cerfor2(3,5,3),cerfor3(3,5,3), 
     .          cerfor4(3,5,3),cerforw(3,5,3),cerford(3,5,3),
     .          fcfrac1(5,2),fcfrac2(5,2),fcfrac3(5,2),fcfrac4(5,2),
     .          fcfracw(5,2),fcfracd(5,2),
     .          leafdr1(12),leafdr2(12),leafdr3(12),leafdr4(12),
     .          leafdrw(12),leafdrd(12), 
     .          btolai1,btolai2,btolai3,btolai4,btolaiw,btolaid,
     .          maxlai1,maxlai2,maxlai3,maxlai4,maxlaiw,maxlaid,
     .          wooddr1(5),wooddr2(5),wooddr3(5),wooddr4(5), 
     .          wooddrw(5),wooddrd(5),
     .          sitpot1,sitpot2,sitpot3,sitpot4,sitpotw,sitpotd,
     .          klai1,klai2,klai3,klai4,klaiw,klaid,
     .          needled,needlew,
     .          suprtEnd

      real      prdx1,prdx2,prdx3,prdx4,prdxw,prdxd,
     .          ppdf1,ppdf2,ppdf3,ppdf4,ppdfw,ppdfd,
     .          cerfor1,cerfor2,cerfor3,cerfor4,cerforw,cerford,
     .          fcfrac1,fcfrac2,fcfrac3,fcfrac4,fcfracw,fcfracd,
     .          leafdr1,leafdr2,leafdr3,leafdr4,leafdrw,leafdrd, 
     .          btolai1,btolai2,btolai3,btolai4,btolaiw,btolaid,
     .          maxlai1,maxlai2,maxlai3,maxlai4,maxlaiw,maxlaid,
     .          wooddr1,wooddr2,wooddr3,wooddr4,wooddrw,wooddrd, 
     .          sitpot1,sitpot2,sitpot3,sitpot4,sitpotw,sitpotd,
     .          klai1,klai2,klai3,klai4,klaiw,klaid,
     .          needled,needlew


      real suprtStart, suprtEnd
      equivalence (suprtStart, prdx1(1))
      save /suprt/
