c*****************************************************************
c	$RCSfile: parfs.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: parfs.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:48:05  conklind
cInitial revision
c
c******************************************************************
       common/parfs/basfct,btolai,ccefor(2,5,3),cerfor(3,5,3),
     - decid,decw1,decw2,decw3,
     - fcfrac(5,2),forgrw,forrtf(3),klai,laitop,ldrmlt,leafdr(12),
     - lait,maxlai,maxldr,sapk,swold,
     - wdlig(5),wooddr(5),mf_dr(3),
     - parfsEnd

       integer decid,forgrw

       real basfct,btolai,ccefor,cerfor,decw1,decw2,decw3
     -,fcfrac,forrtf,klai,laitop,ldrmlt,leafdr
     -,lait,maxlai,maxldr,sapk,swold
     -,wdlig,wooddr,mf_dr

      real parfsStart, parfsEnd
      equivalence (parfsStart, basfct)
      save /parfs/
