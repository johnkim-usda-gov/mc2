c*****************************************************************
c	$RCSfile: chrvar.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: chrvar.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:38:40  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:42:32  conklind
cInitial revision
c
c******************************************************************
c**** Character values needed in common block:
c     cmdary    array holding command for each event in the block
c     curcrp    the current crop
c     curtre   the current tree
c     initcp	the initial crop
c     initre	the initial tree
c     sitnam	name of site file
c     typary    array holding command type for each event in block
c     wlabel	string of weather labels "prec", "tmin", "tmax"
c     wthnam    if wthr = F, the weather file name
c     wthr      the weather choice

      common/chrvar/cmdary(2500), curcrp, curtre, initcp, initre, 
     $              sitnam,typary(2500), wlabel(3), wthnam, wthr,
     $              chrvarPad2, chrvarEnd

      character*20 cmdary, wthnam, wthr
      character*5 curcrp, curtre, initcp, initre, typary
      character*50 sitnam
      character*4 wlabel

      real chrvarStart, chrvarEnd
      character*2 chrvarPad2
      equivalence (chrvarStart, cmdary(1))
      save /chrvar/
