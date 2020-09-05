c*****************************************************************
c	$RCSfile: seq.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: seq.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:48:48  conklind
cInitial revision
c
c******************************************************************
      common/seq/cursys, decsys,
     -  seqEnd

      integer cursys, decsys

c...cursys = 1 grass/crop
c          = 2 forest
c          = 3 savanna

c...decsys = 1 No wood is present to decompose
c          = 2 Wood is present to decompose

      real seqStart, seqEnd
      equivalence (seqStart, cursys)
      save /seq/
