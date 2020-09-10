c*****************************************************************
c	$RCSfile: timvar.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: timvar.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:49:21  conklind
cInitial revision
c
c******************************************************************
c****  Time variables not declared by Time-0
c      blktnd = ending time of the block
c      decodt = time step used for the decomposition model
c      month = current month
c      strtyr = starting year of the block
c      tplt = next time (years) when variables should be written to the
c             plot file

       common/timvar/blktnd, decodt, month, strtyr, tplt,
     -   timvarEnd
     
       integer month, strtyr
       real blktnd, decodt, tplt

      real timvarStart, timvarEnd
      equivalence (timvarStart, blktnd)
      save /timvar/
