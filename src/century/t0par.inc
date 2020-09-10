c*****************************************************************
c	$RCSfile: t0par.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: t0par.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:49:14  conklind
cInitial revision
c
c******************************************************************
c**** Time variables declared by Time-0; order decided by Time-0
c
c      dt = model time step (currently one month)
c      dtpl = how often (years) the output should be written out
c      tend = simulation end time

      common /t0par/ tend, dtpl, dt,
     -  t0parEnd
           
      real tend, dtpl, dt

      real t0parStart, t0parEnd
      equivalence (t0parStart, tend)
      save /t0par/
