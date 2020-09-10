
c*****************************************************************
c	$RCSfile: gpdf.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: gpdf.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:31:31  bachelet
cheader
c
c******************************************************************
C...GPDF.F

        real function gpdf(x,a,b,c,d) 

	real x, a, b, c, d


c******************** flowlib ********************* 

c...(run-time sub-set of modaid, exclusive of modctl) 

c...Release 1.0  (first formal release of modaid) 

c       james m. vevea
c       natural resource ecology lab
c       colorado state university 
c       fort collins, colorado  80523 

c...This routine is functionally equivalent to the routine of the
c     same name, described in the publication: 

c       Some Graphs and their Functional Forms 
c       Technical Report No. 153 
c       William Parton and George Innis (1972) 
c       Natural Resource Ecology Lab.
c       Colorado State University
c       Fort collins, Colorado  80523

c...12/90 Corrected by McKeown - exponent on frac changed from d to c

 
c...Local variables
        integer, parameter :: dp = selected_real_kind(15,307)
        real(dp)    frac

        frac = (b-x) / (b-a) 
        gpdf = 0.
        if (frac .gt. 0.) then
          gpdf = exp(c/d * (1. - frac**d)) * (frac**c)
        endif

        return 
        end
