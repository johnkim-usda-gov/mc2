
c*****************************************************************
c	$RCSfile: init_climate.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: init_climate.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:46:19  bachelet
cheader
c
c******************************************************************

      subroutine cen_init_climate(ppt, tmax, tmin)      

      include 'wth.inc'
      include 'param.inc'

      real ppt(12), tmax(12), tmin(12)
      integer m

      do 50 m = 1, 12
        precip(m) = ppt(m) * 0.1
		if (precip(m) .ne. precip(m)) then
		  precip(m) = ppt(m) * 0.1
	      if (precip(m) .ne. precip(m)) then
		    print *, 
     $          '*** cen_init_climate(): m, precip(m), ppt(m) = ', 
     $          m, precip(m), ppt(m)
            stop
          endif
		  print *, '*** cen_init_climate(): ',
     $        'NaN in precip(m) has been fixed. m, precip(m) = ', 
     $        m, precip(m)
	    endif		
        tmx2m(m) = tmax(m)
        tmn2m(m) = tmin(m)
c-------More climate arrays to satisfy...
        mintmp(m) = tmin(m)
        maxtmp(m) = tmax(m)
        prcurr(m) = ppt(m) * 0.1
        prcnxt(m) = ppt(m) * 0.1
		if (ppt(m) .ne. ppt(m)) then
		  print *, 'cen_init_climate: ppt(m) is a NaN. m = ', m
		  stop
		endif
50    continue
      return
      end
