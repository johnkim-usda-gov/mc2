      
c*****************************************************************
c	$RCSfile: wthini.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: wthini.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:56:04  bachelet
cheader
c
c******************************************************************
	subroutine wthini(precip,prcstd,prcskw,mintmp,maxtmp)

      real      precip(12), prcstd(12), prcskw(12), mintmp(12), 
     $          maxtmp(12)

c...Determine what weather data will be used.
 
      include 'chrvar.inc'
      include 'const.inc'
      include 'wth.inc'


c...Initialize mintmp and maxtmp
      do 10 m = 1, MONTHS
        mintmp(m) = tmn2m(m)
        maxtmp(m) = tmx2m(m)
10    continue

c...Use values in precip
        do 20 m = 1, MONTHS
          prcnxt(m) = precip(m)
20      continue


      return
      end
