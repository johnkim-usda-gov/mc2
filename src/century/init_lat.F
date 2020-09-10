
c*****************************************************************
c	$RCSfile: init_lat.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: init_lat.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:46:46  bachelet
cheader
c
c******************************************************************

      subroutine cen_init_lat(lat,elev)      

      include 'site.inc'

      real lat,elev

      sitlat = abs(lat)
      sitelev = elev

      return
      end
