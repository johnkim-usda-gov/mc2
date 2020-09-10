
c*****************************************************************
c	$RCSfile: lacalc2.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: lacalc2.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:51:18  bachelet
cheader
c
c******************************************************************
C...LACALC2.F
c...CALLED FROM:  potfor
c...LENIHAN 9/10/97


	subroutine lacalc2(lait, rlwodc, maxlai, klai)

	real      lait, rlwodc, maxlai, klai

        include 'vemap2.inc'


c...Local variables

        real rlai, tlai

c...Calculate tree leaf lai

        rlai = lait

        tlai = maxlai * rlwodc / (klai + rlwodc)

c...Choose the LAI reducer on production.  I don't really understand
c     why we take the average in the first case, but it will probably
c     change...

        if (rlai .lt. tlai) then
           lait = (rlai + tlai) / 2.0
        else
           lait = tlai
        endif

        mapss_lai = lait

        return
        end
