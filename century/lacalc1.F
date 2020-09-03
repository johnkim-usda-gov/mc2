
c*****************************************************************
c	$RCSfile: lacalc1.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: lacalc1.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:50:46  bachelet
cheader.F
c
c******************************************************************
C...LACALC1.F
c...CALLED FROM:  potfor
c...LENIHAN 9/10/97
c...LENIHAN 7/16/98 Remove fudge factor from lai calc

c---***The constants set in this function should match those set in
c---CenturyFCarbToMapssLaiSLA on the Mapss side. ***

	subroutine lacalc1(lait, rleavc)

	real      lait, rleavc

        include 'plot4.inc'

c...Local variables

        real factor
        real sla, sla_needle, sla_broad

c...Assign values

        factor = 2.2
        sla_needle = 0.003
        sla_broad = 0.012
 
c...Calculate sla and conversion factor based on mix index **/
 
        sla = ((mx_index*0.01)*sla_needle) + 
     .  ((1.0-(mx_index*0.01))*sla_broad)
 
c...Calculate tree leaf lai

        lait = rleavc * factor * sla

        return
        end
