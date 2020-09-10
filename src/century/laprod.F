	
c*****************************************************************
c	$RCSfile: laprod.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: laprod.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:52:03  bachelet
cheader
c
c******************************************************************
	real function laprod(lait,laitop)
	real lait,laitop
      include 'schvar.inc'

c****	This function relates aboveground wood production to leaf area
c	index.  Wood production is in g/m2/mo and then normalized between
c	0 and 1 and then fit to the natural growth function.
c	The natural growth function described in:
c		Some Graphs and Their Functional Forms
c			Parton and Innis
c
c       REF:    Efficiency of Tree Crowns and Stemwood
c               Production at Different Canopy Leaf Densities
c               by Waring, Newman, and Bell
c               Forestry, Vol. 54, No. 2, 1981
c
c	CALLED FROM:  potfor

c****	The minimum LAI to calculate effect is .2.

	if (lait .lt. lait_low_lim) then
          lait = lait_low_lim
        endif
	laprod = 1.0 - exp(laitop * lait)

	return
	end

