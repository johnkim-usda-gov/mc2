
c*****************************************************************
c	$RCSfile: tcalc.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: tcalc.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:48:47  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


	real function tcalc(stemp, teff)

	real stemp, teff(3), result

c...This function computes the effect of temperature on
c   decomposition.  It is an exponential function.  Older
c   versions of Century used a density function.
c   Created 10/95 - rm
c
c...Called From:	CYCLE
c
c...Variables
c	STEMP:		Soil Temperature
c	TEFF(1):	Intercept
c	TEFF(2):	Slope
c	TEFF(3):	Exponent

	result = teff(1) + teff(2) * exp(teff(3) * stemp)
	if (result .ne. result) then
	  print *, 'tcalc is a nan. stemp, teff(1:3) = ',
     $    stemp, teff(1), teff(2), teff(3)
	  stop
	endif
	
	tcalc = result	
	return
	end

