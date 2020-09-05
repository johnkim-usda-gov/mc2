	
c*****************************************************************
c	$RCSfile: soilos.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: soilos.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:46:03  bachelet
cheader
c
c******************************************************************
	subroutine soilos(time,nelem,nlr,flost,somc,somci,
     $ 	                  csrsnk,some,esrsnk)

	include 'const.inc'

	integer   nelem, nlr
	real      time, flost, somc, somci(nlr,ISOS), csrsnk(ISOS), 
     $            some(nlr,MAXIEL), esrsnk(MAXIEL)

c...Compute soil loss for som1, som2, or som3.
c...nlr tells the number of layers being simulated for this som
c     component.  It will be 2 for som1; it will be 1 for som2
c     and som3.  Only the soil layer is considered here, so nlr
c     can be used as the index.   vek  08-91

c...Local variables
	integer   ielem
	real      accum(ISOS), closs, eloss


      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0


c...Loss of carbon isotopes
	closs = somc * flost
	call csched(closs,somci(nlr,LABELD),somc,
     $            somci(nlr,UNLABL),csrsnk(UNLABL),
     $            somci(nlr,LABELD),csrsnk(LABELD),
     $            1.0,accum)

c...Loss for each other element is based on element/carbon ratio
	do 10 ielem = 1, nelem
  	  eloss = closs * some(nlr,ielem)/somc
	  call flow(some(nlr,ielem),esrsnk(ielem),time,eloss)
10	continue

	return
	end
