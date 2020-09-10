
c*****************************************************************
c	$RCSfile: dshoot.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: dshoot.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.3  2005/01/12 21:20:18  bachelet
cto prevent extremely small pools into aglivc I changed the original if statement
cso that aglivc > 0.0001 rather than just 0 before death can occur.
c
cRevision 1.2  2004/09/24 20:18:31  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine dshoot(wfunc)

      real wfunc

c...Simulate death of shoots for the month.
c
      include 'const.inc'
      include 'dovars.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'
      include 'vemap2.inc'

c...Local variables
      integer   iel, index
      real      accum(ISOS), dthppt, fdeth, sdethe, tostore

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0


c...Death of shoots.  Most shoots die during the month of senescence.
c     During other months, only a small fraction of the shoots die.
c... DB 10/98
      if (aglivc .gt. 0.0001) then
        index = 1
        dthppt = 1. - wfunc
        if (dosene) then
          index = 2
          dthppt = 1.0
        endif
        fdeth = fsdeth(index) * dthppt

c...Increase the death rate of shoots to account for effect of shading.
c     This is not done during senescence (when the death rate is greater
c     than or equal to .4)
        if (fsdeth(index) .lt. 0.4 .and. aglivc .gt. fsdeth(4)) then
      	   fdeth = fdeth + fsdeth(3)
        endif

c...Constrain the fraction
      	if (fdeth .gt. 0.95) then
          fdeth = 0.95
        endif

c...Calculate the amounts and flow
        sdethc = aglivc * fdeth


c-------LENIHAN 9/15/98: Assign sdethc to grass litterfall var 
c-------for VEMAP2 output

        lgrass = sdethc

	call csched(sdethc,aglcis(LABELD),aglivc,
     $              aglcis(UNLABL),stdcis(UNLABL),
     $              aglcis(LABELD),stdcis(LABELD),
     $              1.0,accum)
        do 10 iel = 1, nelem
          sdethe = fdeth * aglive(iel)
          if (iel .eq. N) then
             vlosse = vlossp * sdethe

c... DB 3-16-99 trying to catch underflows and Nans, vlosse gets negative 

	     if(vlosse .gt. 0.0) then
             	call flow(aglive(iel),esrsnk(iel),time,vlosse)
             else
		vlosse = 0.0
             endif	
             volpla = volpla + vlosse
             sdethe = sdethe - vlosse
          endif
          tostore = sdethe * crprtf(iel)
c... adding trap DB 10-14/98
	  if(tostore .gt. aglive(iel)) then
		tostore = aglive(iel)
		crprtf(iel) = tostore / sdethe
          endif
          call flow(aglive(iel),crpstg(iel),time,tostore)
          sdethe = sdethe * (1 - crprtf(iel))
          call flow(aglive(iel),stdede(iel),time,sdethe)
10      continue
      endif

      return
      end
