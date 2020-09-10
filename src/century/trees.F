C*****************************************************************
c       $RCSfile: trees.F,v $
c       $Revision: 1.1 $
c       $Date: 2007/01/04 22:26:56 $
c       $Locker:  $
c
c$Log: trees.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.4  2004/12/13 21:21:09  bachelet
ccommented out equation for ratio, deleted from declaration stmt
c
cRevision 1.3  2004/09/24 23:51:05  bachelet
cheader/clean
c
cRevision 1.2  2004/05/19 00:08:12  bachelet
cupdating rimpct to be consistent with growth.F
c
cRevision 1.1  2004/05/18 23:51:41  bachelet
cInitial revision
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

      subroutine trees (month, cisotf, wfunc)

      integer   month
      real      cisotf, wfunc

c...Simulate forest production for the month.

      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'zztim.inc'

      include 'mapssvar.inc'

c...Function declarations
      real      fsfunc, rtimp
      INTEGER	GROCHk
      external  fsfunc, rtimp, grochk

c...Local variables
      integer   iel, ipart, iptr, lyr
      real      amt, calcup, cfrac(FPARTS),
     $          euf(FPARTS), fsol, rimpct,
     $          remc, totcup,
     $          uptake(4,MAXIEL), avefrc(MAXIEL)

c...Death of tree parts
      call wdeath(tave, wfunc)

c...Option to add organic matter; added here so that it won't be
c     done twice in savanna
      if (doomad) then
        call partit(astgc,astrec,1,csrsnk,esrsnk,astlig,astlbl)
      endif

c...Update flows so direct absorption will be accounted for
c     before plant uptake

      call flowup(time)
      call sumcar

c...Calculate impact of root biomass on available nutrients
      rimpct = rtimp(riint, rictrl, frootc)

c...Determine old or new forest
c     ptr points to new forest carbon allocation fractions (ptr = 1) or
c     mature forest carbon allocation fractions (ptr = 2) for each of the
c     tree parts; leaves, fine roots, fine branches, large wood, and
c     coarse roots.  switch from new forest allocation fractions to old
c     forest allocation fractions at time = swold
 
      if (time .le. swold) then
        iptr = 1
      else
        iptr = 2
      endif

c...Calculate carbon fraction in each part
      do 30 ipart = 1, FPARTS
	cfrac(ipart) = fcfrac(ipart, iptr)
30    continue

c...For seasonal deciduous/conifer forest system, reapportion growth
c     if its the first month of the growing season to simulate leaf
c     bloom (80% of growth goes to leaves)

      if (decid .eq. 1) then
c...Check to see if it is the first month of the growing season
	if (grochk(tave) .eq. 1) then

c...Carbon reapportionment . . .
c       Changed to balance the portions  -rm 5/91
c       REMC - the remaining carbon to be partitioned
c       TOTCUP - total of carbon uptakes for eveything but leaves

c---B. McKeown 10-30-97 New carbon reapportionment to leaves
c     Leaves get 80% of the carbon until they achieve the 
c     maximum LAI possible for the amount of large wood.

          if (lait .lt. rlwodc / ((klai+ rlwodc) / (maxlai))) then
            cfrac(LEAF) = cfrac(LEAF) + wooddr(LEAF) * 
     .                  (0.8 - cfrac(LEAF))

            remc = 1.0 - cfrac(LEAF)
            totcup = 0
            do 40 ipart = FROOT, CROOT
              totcup = cfrac(ipart) + totcup
40          continue

            do 50 ipart = FROOT, CROOT
              cfrac(ipart) = remc * (fcfrac(ipart,iptr) / totcup)
50          continue

          endif
        endif
      endif

c...Previous flowup should have updated mineral pools
        do 60 iel = 1, nelem
	  avefrc(iel) = 0.0
          do 55 lyr = 1, nlayer
          if (minerl(lyr, iel) .gt. 0.) then
	          avefrc(iel) = avefrc(iel) + minerl(lyr,iel)
             endif
55        continue
60      continue

c...Determine actual production values, restricting the C/E ratios

      if (forgrw .eq. 1 .and. pforc.gt.0.) then
	  call restrp(nelem,FPARTS,avefrc,ccefor,cfrac,pforc,rimpct,
     $                forstg,snfxmx(FORSYS),cprodf,eprodf,eup,
     $                uptake,elimit,nfix,relyld,snfxac(FORSYS),svuptk)
      else
	  cprodf = 0.
      endif

c...If growth occurs...
      if (cprodf .gt. 0.) then

c...Update accumulators for N, P, and S uptake
	do 80 iel = 1, nelem
	  eupacc(iel) = eupacc(iel) + eprodf(iel)
	  do 70 ipart = 1, FPARTS
   	    eupprt(ipart,iel) = eupprt(ipart,iel) + eup(ipart,iel)
70        continue
80	continue

c...C/N ratio for production
c	tcnpro = cprodf/eprodf(1)

c...Calculate production for each tree part
	do 85 ipart = 1, FPARTS
	  cfrac(ipart) = cfrac(ipart) * cprodf
	  if (cfrac(ipart) .ne. cfrac(ipart)) then
	    print *, 'trees: cfrac(ipart) is a nan: ipart = ', ipart
		stop
	  endif
85	continue

c...Growth of leaves; split into labeled & unlabeled parts
        call csched(cfrac(LEAF),cisotf,1.0,
     $              csrsnk(UNLABL),rlvcis(UNLABL),
     $              csrsnk(LABELD),rlvcis(LABELD),
     $              1.0,alvcis)

c...Growth of fine roots; split into labeled & unlabeled parts
        call csched(cfrac(FROOT),cisotf,1.0,
     $              csrsnk(UNLABL),frtcis(UNLABL),
     $              csrsnk(LABELD),frtcis(LABELD),
     $              1.0,afrcis)

c...Growth of fine branches; split into labeled & unlabeled parts
        call csched(cfrac(FBRCH),cisotf,1.0,
     $              csrsnk(UNLABL),fbrcis(UNLABL),
     $              csrsnk(LABELD),fbrcis(LABELD),
     $              1.0,afbcis)

c...Growth of large wood; split into labeled & unlabeled parts
        call csched(cfrac(LWOOD),cisotf,1.0,
     $              csrsnk(UNLABL),rlwcis(UNLABL),
     $              csrsnk(LABELD),rlwcis(LABELD),
     $              1.0,alwcis)

c...Growth of coarse roots; split into labeled & unlabeled parts
        call csched(cfrac(CROOT),cisotf,1.0,
     $              csrsnk(UNLABL),crtcis(UNLABL),
     $              csrsnk(LABELD),crtcis(LABELD),
     $              1.0,acrcis)

c...Actual Uptake
	do 90 iel = 1, nelem
          euf(LEAF) = eup(LEAF,iel) / eprodf(iel)
	  euf(FROOT) = eup(FROOT,iel) / eprodf(iel)
	  euf(FBRCH) = eup(FBRCH,iel) / eprodf(iel)
	  euf(LWOOD) = eup(LWOOD,iel) / eprodf(iel)
	  euf(CROOT) = eup(CROOT,iel) / eprodf(iel)

c...Take up nutrients from internal storage pool
	  amt = uptake(ESTOR,iel) * euf(LEAF)
          call flow(forstg(iel),rleave(iel),time,amt)
	  amt = uptake(ESTOR,iel) * euf(FROOT)
          call flow(forstg(iel),froote(iel),time,amt)
	  amt = uptake(ESTOR,iel) * euf(FBRCH)
          call flow(forstg(iel),fbrche(iel),time,amt)
	  amt = uptake(ESTOR,iel) * euf(LWOOD)
          call flow(forstg(iel),rlwode(iel),time,amt)
	  amt = uptake(ESTOR,iel) * euf(CROOT)
          call flow(forstg(iel),croote(iel),time,amt)


c...Take up nutrients from soil
          do 100 lyr = 1, nlayer
            if (minerl(lyr,iel) .gt. 0.) then
              fsol = 1.0
              if (iel .eq. P) then
                fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
              endif
              calcup =
     $   uptake(ESOIL,iel)*minerl(lyr,iel)*fsol/avefrc(iel)

	      amt = calcup * euf(LEAF)
              call flow(minerl(lyr,iel),rleave(iel),time,amt)
	      amt = calcup * euf(FROOT)
              call flow(minerl(lyr,iel),froote(iel),time,amt)
	      amt = calcup * euf(FBRCH)
              call flow(minerl(lyr,iel),fbrche(iel),time,amt)
	      amt = calcup * euf(LWOOD)
              call flow(minerl(lyr,iel),rlwode(iel),time,amt)
	      amt = calcup * euf(CROOT)
              call flow(minerl(lyr,iel),croote(iel),time,amt)
             endif
100        continue


c...Take up nutrients from nitrogen fixation
          if (iel .eq. 1 .and. nfix .gt. 0.) then
            amt = uptake(ENFIX,iel) * euf(LEAF)
            call flow(esrsnk(iel),rleave(iel),time,amt)
            amt = uptake(ENFIX,iel) * euf(FROOT)
            call flow(esrsnk(iel),froote(iel),time,amt)
            amt = uptake(ENFIX,iel) * euf(FBRCH)
            call flow(esrsnk(iel),fbrche(iel),time,amt)
            amt = uptake(ENFIX,iel) * euf(LWOOD)
            call flow(esrsnk(iel),rlwode(iel),time,amt)
            amt = uptake(ENFIX,iel) * euf(CROOT)
            call flow(esrsnk(iel),croote(iel),time,amt)
          endif
 
c...Take up nutrients from fertilizer
          if (aufert .ne. 0) then
            if (uptake(EFERT,iel) .gt. 0.) then

c...Automatic fertilizer added to plant pools
              amt = uptake(EFERT,iel) * euf(LEAF)
              call flow(esrsnk(iel),rleave(iel),time,amt)
              amt = uptake(EFERT,iel) * euf(FROOT)
              call flow(esrsnk(iel),froote(iel),time,amt)
              amt = uptake(EFERT,iel) * euf(FBRCH)
              call flow(esrsnk(iel),fbrche(iel),time,amt)
              amt = uptake(EFERT,iel) * euf(LWOOD)
              call flow(esrsnk(iel),rlwode(iel),time,amt)
              amt = uptake(EFERT,iel) * euf(CROOT)
              call flow(esrsnk(iel),croote(iel),time,amt)

c...Automatic fertilizer added to mineral pool
              amt = uptake(EFERT,iel) * (1./favail(iel) - 1.)
              fertot(iel) = fertot(iel) + uptake(EFERT,iel) + amt
              call flow(esrsnk(iel),minerl(SRFC,iel),time,amt)
            endif
	  endif
90      continue
 
c...Else there is no production this month
      else  
	cprodf = 0.
	do 140 iel = 1, 3
	   eprodf(iel) = 0.
	   do 130 ipart = 1, FPARTS
             eup(ipart,iel) = 0.
130	   continue
140	continue
      endif

c...Forest removal option
c...Moved to SIMSOM for consistency when running savanna model (rm 7/94)

c     if (dotrem .or. dofire(FORSYS)) then
c       call frem()
c     endif

      return
      end
