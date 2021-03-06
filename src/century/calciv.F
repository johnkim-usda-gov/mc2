
c*****************************************************************
c	$RCSfile: calciv.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: calciv.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:00:21  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


        subroutine calciv

c...Calculate initial values for temperature, water and live root
c     carbon variables.
c...Called from detiv.
c...Note that variables which are functions of changeable parameters
c     (i.e. read from 'site'.par) should be computed in prelim instead
c     of calciv.

	include 'chrvar.inc'
	include 'const.inc'
	include 'ligvar.inc'
	include 'param.inc'
	include 'parcp.inc'
	include 'parfs.inc'
	include 'plot1.inc'
	include 'plot3.inc'
	include 'site.inc'
	include 'zztim.inc'

c...Local variables
        integer   iel, ilayer, iso, m
        real      avtemp, arain, dumye(MAXIEL), dumyc(ISOS), frc14,
     $            reclt1(MAXIEL), reclt2(MAXIEL), tcg, tcl
	character string*80, char1*1


c...Initialize soil C pools using Burke's equations.
c       ivauto = 0	the user has supplied the initial values
c       ivauto = 1	initialize using the grassland soil parameters
c       ivauto = 2	initialize using the crop soil parameters

c...Initialize dumyc and dumye variables.
	dumyc(LABELD) = 1000.
	dumyc(UNLABL) = 1000.
	dumye(N) = 100.
	dumye(P) = 100.
	dumye(S) = 100.

c...Compute mean annual temperature (avtemp) and mean annual 
c     precipitation (arain)
	avtemp = 0.
	arain = 0.
	do 10 m = 1, MONTHS
	  avtemp = avtemp + (mintmp(m) + maxtmp(m))/2.
	  arain = arain + precip(m)
10	continue
	avtemp = avtemp/12.
	if (avtemp .gt. 23.) then
          avtemp = 23.
        endif
	if (arain .gt. 120.) then
          arain = 120.
        endif

c...Initialize soil C pools for a grassland soil
        if (ivauto .eq. 1) then

c...tcg = total soil carbon in grams (som1c + som2c(1) + som3c(1))
          tcg = (-8.27E-01 * avtemp + 2.24E-02 * avtemp * avtemp
     $		      + arain * 1.27E-01 - 9.38E-04 * arain * arain
     $		      + arain * silt * 8.99E-02
     $		      + arain * clay * 6.00E-02 + 4.09) *1000.

c...Assign a fixed value to surface som1.   vek  08-91
          som1ci(SRFC,UNLABL) = 10.
          som1ci(SRFC,LABELD) = 0.

c...Burke's equations only apply to soil compartments.
          som1ci(SOIL,UNLABL) = tcg * .02
          som1ci(SOIL,LABELD) = 0.
	  som2ci(UNLABL) = tcg * .64
	  som2ci(LABELD) = 0.
	  som3ci(UNLABL) = tcg * .34
	  som3ci(LABELD) = 0.
	  stdcis(UNLABL) = 80.
	  stdcis(LABELD) = 0.
	  stdede(N) = 1.6
	  stdede(P) = .3
	  stdede(S) = .3
	  bglcis(UNLABL) = 200.
	  bglcis(LABELD) = 0.
	  bglive(N) = 3.0
	  bglive(P) = .5
	  bglive(S) = .5
	  clittr(SRFC,UNLABL) = 100.
	  clittr(SRFC,LABELD) = 0.
	  clittr(SOIL,UNLABL) = 100.
	  clittr(SOIL,LABELD) = 0.
       endif

c...Initialize soil C pools for cultivated soils
       if (ivauto .eq. 2) then

c...tcg = total soil carbon in grams (som1c + som2c + som3c)
	 tcg = (-7.50E-01 * avtemp + 2.10E-02 * avtemp * avtemp
     $		      +5.81E-02 * arain -4.58E-04 * arain * arain
     $		      +arain * silt * 4.94E-02 
     $		      + arain * 5.82E-02 * clay + 5.15) * 1000.
c...Assign a fixed value to surface som1.   vek  08-91
         som1ci(SRFC, UNLABL) = 10.
         som1ci(SRFC, LABELD) = 0.

c...Burke's equations only apply to soil compartments. vek  08-91
         som1ci(SOIL,UNLABL) = tcg * .02
         som1ci(SOIL,LABELD) = 0.
	 som2ci(UNLABL) = tcg * .54
	 som2ci(LABELD) = 0.
	 som3ci(UNLABL) = tcg * .44
         som3ci(LABELD) = 0.
	 stdcis(UNLABL) = 20.
	 stdcis(LABELD) = 0.
	 stdede(N) = .40
	 stdede(P) = .075
	 stdede(S) = .075
	 clittr(SRFC,UNLABL) = 10.
	 clittr(SRFC,LABELD) = 0.
	 clittr(SOIL,UNLABL) = 10.
	 clittr(SOIL,LABELD) = 0.
      endif

c...End of soil C pool initialization


c...Starting values for nitrogen, phosphorous, and sulfur depend on
c     carbon values and the ratios of carbon to each other element.
c     Initialize structural and metabolic pools C, N, P, and S.
c     First set them to zero and calculate N/C, P/C, & S/C ratios.

      do 40 ilayer = SRFC, SOIL
        do 20 iso = 1, ISOS
	  strcis(ilayer,iso) = 0.
	  metcis(ilayer,iso) = 0.
20      continue
        do 30 iel = 1, MAXIEL
	  struce(ilayer,iel) = 0.
	  metabe(ilayer,iel) = 0.
30      continue
40    continue

c...Compute N/C, P/C, and S/C ratios from C/N, C/P, and C/S.
c...This is for use in partit.
      do 50 iel = 1, MAXIEL
        if (rcelit(SRFC, iel) .gt. 0.) then
          reclt1(iel) = 1. / rcelit(SRFC, iel)
        endif
50    continue
      do 55 iel = 1, MAXIEL
        if (rcelit(SOIL, iel) .gt. 0.) then
          reclt2(iel) = 1. / rcelit(SOIL, iel)
        endif
55    continue

c...Sum carbon isotopes for use in partit.
c	  print *, 'calling sumcar from calciv.F l177'
      call sumcar

c...Split litter C content into structural/metabolic based upon 
c     litter C and litter lignin content and compute structural and
c     metabolic N, P, & S based upon amount of C and the ratios
c     computed above.
      if (initcp .ne. ' ' .and. initre .ne. ' ') then
        pltlig(SRFC) = (wdlig(LEAF)+fligni(INTCPT,ABOVE) +
     +                  fligni(SLOPE,ABOVE) * arain) / 2.0
        pltlig(SOIL) = (wdlig(FROOT)+fligni(INTCPT,BELOW) +
     +                  fligni(SLOPE,BELOW) * arain) / 2.0
      else if (initcp .ne. ' ') then
        pltlig(ABOVE) = fligni(INTCPT,ABOVE)+fligni(SLOPE,ABOVE)*arain
        pltlig(BELOW) = fligni(INTCPT,BELOW)+fligni(SLOPE,BELOW)*arain
      else if (initre .ne. ' ') then
        pltlig(SRFC) = wdlig(LEAF)
        pltlig(SOIL) = wdlig(FROOT)
      endif

c...Total C in litter
      tcl = clittr(SRFC,UNLABL)+clittr(SRFC,LABELD)
      frc14 = clittr(SRFC,LABELD)/tcl
      call partit(tcl,reclt1,1,dumyc,dumye, 
     -                     pltlig(SRFC),frc14)
      tcl = clittr(SOIL,UNLABL)+clittr(SOIL,LABELD)
      frc14 = clittr(SOIL,LABELD)/tcl
      call partit(tcl,reclt2,2,dumyc,dumye, 
     -                     pltlig(SOIL),frc14)

      call flowup(time)
c	  print *, 'calling sumcar from calciv.F l208'
      call sumcar

      do 70 iel=1,MAXIEL
c...Compute N, P, and S for surface and soil som1, as well as for
c     som2 and som3.   vek  08-91
          if (rces1(SRFC,iel) .gt. 0.) then
                  som1e(SRFC,iel)=som1c(SRFC)/rces1(SRFC,iel)
          endif
          if (rces1(SOIL,iel) .gt. 0.) then
                  som1e(SOIL,iel)=som1c(SOIL)/rces1(SOIL,iel)
          endif
          if (rces2(iel) .gt. 0.) then
                  som2e(iel)=som2c(1)/rces2(iel)
          endif
          if (rces3(iel) .gt. 0.) then
                  som3e(iel)=som3c(1)/rces3(iel)
          endif
70    continue

      if (initre .ne. ' ') then
        do 80 iel = 1, MAXIEL
          if (cerfor(IVAL,FBRCH,iel) .gt. 0.) then
                    wood1e(iel)=wood1c/cerfor(IVAL,FBRCH,iel)
          endif
          if (cerfor(IVAL,LWOOD,iel) .gt. 0.) then
                    wood2e(iel)=wood2c/cerfor(IVAL,LWOOD,iel)
          endif
          if (cerfor(IVAL,CROOT,iel) .gt. 0.) then
                    wood3e(iel)=wood3c/cerfor(IVAL,CROOT,iel)
          endif
80      continue
      endif

c...Surface temperature and soil temperature
      tave = (mintmp(1) + maxtmp(1)) / 2.0
      stemp = tave

c...Make sure there is N, P, and S for roots
      if (bglcis(UNLABL)+bglcis(LABELD) .gt. 0.0) then
        do 90 iel = 1, nelem
	  if (bglive(iel) .le. 0.) then
	    char1 = char(ichar(char(iel)) + ichar('0'))
            string = '   Value for bglive(' // char1
     +               // ') must be greater than 0.'
            call message(string)
	    STOP
          endif
90      continue
      endif

c...Initialize grain pools
      cgrain = 0.0
      do 100 iel = 1, MAXIEL
        egrain(iel) = 0.0
100   continue

      return
      end
