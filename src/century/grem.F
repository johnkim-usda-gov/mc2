
c*****************************************************************
c	$RCSfile: grem.F,v $
c	$Revision: 1.2 $
c	$Date: 2004/09/24 20:32:30 $
c	$Locker:  $
c	
c$Log: grem.F,v $
cRevision 1.2  2004/09/24 20:32:30  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine grem()

c...Simulate removal of crop/grass by fire or grazing for the month.

      include 'const.inc'
      include 'dovars.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'
      include 'seq.inc'

c...Local variables
      integer  i, iel, lyr
      real     flrem, fdrem(2), fcret, shremc, shreme(MAXIEL),
     $         sheret, litrme(MAXIEL), sdremc, sdreme(MAXIEL),
     $         sderet, cret, eret(MAXIEL), feces, urine,
     $         recres(MAXIEL), ciso, friso

C_NOTES: 
c...ciso tells what fraction of the C returned is labeled (C14).
c...Initialize flrem, fdrem, and fcret based on fire or grazing.
c...Mod. removed FCRET=FFCRET, FCRET=0 as initialized. -rm 8-90.
c...Mod. fire routines, created 'litburn() subroutine  -mse 4-94.
c...If dofire(cursys)=1 -> CRPSYS burns standing dead and litter.
 
      cret = 0.0

c-----LENIHAN 4/28/97

      if ( did_burn ) then
        flfrem    = mf_flfrem
        fdfrem(1) = mf_fdfrem(1)
        fdfrem(2) = mf_fdfrem(2)
      endif
  
      if (dofire(CRPSYS).or.dofire(SAVSYS)) then
        flrem = flfrem
        fdrem(1) = fdfrem(1)
        fdrem(2) = fdfrem(2)
c       fcret = ffcret
      else if (dograz) then
        flrem = flgrem
        fdrem(1) = fdgrem
        fcret = gfcret
      endif

 
c...Added for local initialization of variables which may not
c     get initialized during a run. 8-31-90 -rm
 
      ciso = 0
      shremc = 0.0
      sdremc = 0.0
      do 10 i = 1, MAXIEL
        sdreme(i) = 0.0
        shreme(i) = 0.0
        litrme(i) = 0.0
 10   continue
 

c...Shoots removed
c...increased the lower limits to avaoid NaN when pools get low DB 11-13-98
      if (aglivc .gt. 0.0) then

c...carbon
        shremc = flrem*aglivc
        shrema = shrema + shremc
        call csched(shremc,aglcis(LABELD),aglivc,
     $              aglcis(UNLABL),csrsnk(UNLABL),
     $              aglcis(LABELD),csrsnk(LABELD),
     $              1.0,shrmai)
        ciso = ciso + (shremc*aglcis(LABELD)/aglivc)

c...elements
        do 20 iel = 1, nelem
          shreme(iel) = shremc*aglive(iel)/aglivc
          shrmae(iel) = shrmae(iel) + shreme(iel)
          call flow(aglive(iel),esrsnk(iel),time,shreme(iel))
20      continue
      endif

c...Standing dead removed
      if (stdedc .gt. 0.) then

c...carbon
        sdremc = fdrem(1)*(stdedc)
        sdrema = sdrema + sdremc
        call csched(sdremc,stdcis(LABELD),stdedc,
     $              stdcis(UNLABL),csrsnk(UNLABL),
     $              stdcis(LABELD),csrsnk(LABELD),
     $              1.0,sdrmai)
        ciso = ciso + (sdremc*stdcis(LABELD)/stdedc)

c...elements
        do 30 iel = 1, nelem
          sdreme(iel) = sdremc*stdede(iel)/stdedc
          sdrmae(iel) = sdrmae(iel) + sdreme(iel)
          call flow(stdede(iel),esrsnk(iel),time,sdreme(iel))
30      continue
      endif

c...FIRE
      if (dofire(CRPSYS).or. dofire(SAVSYS)) then

c...Residue (surface litter) removed by fire       vek 5/26/90

        call litburn(litrme)

        ciso = ciso + (fdrem(2)*strcis(SRFC,LABELD))
        ciso = ciso + (fdrem(2)*metcis(SRFC,LABELD))



c...Carbon and nutrient return following removal values for
c     the following variables set in greopt :
c         fcret - fraction of carbon returned 
c         fret(iel) - fraction of element returned by fire
c         gret(iel) - fraction of element returned by grazing
c
c...The following variables have units g/m**2/month and are:
c         cret - C returned to system
c         sheret - elemental return for shoots
c         sderet - elemental return for standing dead and litter
c         eret(iel) - total elemental return for aboveground removal
c         urine - amount of urine returned
c         feces- amount of fecal material returned (N, P, S)
c
c...Assume no Carbon return from fire.  If the next 2 statements are
c     reactivated, it will be necessary to partition the flow into
c     strcis(1) and strcis(2) instead of strucc(1):
c         cret = fcret * (shremc + sdremc)
c         call flow(csrsnk(1),strucc(1),time,cret)

        do 60 iel = 1, nelem
          sheret = fret(iel) * shreme(iel)
          sdreme(iel) = sdreme(iel) + litrme(iel)
          sderet = fret(iel) * sdreme(iel)
          eret(iel) = sheret + sderet
          call flow(esrsnk(iel),minerl(1,iel),time,eret(iel))
60      continue

c...END FIRE
 
c...GRAZE
      else
 
C...NOTES:
c...Grazing return with feces and urine explicitly separated.
c...All carbon returned by grazing is in the form of feces.
c...cret is the amount returned.
c...To adjust for the changing lignin content of added material
c     strucc(1) and strlig are recomputed.

        cret = fcret * (shremc + sdremc)
        if (cret .le. 0.0) then
          cret = 0.0
          do 70 iel = 1, nelem
            eret(iel) = 0.0
70        continue
        else
          do 80 iel = 1, nelem
            sheret = gret(iel) * shreme(iel)
            sderet = gret(iel) * sdreme(iel)
            eret(iel) = sheret + sderet
            urine= (1-fecf(iel)) * eret(iel)
            feces= fecf(iel) * eret(iel)
            recres(iel) = feces/cret
            call flow(esrsnk(iel),minerl(1,iel),time,urine)
80        continue
          lyr = 1

c...Mod. to add structural & metabolic C into labeled (numerator)
c     and total (denominator) C removed.  (vek  05/26/90)
c...friso tells what fraction of the C returned is labeled

          friso = ciso / (shremc + sdremc)
          call partit(cret,recres,lyr,csrsnk,esrsnk,feclig,friso)
        endif

c**** END GRAZE
      endif

c**** Accumulate amounts returned
      creta= creta + cret
      do 90 iel = 1, nelem
        ereta(iel)= ereta(iel) + eret(iel)
90    continue

      return
      end
