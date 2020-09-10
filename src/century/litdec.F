
c*****************************************************************
c	$RCSfile: litdec.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: litdec.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:53:51  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C...LITDEC.F

      subroutine litdec(dtm)

      real      dtm

c...Litter Decomposition
c     Decompose structural and metabolic material for surface and soil.
c						written by vek 04/91

      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'


c...Function declarations
      real agdrat, bgdrat
      external agdrat, bgdrat

c...Local variables
      integer iel, lyr
      real accum(ISOS), biocnv, cfmes1, co2los, mnrflo, rceto1(3), 
     $     tcflow, ligcon
      logical candec


c...Factor to convert C to biomass is 2.5 for everything but wood.
      parameter (biocnv = 2.5)

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c...Surface STRUCTURAL Material
      if(strucc(SRFC) .gt. 1.e-07) then

c...     Compute total C flow out of structural in layer SRFC
c           print *, 'strucc(SRFC), strmax(SRFC), defac, dec1(SRFC) =',
c     $         strucc(SRFC), strmax(SRFC), defac, dec1(SRFC)
c	       print *, 'pligst(SRFC), strlig(SRFC), dtm = ',
c     $         pligst(SRFC), strlig(SRFC), dtm
           tcflow = min(strucc(SRFC),strmax(SRFC)) * defac * dec1(SRFC) 
     1                     * exp(-pligst(SRFC)*strlig(SRFC)) * dtm
c               where
c                   tcflow is grams C
c                   strucc(SRFC)  is the current value for total structural C
c                               (grams C) in layer SRFC
c                   strmax(SRFC) is the maximum amount of structural C in 
c                              layer SRFC that will decompose (grams C).
c                   defac is the decomposition factor based on water and 
c                         temperature computed in prelim and in cycle
c                   dec1(SRFC) is the intrinsic decomposition rate of 
c                            structural C (a fixed parameter)
c                   pligst is a fixed parameter that represents the effect 
c                          of lignin-to-structural-ratio on structural 
c                          decomposition
c                   strlig(SRFC) is the lignin content of surface structural 
c                             residue (grams lignin/grams biomass)
c                   dtm is the time step in years

c...Decompose structural into som1 and som2 with CO2 loss.
c     Changed csrsnk to st1c2 (10/92)
c         print *, 'tcflow = ', tcflow
         ligcon = strlig(SRFC)
         call declig(aminrl,ligcon,SRFC,nelem,2,ps1co2,rnewas,
     $        rsplig,tcflow,strucc,st1c2,strcis,struce,
     $        gromin,minerl,strmnr,resp,som1ci,som1e,som2ci,som2e)
         strlig(SRFC) = ligcon

      endif

c...Soil STRUCTURAL Material
      if (strucc(SOIL) .gt. 1.e-07) then

c...Compute total C flow out of structural in layer SOIL
c...Added impact of soil anerobic conditions -rm 12/91
           tcflow = min(strucc(SOIL),strmax(SOIL)) * defac * dec1(SOIL) 
     $              * exp(-pligst(SOIL)*strlig(SOIL)) * cltfac(4) * 
     $			anerb * dtm
c         print *, '2nd tcflow = ', tcflow


c               where
c                   tcflow is grams C
c                   strucc(SOIL)  is the current value for total structural C
c                               (grams C) in layer SOIL
c                   strmax(SOIL) is the maximum amount of structural C in 
c                              layer SOIL that will decompose (grams C).
c                   defac is the decomposition factor based on water and 
c                         temperature computed in prelim and in cycle
c                   dec1(SOIL) is the intrinsic decomposition rate of 
c                            structural C (a fixed parameter)
c                   pligst is a fixed parameter that represents the effect 
c                          of lignin-to-structural-ratio on structural 
c                          decomposition
c                   strlig(SOIL) is the lignin content of soil structural 
c                             residue (grams lignin/grams biomass)
c                   cltfac(4) is the cultivation factor for soil
c                             structural material (set in cycle)
c                   dtm is the time step in years (= dt/ntspm)

c...Decompose structural into som1 and som2 with CO2 loss.
c...Changed csrsnk to st2c2 (10/92)

          call declig(aminrl,strlig(SOIL),SOIL,nelem,2,ps1co2,rnewbs,
     $        rsplig,tcflow,strucc,st2c2,strcis,struce,
     $        gromin,minerl,strmnr,resp,som1ci,som1e,som2ci,som2e)


      endif
c...End of Structural Decomposition


c...METABOLIC Material

c...Process each layer
      do 30 lyr = SRFC, SOIL
c          print *, 'litdec: lyr = ', lyr
          if (metabc(lyr) .gt. 1.e-07) then

c...Determine C/E ratios for flows to SOM1
              do 10 iel=1,nelem

c...Compute ratios for surface metabolic residue
                  if (lyr .eq. SRFC) then
                      rceto1(iel) = agdrat(metabe(lyr,iel),metabc(lyr),
     $                              biocnv,iel,pcemic,cemicb)

c...Compute ratios for soil metabolic residue
                  else
                      rceto1(iel) = bgdrat(aminrl,varat1,iel)
                  endif
10            continue

c...Compute total C flow out of metabolic in layer lyr
              tcflow = metabc(lyr) * defac * dec2(lyr) * dtm

c...Added impact of soil anerobic conditions -rm 12/91
              if (lyr .eq. SOIL) then
                tcflow = tcflow * anerb
              endif

c       where:
c           tcflow is grams C
c           metabc(lyr) is the current value for total metabolic C
c             (grams C) in layer lyr
c           defac is the decomposition factor
c           dec2(lyr) is the intrinsic decomposition rate of 
c             metabolic C
c           dtm is the time step in years


c...Make sure metab does not go negative.
              if (tcflow .gt. metabc(lyr)) then
                tcflow = metabc(lyr)
              endif

c...If decomposition can occur,
              if (candec(nelem,aminrl,metabc(lyr),metabe,2,lyr,rceto1))
     $          then

c...CO2 loss
                  co2los = tcflow * pmco2(lyr)

c...Changed csrsnk to mt1c2, mt2c2 (10/92)
                  if (lyr .eq. SRFC) then
                    call respir(co2los,2,lyr,metabc,metcis,mt1c2,
     $                      resp,metabe,minerl,gromin,metmnr)
                  else
                    call respir(co2los,2,lyr,metabc,metcis,mt2c2,
     $                      resp,metabe,minerl,gromin,metmnr)
                  endif

c...Decompose metabolic into som1
                  cfmes1 = tcflow - co2los

c...Partition and schedule C flows by isotope
                  call csched (cfmes1,metcis(lyr,LABELD),metabc(lyr),
     $                         metcis(lyr,UNLABL),som1ci(lyr,UNLABL),
     $                         metcis(lyr,LABELD),som1ci(lyr,LABELD),
     $                         1.0,accum)

c...Compute and schedule N, P, and S flows and update mineralization accumulators.
                  do 20 iel = 1, nelem
             
                      call esched(cfmes1,metabc(lyr),rceto1(iel),
     $                          metabe(lyr,iel),som1e(lyr,iel),
     $                          minerl(SRFC,iel),mnrflo)
                      call mnracc(mnrflo,gromin(iel),metmnr(lyr,iel))
20                continue
              endif
          endif

c...End of Metabolic Decomposition


c...Next layer
30    continue

      return
      end
