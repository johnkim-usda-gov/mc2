


      subroutine potfor(month)
      integer month

c...Compute monthly potential production for forest

      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'potent.inc'
c...schvar.inc is included here because it has the code_flags
      include 'schvar.inc'
      include 'wth.inc'
      include 'vemap2.inc'
c      include 'mapssvar.inc'

c...RESPPT is an array of respiration values for the forest system
c     production components.

c...Added savanna model, pcropc now pforc
c                        tgprod now tfprod (BO)

c...Function declarations
      real     gpdf, frespr, laprod, pevap, pprdwc, pprdwc_old
      external gpdf, frespr, laprod, pevap, pprdwc, pprdwc_old

c...Local variables
      integer  i,mode
      real     fcmax, frlive, potprd, pptprd, resppt(5), 
     $         tfprod, water_yield,prodclim,ronsratio,sumdd
      save sumdd
      real pprdwcA, pprdwcB, pprdwcC 

            
	water_yield = 0.0
	ronsratio = 0.0
	prodclim = 0.0
	if(month .eq. 1) sumdd = 0.0

c...	initialization just in case Bachelet 10-8-97

c...Estimate potential production based on temp & h2o

c---C. Daly 12-19-97.  Moved respiration outside temperature loop.
c...Calculate monthly maintenance respiration
        frlive = sapk / (sapk + rlwodc)

        resppt(LEAF) = frespr(tave,rleave(N))
        resppt(FROOT) = frespr(stemp,froote(N))
        resppt(FBRCH) = frespr(tave,fbrche(N))
        resppt(LWOOD) = frespr(tave,frlive*rlwode(N))
        resppt(CROOT) = frespr(stemp,frlive*croote(N))

        sumrsp = 0.0
        do 10 i = 1, FPARTS
          sumrsp = sumrsp + resppt(i)
10      continue

      if (stemp .gt. 0.0) then
c...Calculate temperature effect on growth.
        potprd = gpdf(stemp, ppdf(1,2), ppdf(2,2), ppdf(3,2), 
     $                ppdf(4,2))
 
        fprd_tmp = potprd

c...Calculate moisture effect on growth -
c     Value for potential plant production is now calculated from the
c     equation of a line whose intercept changes depending on water
c     content based on soil type.  The function PPRDWC contains the
c     equation for the line.  -rm 9/91
c...    distributing monthly rainfall as a function of canopy size
c...    DB 11-20-98
 
        pet = pevap(month)
        if (pet .ge. .01) then
c	allow trees access to water in the profile -mdh 11/96

c...	testing sensitivity of the model to water availability 
c...	between grasses and trees  DB/2-11-98

c...      this will be run stds (standard code with shallow trees)

        pptprd = (avh2o(2) + prcurr(month) + irract)/pet
c      if (diags) then
c        write(*,*) 'month = ', month
c        write(*,*) 'avh2o(2), prcurr(month), irract, pet =',
c     $    avh2o(2), prcurr(month), irract, pet
c      endif
        else
         pptprd = .01
        endif
      if (pptprd .gt. 1.0) pptprd = 1.0
c...LENIHAN 4/28/05

        mi_tree = pptprd

c...LENIHAN 5/30/97

        mf_lmc(1) = pptprd
        if ( mf_lmc(1) .gt. 1.0 ) mf_lmc(1) = 1.0
        if ( pet .le. .01 ) mf_lmc(1) = 1.0

        call get_pprdwc3(pprdwcA, pprdwcB, pprdwcC, 1) 
        pptprd = pprdwc(wc, pptprd, pprdwcA, pprdwcB, pprdwcC)

c...LENIHAN 4/28/05 
 
        fprd_ppt = pptprd
 
c...For large wood, calculate the percentage which is live (sapwood)
c FERSCHWEILER 11/18/10 - moved this to before frlive gets used
c       frlive = sapk / (sapk + rlwodc)
 
c----C. Daly 12-19-97.  Moved new laicalc outside of temp loop.

c...calculating water yield
c   Bachelet 10-8-97

        if(fcacc .eq. 0.0) then
                water_yield = 1.0
c	write(*,*) 'fcacc=0 in potfor'
        else
                water_yield = tran/fcacc
        endif
c
c...calculating other diagnostic ratios:
c
	if(lait .ne. 0.0) then
		ronsratio = (rlvacc+fbracc+rlwacc)/ lait
	else
		ronsratio = 0.0
c	write(*,*) 'lait=0 in potfor'
	endif
c
c...check production vs degree days

c	
	if (maxtmp(month) .gt. 0.0 .and. mintmp(month) .gt. 0.0) then
		sumdd = sumdd + (maxtmp(month)+mintmp(month))/2.0	
		prodclim = fcacc *1.0/sumdd
	else
		sumdd = 0.0
		prodclim = 0.0
	endif
	
        tfprod = prdx(2)*potprd*pptprd
        tfprod = tfprod*laprod(lait,laitop)*co2cpr(FORSYS)

c----C. Daly 12-19-97.  Moved respiration calcs outside of temp loop.

c...Use 2.0 to convert from biomass to carbon in forest system
c     - Mike Ryan & Dennis Ojime
c...Added calculation of fcmax  -rm  11/91

        pforc = (tfprod / 2.0) - sumrsp
        fcmax = prdx(3) * potprd * pptprd * laprod(lait,laitop)

c...David King's fix
c...Added direct CO2 effect in NPP mode 7/2/12
        fcmax = fcmax * co2cpr(FORSYS)

        if (pforc .lt. 0.0) then
          pforc = 0.0
        endif
        if (pforc .gt. fcmax) then
          pforc = fcmax
        endif
        
c end of the if (stemp .gt. 0) block        
      else
        tfprod = 0.
        pforc = 0.

c-------C. Daly 12-19-97.  More initialized vars when stemp < 0.0
        potprd = 0.0
        pptprd = 0.0
        tfprod = 0.0
        
      endif

      return
      end
