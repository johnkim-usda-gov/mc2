	subroutine droot(pltlig)

      real      pltlig(2)

c...Simulate death of roots for the month.
c
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfx.inc'
      include 'plot1.inc'

c...Local variables
      integer   iel
      real      fr14, recres(MAXIEL), rdeath, rtdh
	  real mortality


c...Death of roots

      if (bglivc .le. 0.) return 

      if (stemp .le. -5.) then 
c       root death in cold soil
        mortality = -.25 - .06 * stemp
      else
c       death of fine surface roots due to the top soil layers drying out
        rtdh = 1.0 - avh2o(1)/(deck5+avh2o(1))
        mortality = rdr * rtdh
      endif
      mortality = MIN(mortality, 0.95)	  
      rdeath = mortality * bglivc

      do 10 iel = 1, nelem
        recres(iel) = bglive(iel)/bglivc
10    continue
      fr14 = bglcis(LABELD)/bglivc
      call partit(rdeath,recres,2,bglcis,bglive,pltlig(BELOW),fr14)

      return
      end
