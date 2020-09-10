	real function pprdwc(wc,x,A,B,C)

       real     wc, x, A, B, C

c...This function returns a scalar for potential plant production
c     due to water content.  Basically you have an equation of a
c     line with a moveable x-intercept depending on the soil type.
c     x is the ratio of available water to potential evapotranspiration
c     wc is the soil water holding capacity in the top soil layer,
c     normalized by the depth of the layer

c     A:  The maximum ratio of available water to pet which
c                 would completely limit production assuming wc=0.
c     B:  The effect of wc on the intercept, allows the
c                 user to increase the value of the intercept and
c                 thereby increase the slope of the line.
c     C:  The lowest ratio of available water to pet at which
c                 there is no restriction on production.

c...Local variabes
       real     intcpt
              

c...The equation for the x-intercept (intcpt) is A+B*wc.  A and B
c     determine the effect of soil texture on plant production based
c     on moisture.
c...So the first point on the line is (A+B*wc, 0.0).
c...The second point on the line is (C,1.0)

       intcpt = A + B*wc
       
       if (x .le. intcpt) then
         pprdwc = 0.
       elseif (x .ge. C) then 
         pprdwc = 1.
       else
         pprdwc = (x - intcpt)/(C - intcpt)
       endif

       return
       end
       


	real function pprdwc_old(wc,x,pprpts,mode)

       integer mode
       real     wc, x, pprpts(3),pppG,pppT
       real pprdwc

c...This funtion returns a value for potential plant production
c     due to water content.  Basically you have an equation of a
c     line with a moveable y-intercept depending on the soil type.
c     The value passed in for x is ((avh2o(1) + prcurr(month) + irract)/pet)

c     pprpts(1):  The minimum ratio of available water to pet which
c                 would completely limit production assuming wc=0.
c     pprpts(2):  The effect of wc on the intercept, allows the
c                 user to increase the value of the intercept and
c                 thereby increase the slope of the line.
c     pprpts(3):  The lowest ratio of available water to pet at which
c                 there is no restriction on production.

c...Local variabes
       real     intcpt, slope, tmp

c...The equation for the y-intercept (intcpt) is A+B*WC.  A and B
c     determine the effect of soil texture on plant production based
c     on moisture.

c...Old way:
c      intcpt = 0.0 + 1.0 * wc
c      The second point in the equation is (.8,1.0)
c      slope = (1.0-0.0)/(.8-intcpt)
c      pprdwc = 1.0+slope*(x-.8)

c...New way:
c... pprpts(3) is changed to .70 with tmp = .5  DB  10-8-98
c... pprpts(3) is changed to .82 with tmp = .75 DB 11-11-98

       intcpt = pprpts(1) + (pprpts(2) * wc)

c... DB 10/8/98 separate trees afrom grasses using mode (0 is grasses)
c... modifying the intercept of that funky curve

	tmp = .5
	pppT = .9
	pppG = .6
c       tmp = 1.0
c	pppT = 0.8
c	pppG = 0.8

	if(mode .eq. 0) then 
		slope = 1.0 / (pppG - intcpt)
       		pprdwc = 1.0 + slope * (x - pppG)
	elseif (mode .eq. 1) then
		slope = 1.0 / (pppT - intcpt)
		pprdwc = tmp + slope * (x - pppT)
	endif

       if (pprdwc .gt. 1.0) then
         pprdwc = 1.0
       elseif (pprdwc .lt. 0.01) then
         pprdwc = 0.01
       endif
       
       pprdwc_old = pprdwc

       return
       end
