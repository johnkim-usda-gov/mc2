c--passvars.F


      subroutine statein(state)

c-----Routine disburses Century state vars into proper global vars

      include 'const.inc'
      include 'parfx.inc'
      include 'parfs.inc'
      include 'forrem.inc'
      include 'parcp.inc'
      include 'param.inc'
c...schvar.inc is included here because it has the code_flags
      include 'schvar.inc'

      real state(200) 
      integer i

c-----LENIHAN 4/28/97
 
      mf_remf(1)   = state(1)
      mf_remf(2)   = state(2)
      mf_remf(3)   = state(3)
      mf_remf(4)   = state(4)
      mf_remf(5)   = state(5)
      mf_flfrem    = state(6)
      mf_fdfrem(1) = state(7)
      mf_fdfrem(2) = state(8)
      mf_dr(1)     = state(9)
      mf_dr(2)     = state(10)
      mf_dr(3)     = state(11)
      mf_fd(1)     = state(12)
      mf_fd(2)     = state(13)
c      did_burn     = state(14)
      if (state(14).eq.0.0) then
       	did_burn = .false.
      else
       	did_burn = .true.
      endif

      do 10 i=15,187
        state(i) = 0.0
10    continue

      runoff_x_intercept = state(188)
      runoff_slope = state(189)
      lait_low_lim = state(190)
      
      do 20 i = 191, 200
c       code_flags(i - 190) = state(i)
c-----FERSCHWEILER 10/4/10
      if (state(i) .eq. 0.) then
        code_flags(i - 190) = .false.
      else
        code_flags(i - 190) = .true.
      endif
20    continue

      return
      end

c--------------------------------------------------------------------
      subroutine stateout(state)

c-----Routine collects Century state vars for pass to MAPSS

      include 'const.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'forrem.inc'
      include 'parcp.inc'
      include 'param.inc'
      include 'vemap2.inc'
      include 'plot3.inc'

      real state(200) 
      integer i

c-----LENIHAN 4/28/97
 
cx      state(1) = mf_remf(1)
cx      state(2) = mf_remf(2)
cx      state(3) = mf_remf(3)
cx      state(4) = mf_remf(4)
cx      state(5) = mf_remf(5)
cx      state(6) = mf_flfrem
cx      state(7) = mf_fdfrem(1)
cx      state(8) = mf_fdfrem(2)
cx      state(9) = mf_dr(1)
cx      state(10) = mf_dr(2)
cx      state(11) = mf_dr(3)
cx      state(12) = mf_fd(1)
cx      state(13) = mf_fd(2)
cx      state(14) = did_burn

      do 2 i=1,14
        state(i) = 0.0
2     continue

      state(15) = wooddr(3)
      state(16) = wooddr(4)
      state(17) = decw1
c     changed it to decw2 DB 10-10-98
      state(18) = decw2
      state(19) = pligst(SRFC)

      state(20) = mf_lmc(1)
      state(21) = mf_lmc(2)

      state(22) = nlaypg
      state(23) = klai

      state(24) = nppx
      state(25) = rspx
      state(26) = nepx
      state(27) = aetx
c state(28) is a spare
      state(28) = 0
      state(29) = vegc
      state(30) = tslc
      state(31) = minx
      state(32) = litf
      state(33) = clsl 
      state(34) = needle_index
      state(35) = everg_index
      state(36) = mapss_lai
      state(37) = fcaccx
      state(38) = agcaccx
      state(39) = bgcaccx
      state(40) = frstcx
      state(41) = aglivcx
      state(42) = bglivcx
      state(43) = soilnmx
      state(44) = afcaccx
      state(45) = awhc
      state(46) = frootcx
      state(47) = crootcx
      state(48) = adeadc
      state(49) = bdeadc
      state(50) = tot_whc
      state(51) = fprd_ppt
      state(52) = fprd_tmp
      state(53) = gprd_ppt
      state(54) = gprd_tmp
      state(55) = mi_tree
      state(56) = mi_grass
      state(57) = bio_consume
      state(58) = rleavc
      state(59) = ddecid
      state(60) = nlayer
      state(61) = lfroot

      do 10 i=62,200
        state(i) = 0.0
10    continue

 
      return
      end

c--------------------------------------------------------------------

      subroutine stateinit(state)

c-----Initializes all state vars to zero.  Call statein()
c-----afterwards to disburse zeros to global state vars.

      real state(200) 
      integer i

      do 25 i=1,200
        state(i) = 0.0
25    continue

      return
      end

c--------------------------------------------------------------------

      subroutine varsout(outvars)

c-----Routine collects Century output vars for pass to MAPSS
      
c      include 'plot1.inc' /* 280 */
c      include 'plot2.inc' /* 107 */
c      include 'plot3.inc' /* 111 */
c      include 'plot4.inc' /* 15 */	

      include 'outval.inc'
      include 'param.inc'	/* outvars(550) = annet */

      real outvars(600)
      integer i, ii

      ii = 0

      do 10 i=1,280
        outvars(ii+i) = vals1(i)
10    continue

      ii = 280
      do 20 i=1,107
        outvars(ii+i) = vals2(i)
20    continue

      ii = 387
      do 30 i=1,111
        outvars(ii+i) = vals3(i)
30    continue

      ii = 498
      do 40 i=1,15
        outvars(ii+i) = vals4(i)
40    continue

      outvars(550) = annet
      
      return
      end

c--------------------------------------------------------------------
      subroutine varsin(outvars)

c-----Routine distributes Century output vars after pass from MAPSS
      
c      include 'plot1.inc' /* 280 */
c      include 'plot2.inc' /* 107 */
c      include 'plot3.inc' /* 111 */
c      include 'plot4.inc' /* 7 */

      include 'outval.inc'
      include 'param.inc'	/* outvars(550) = annet */

      real outvars(600)
      integer i, ii

      ii = 0

      do 10 i=1,280
        ii = ii + 1
        vals1(i) = outvars(ii) 
10    continue

      do 20 i=1,107
        ii = ii + 1
        vals2(i) = outvars(ii)
20    continue

      do 30 i=1,111
        ii = ii + 1
        vals3(i) = outvars(ii)
30    continue

      do 40 i=1,15
        ii = ii + 1
        vals4(i) = outvars(ii)
40    continue

      annet = outvars(550)

      return
      end

c--------------------------------------------------------------------
      subroutine varsinit(outvars)

c-----Routine initializes Century output vars.  Call varsin()
c-----to distribute zeros to global output vars.
      
      include 'outval.inc'

      real outvars(600)
      integer i

      do 10 i=1,600
        outvars(i) = 0.0
10    continue

      return
      end
c--------------------------------------------------------------------

      subroutine divcar

c-----Divides total carbon back into labeled and unlabeled
c-----components (labeled = 0).  Reverse of sumcar.

      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

      strcis(1,1) = strucc(1)
      strcis(1,2) = 0.0
      strcis(2,1) = strucc(2)
      strcis(2,2) = 0.0
      metcis(1,1) = metabc(1)
      metcis(1,2) = 0.0
      metcis(2,1) = metabc(2)
      metcis(2,2) = 0.0

      som1ci(1,1) = som1c(1)
      som1ci(1,2) = 0.0
      som1ci(2,1) = som1c(2)
      som1ci(2,2) = 0.0
      som2ci(1) = som2c(1)
      som2ci(2) = 0.0
      som3ci(1) = som3c(1)
      som3ci(2) = 0.0

      wd1cis(1) = wood1c
      wd1cis(2) = 0.0
      wd2cis(1) = wood2c
      wd2cis(2) = 0.0
      wd3cis(1) = wood3c
      wd3cis(2) = 0.0

      aglcis(1) = aglivc
      aglcis(2) = 0.0
      stdcis(1) = stdedc
      stdcis(2) = 0.0
      bglcis(1) = bglivc
      bglcis(2)= 0.0
      agcisa(1) = agcacc
      agcisa(2) = 0.0
      bgcisa(1) = bgcacc
      bgcisa(2) = 0.0

      rlvcis(1) = rleavc
      rlvcis(2) = 0.0
      frtcis(1) = frootc
      frtcis(2) = 0.0
      fbrcis(1) = fbrchc
      fbrcis(2) = 0.0
      rlwcis(1) = rlwodc
      rlwcis(2) = 0.0
      crtcis(1) = crootc
      crtcis(2) = 0.0

      alvcis(1) = rlvacc
      alvcis(2) = 0.0
      afrcis(1) = frtacc
      afrcis(2) = 0.0
      afbcis(1) = fbracc
      afbcis(2) = 0.0
      alwcis(1) = rlwacc
      alwcis(2) = 0.0
      acrcis(1) = crtacc
      acrcis(2) = 0.0

      return                                                                    
      end                                           



