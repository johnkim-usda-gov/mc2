      
c*****************************************************************
c	$RCSfile: woodec.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: woodec.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:54:43  bachelet
cheader
c
c******************************************************************
	subroutine woodec (dtm)

      real      dtm

c...Wood decomposition 			written by vek 04/91

c     defac = decomposition factor based on water and temperature
c             (computed in prelim and in cycle)
c     pligst = fixed parameter that represents the effect of
c              of lignin-to-structural-ratio on structural
c              decomposition

      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

c...Local variables
      real      tcflow
      real wood1c_array(1), wood2c_array(1), wood3c_array(1)


c...FINE BRANCHES
c     wood1c = C in dead fine branch component of forest system (g/m2)
c     decw1 = intrinsic rate of decomposition of dead fine branches
c     wdlig(FBRCH) = lignin fraction for fine branches

      if (wood1c .gt. 1.e-07) then

c...Compute total C flow out of fine branches
          tcflow = wood1c * defac * decw1 * exp(-pligst(SRFC) * 
     $             wdlig(FBRCH)) * dtm

c...Decompose fine branches into som1 and som2 with CO2 loss.

          wood1c_array(1) = wood1c
          call declig(aminrl,wdlig(FBRCH),SRFC,nelem,1,ps1co2,rneww1,
     $        rsplig,tcflow,wood1c_array,csrsnk,wd1cis,wood1e,
     $        gromin,minerl,w1mnr,resp,som1ci,som1e,som2ci,som2e)
          wood1c = wood1c_array(1)
      endif


c...LARGE WOOD
c     wood2c = C in dead large wood component of forest system (g/m2)
c     decw2 = intrinsic rate of decomposition of dead large wood
c     wdlig(LWOOD) = lignin fraction for large wood

      if (wood2c .gt. 1.e-07) then

c...Compute total C flow out of large wood
          tcflow = wood2c * defac * decw2 * exp(-pligst(SRFC) * 
     $             wdlig(LWOOD)) * dtm

c...Decompose large wood into som1 and som2 with CO2 loss.

          wood2c_array(1) = wood2c
          call declig(aminrl,wdlig(LWOOD),SRFC,nelem,1,ps1co2,rneww2,
     $        rsplig,tcflow,wood2c_array,csrsnk,wd2cis,wood2e,
     $        gromin,minerl,w2mnr,resp,som1ci,som1e,som2ci,som2e)
         wood2c = wood2c_array(1)
      endif


c...COARSE ROOTS
c     wood3c = C in dead coarse root component of forest system (g/m2)
c     decw3 = intrinsic rate of decomposition of dead coarse roots
c     wdlig(CROOT) = lignin fraction for coarse roots

      if (wood3c .gt. 1.e-07) then

c...Compute total C flow out of coarse roots.
          tcflow = wood3c * defac * decw3 * exp(-pligst(SOIL) *
     $             wdlig(CROOT)) *  anerb * dtm


c...Decompose coarse roots into som1 and som2 with CO2 loss.
          wood3c_array(1) = wood3c
          call declig(aminrl,wdlig(CROOT),SOIL,nelem,1,ps1co2,rneww3,
     $        rsplig,tcflow,wood3c_array,csrsnk,wd3cis,wood3e,
     $        gromin,minerl,w3mnr,resp,som1ci,som1e,som2ci,som2e)
          wood3c = wood3c_array(1)

      endif

      return
      end
