
c*****************************************************************
c	$RCSfile: prelim.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: prelim.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.4  2005/01/20 22:52:44  lenihan
cchanged tot_whc assignment to water in rooted layers
c
cRevision 1.3  2004/12/03 18:38:20  lenihan
cadded include vemap2.inc
cadded tot_whc assignment
c
cRevision 1.2  2004/09/24 23:27:21  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine prelim

c...Initialize variables and parameters
 
      include 'comput.inc'
      include 'const.inc'
      include 'mapssvar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
c /plot4/ is included for asmos_subsoil
      include 'plot4.inc'
      include 'seq.inc'
      include 'site.inc' 
      include 't0par.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'zztim.inc'
      include 'vemap2.inc'

c...Function declarations
      real      atanf, tcalc
      external  atanf, tcalc

c...Local variables
      integer   i, iel, iso, lyr, m
      real      dely, delx, fcbd(6), fccl(6), fcin(6), fcom(6), 
     $          fcsa(6), fcsi(6), fcwp(6), ompc, xslope, 
     $          textur, tfunc, wfunc, wpbd(6), wpcl(6), wpin(6),
     $          wpom(6), wpsa(6), wpsi(6), wpwp(6), yint

c...swflag lets the model user choose between using actual data 
c     for awilt and afiel or equations from Gupta and Larson (1979) 
c     or Rawls et al (1982).

c     swflag=0 Use actual data
c     swflag=1 Use G&L for both awilt (-15 bar) and afiel (-0.33 bar)
c     swflag=2 Use G&L for both awilt (-15 bar) and afiel (-0.10 bar)
c     swflag=3 Use Rawls for both awilt (-15 bar) and afiel (-0.33 bar)
c     swflag=4 Use Rawls for both awilt (-15 bar) and afiel (-0.10 bar)
c     swflag=5 Use Rawls for afiel (-0.33 bar) and actual data for awilt
c     swflag=6 Use Rawls for afiel (-0.10 bar) and actual data for awilt

c...       swflag   1          2          3        4       5       6
      data fcsa / 0.3075,    0.5018,   -0.20,   -0.30,  -0.19,   0.31/
      data fcsi / 0.5886,    0.8548,    0.0,     0.0,    0.0,    0.0/
      data fccl / 0.8039,    0.8833,    0.36,    0.23,   0.0,    0.0/
      data fcom / 2.208E-03, 4.966E-03, 0.0299,  0.0317, 0.0210, 0.0260/
      data fcbd /-0.1434,   -0.2423,    0.0,     0.0,    0.0,    0.0/
      data fcwp / 0.0,       0.0,       0.0,     0.0,    0.72,   0.41/
      data fcin / 0.0,       0.0,       0.2576,  0.4118, 0.2391, 0.4103/
      data wpsa /-0.0059,   -0.0059,    0.0,     0.0,    0.0,    0.0/
      data wpsi / 0.1142,    0.1142,    0.0,     0.0,    0.0,    0.0/
      data wpcl / 0.5766,    0.5766,    0.50,    0.50,   0.0,    0.0/
      data wpom / 2.228E-03, 2.228E-03, 0.0158,  0.0158, 0.0,    0.0/
      data wpbd / 0.02671,   0.02671,   0.0,     0.0,    0.0,    0.0/
      data wpwp / 0.0,       0.0,       0.0,     0.0,    1.0,    1.0/    
      data wpin / 0.0,       0.0,       0.0260,  0.0260, 0.0,    0.0/
      
c...Time initializations -  time step is one month
      dt = 1.0/12.0
      time = strtyr
      month = 0
      dayhrs = 0.0

c...Allow for time step < 1 month for running decomp
c     ntspm is the number of time steps per month for decomp
c     (read from the fix.100 file)
c     decodt is the time step used in subroutine decomp
      decodt = dt/real(ntspm)
 
c...Initializations
      crpgrw = 0
      seedl = 0
      forgrw = 0
      falprc = 0
      elimit = 0.

c...Initialize volitalization accumulators
      volgma = 0.0
      volexa = 0.0
      volpla = 0.0
 
c...Initialize erosion variables
      scloss = 0.0
      sclosa = 0.0
 
c...Initialize total co2 loss
      totco2 = 0.0
	
c...Initialize accumulators
      call annacc
      call inprac

c...Open the c14 data file 
      if (labtyp .eq. 1) then
c        open(unit=10,file='c14data',status='OLD')
         stop 'Unexpected code for MC1'
      endif

c...Calculate C,N,P,S in lower horizon soil pools for use as soil
c     replacement with erosion events
      do 3 iso = 1, 2
        lhzci(1,iso) = som1ci(SOIL,iso)*lhzf(1)
        lhzci(2,iso) = som2ci(iso)*lhzf(2)
        lhzci(3,iso) = som3ci(iso)*lhzf(3)
3     continue
      do 5 iel = 1, nelem
        lhze(1,iel) = som1e(SOIL,iel)*lhzf(1)
        lhze(2,iel) = som2e(iel)*lhzf(2)
        lhze(3,iel) = som3e(iel)*lhzf(3)
5     continue

c...Field capacity and wilting point.  Computations based on
c     Gupta and Larson 1979, 'Estimating soil and water retention
c     characteristics from particle size distribution, organic 
c     matter percent and bulk density'. Water Resources Research 15:1633
c     or Rawls et al (1982) 'Estimation of soi water properties'
c     Trans. ASAE ???:1316
c...Field capacity options of -0.1 or -0.33 bar.
c...Wilting point assumed to be water content at -15 bars.
c...Calculate organic matter from initial conditions, ivauto or 
c     or value at the beginning of an extend
c...Note that Gupta and Larson and Rawls use % for texture
c     but values here are fractions.

      if (swflag .ne. 0) then
	if(bulkd .eq. 0.) then
		ompc = 0.0
	else
        	ompc = somsc*1.724/(10000*bulkd*edepth)
	endif
      if (diags) then
      write(*,*) 'sand, clay, silt, bulkd = ', 
     .    sand, clay, silt, bulkd
      write(*,*) 'somsc, edepth, ompc = ', somsc, edepth, ompc
      endif
        do 10 lyr = 1, nlayer
          afiel(lyr) =   fcsa(swflag)*sand  + fcsi(swflag)*silt
     $                 + fccl(swflag)*clay  + fcom(swflag)*ompc
     $                 + fcbd(swflag)*bulkd + fcwp(swflag)*awilt(lyr)
     $                 + fcin(swflag)
          awilt(lyr) =   wpsa(swflag)*sand  + wpsi(swflag)*silt
     $                 + wpcl(swflag)*clay  + wpom(swflag)*ompc
     $                 + wpbd(swflag)*bulkd + wpwp(swflag)*awilt(lyr)
     $                 + wpin(swflag)
          ompc = ompc * 0.85
          xtemp = afiel(lyr)*1.e5 + 0.5
          itemp = anint(xtemp)
          afiel(lyr) = float(itemp)/1.e5
          xtemp = awilt(lyr)*1.e5 + 0.5
          itemp = anint(xtemp)
          awilt(lyr) = float(itemp)/1.e5
      if (diags) then
      write(6,9755) 
     .    lyr, afiel(lyr), awilt(lyr)
 9755 format('lyr, afiel(lyr), awilt(lyr) = ', i2, 2f15.12)
      endif
   10   continue
      endif
        
c...Added calculation for water content which will be used to
c     determine plant production in POTGRS. 10-90 -rm
      wc = afiel(1)-awilt(1)

c...Re-calculate awhc for the first crop (also done in cropin when crops
c     are changed)
      awhc  = 0.0
      do 12 lyr = 1, nlaypg
        awhc = awhc + (afiel(lyr) - awilt(lyr)) * adep(lyr)
12    continue

c      if (???) then
c        print *, '*** prelim.F: awhc = ', awhc
c        print *, 'afiel[] = ', afiel
c        print *, 'awilt[] = ', awilt
c        print *, 'adep[] = ', adep
c        print *, 'sand, silt, clay, bulkd = ', 
c     $      sand, silt, clay, bulkd
c        print *, 'nlaypg = ', nlaypg
c        stop
c      endif

      tot_whc = awhc

c...Calculate total water holding capacity 11/91 lh
      twhc = 0.0
      do 15 lyr = 1, nlayer
        twhc = twhc + (afiel(lyr) * adep(lyr))
15    continue

c...Calculate initial asmos based on relative water content (rwcf)
 
      do 17 lyr = 1, 10
        asmos(lyr) = 0.0
17    continue
      asmos_subsoil = 0.0

      do 18 lyr = 1, nlayer
        asmos(lyr) = (awilt(lyr) + rwcf(lyr) * 
     +                 (afiel(lyr) - awilt(lyr))) * adep(lyr)
18    continue

c...Computations related to decomposition of soil organic matter
c          Added 08/91   vek
      call predec(sand)

c...Initialize the layer beyond the last one the used for safety

      do 16 lyr= 1,11

      do 19 iel = 1, 2
        minerl(lyr, iel) = 0.0
19    continue

16    continue

c...Intercept for the texture equation of secondary P depends upon
c     pH input.  Note that this code changes the value of a
c     'fixed' parameter (texesp(2))
      if (ph .le. phesp(1)) then
	texesp(2) = phesp(2)
      else if (ph .ge. phesp(3)) then
        texesp(2) = phesp(4)
      else
 	dely = phesp(4) - phesp(2)
	delx = phesp(3) - phesp(1)
	xslope = dely / delx
	yint = phesp(2) - (xslope*phesp(1))
	texesp(2) = (xslope*ph) + yint
      endif

      if (micosm .eq. 0) then
c...Preset array which will contain monthly values of defac
        do 20 i = 1, MONTHS
          defacm(i) = -1.
   20   continue
        adefac = 0.
        defac = 0.
      else
	if (cursys .eq. FORSYS .or. cursys .eq. SAVSYS) then
          call message(' ')
          call message('   Microcosm works for cropping system only.')
c...There is currently no check that only Fallow is chosen
          STOP
        endif
        stemp = mctemp
        tfunc = tcalc(stemp, teff)
        wfunc = 1./(1.+4.*exp(-6.*rwcf(1)))
		if (wfunc .ne. wfunc) then
		  print *, 'prelim: wfunc is a NaN. rwcf(1) = ', rwcf(1)
		  stop
		endif
    
        defac = tfunc*wfunc
c      write(*,*) 'tfunc, wfunc, defac = ', tfunc, wfunc, defac
c...Bound defac to >= 0.0 12/21/92
         if (defac .lt. 0.0) defac = 0.0
        do 30 m = 1, MONTHS
          prcurr(m) = mctemp
          prcnxt(m) = mctemp
          defacm(m) = defac
   30   continue

        do 40 lyr = 1,nlayer
          amov(lyr) = 0.0
   40   continue
      endif

c...Effect of soil texture on the microbe decomposition rate
      eftext = peftxa+peftxb*sand

c...Compute parameters which control decomposition of som1
c     p1co2 must be computed for surface and soil.   vek  08/91
c...Note that p1co2b(1) must equal 0 because there is no
c     soil texture effect on the surface.
      p1co2(SRFC) = p1co2a(SRFC)
      p1co2(SOIL) = p1co2a(SOIL)+p1co2b(SOIL)*sand

c...Decomposition of som1 to som3 is a function of clay content
c             vek june90
      fps1s3 = ps1s3(1) + ps1s3(2) * clay
      fps2s3 = ps2s3(1) + ps2s3(2) * clay

      if (texepp(1) .eq. 1.0) then

c...Calculate pparmn(2)
c     Include effect of texture; weathering factor should be per year
        textur = clay + silt
        pparmn(2) = 12.0 * atanf(textur, texepp(2), texepp(3),
     +                           texepp(4), texepp(5))
      endif

      if (texesp(1) .eq. 1.0) then

c...Calculate psecmn(2)
c     Include effect of texture
        psecmn(2) = 12.0 * (texesp(2) + texesp(3) * sand)
      endif

c...Save initial values for printing or plotting
      call savarp

c...Clear the flow stack.
      call floclr

      return
      end
