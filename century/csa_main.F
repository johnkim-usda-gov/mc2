c
c	$RCSfile: csa_main.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: csa_main.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.6  2005/04/16 15:54:51  conklind
cPut the logic to read the co2ramp file into a subroutine read_co2ramp() so that it can be called by cen_restore().
c
cRevision 1.5  2004/08/30 16:47:03  conklind
cImplement /save_vars/.
c
cRevision 1.4  2004/05/19 18:16:46  conklind
cFix frootcx typo.
c
cRevision 1.3  2004/05/19 17:05:59  conklind
cDB's reconciliation.
c
cRevision 1.2  2003/11/07 22:37:36  conklind
cAdd 3 arguments to cen_step calling sequence: warmstart_init_flag, time_in, and rptyrs_in.
cAt beginning of cen_step, when warmstart_init_flag is true, then copy time_in to time,
crptyrs_in to rptyrs, and vveg_in to vveg.
c
cRevision 1.1  2003/11/05 19:52:09  conklind
cInitial revision
c
c
c23456789012345678901234567890123456789012345678901234567890123456789012

c		Copyright 1993 Colorado State University
c			All Rights Reserved





c			DISCLAIMER
c
c	   Neither the Great Plains System Research Unit - USDA (GPSR) nor
c	Colorado State University (CSU) nor any of their employees make
c	any warranty or assumes any legal liability or responsibility for
c	the accuracy, completeness, or usefulness of any information,
c	apparatus, product, or process disclosed, or represents that its
c	use would not infringe privately owned rights.  Reference to any
c	special commercial products, process, or service by tradename,
c	trademark, manufacturer, or otherwise, does not necessarily
c	constitute or imply endorsement, recommendation, or favoring by  
c	the GPSR or CSU.  The views and opinions of the authors do not
c	necessarily state or reflect those of GPSR or CSU and shall not 
c	be used for advertising or product endorsement. 

c...Century Soil Organic Matter Model
c     Simulation of carbon, nitrogen, phosphorous, and sulfur cycling
c     As of Dec. 1991, uses a 1 month time step
c     Project - Soil Fertility in the Great Plains
c     Modeler - Bill Parton
c     Programmers - Vicki Kirchner, Becky McKeown, Laura Harding

c...State variables and flows are grams/m2.
                                           
c...           (unit 1) = plot/print file used by modaid (unformatted)
c...'site'.100 (unit 7) = parameter values and initial values for 
c                           state variables; see subroutine sitein.
c...fix.100    (unit 8) = fixed parameter values values for 
c                           state variables; see subroutine fixin.
c...           (unit 9) = a file of weather data read in subroutines 
c                           wthini, weathr
c...c14data   (unit 10) = a data file specifying years in which labeled
c                         carbon is added via plant growth and what 
c                         fraction of the growth is labeled.

c...If you're getting floating point errors mentioned after you exit
c	Century, uncomment the following lines, recompile, run Century
c	in dbx with the 'catch FPE' option to find the offending code.
c	You can also run Century outside of dbx, in which case you will
c	get messages on your screen giving you an fpe error code (see
c	the Floating Point Programmer's Guide, p20) and a not-very-
c	useful-to-3rd-or-4th-generation-language-programmers location. 
c	The error handler 'mysigal' is a Fortran callable C routine 
c	written by MartinFowler; it can be replaced by any user written
c	handler or any of several library handlers, the most useful 
c	probably being SIGFPE_ABORT.  The calls to ieee_handler won't 
c	compile using poa's binaries.

c      external mysigal
c      ieeer=ieee_handler('set','invalid',SIGFPE_ABORT)
c      ieeer=ieee_handler('set','division',mysigal)
c      ieeer=ieee_handler('set','overflow',mysigal)
c      ieeer=ieee_handler('set','underflow',SIGFPE_ABORT)

c...You probably wont want to uncomment the following line; inexact
c     floating point signals occur all over the place.

c       ieeer=ieee_handler('set','underflow',mysigal)

c...Obtain startup information from user, do initializations based on
c     answers to Modaid questions

      /****
      This is approximately where the normal Century would start.  I
      am simply going to do an entry here to get things going and
      will then call down lower for the annual step.
      *****/
      /****
      This also serves as the contact point between mapss and century
      when running std century to "eq" before beginning one_step.  
      ****/

c-----C. Daly 2-5-96 - passing years to run
c-----C. Daly 2-21-96 - passing ppt, tmaxi, and tmin
c-----C. Daly 11-96 - passing latitude and soil bd, sand, clay

c-----Need to pass number of soil layers (depth)--------

c-----C. Daly 11-22-96
c-----Passing a mix index, describing mix between subtropical
c-----evergreen broadleaf, temperate deciduous, and boreal 
c-----evergreen needleleaf for eastern US.

      subroutine cen_init(yr, vveg, diag_flag, os_flag,
     .                    state, outvars, path, mx, c3c4, initflag, 
     .                    fm, zn, tmpi, ppti, unlimitedNflag, fi,
     .                    co2ramp_path) 
      include 't0par.inc'
      include 'timvar.inc'
      include 'zztim.inc'
      include 'pass.inc'
      include 'plot3.inc'
      include 'plot1.inc'
      include 'const.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'schvar.inc'
      include 'chrvar.inc'
      
      include 'parfx.inc'
      include 'mapssvar.inc'
      include 'plot4.inc'
	  include 'vemap2.inc'

c-----Declare passed vars. 
      integer   yr, vveg, fm, zn
      logical   diag_flag, os_flag, initflag, unlimitedNflag
      real      state(200)
      real      outvars(600)
      character*100 path, co2ramp_path
      real mx, c3c4, tmpi, ppti, fi

      integer getlen
      integer len

c-----Clear filenames
c-----Be careful of fname, which is built in C.  A null character
c-----is needed for termination on strcpy in vetofix.
      sname(1:15) = '               '
      ename(1:15) = '               '
      nname(1:15) = '               '
      rname(1:7)  = '       '
      fname(1:14) = '              '
      fname(14:15) = char(0)

c-----Set global run flags, years to run, file path, tree mix, 
c-----fire mode, and climate zone. 
      unl_N_flag = unlimitedNflag
      diags = diag_flag
      onestep = os_flag
      years_to_run = yr
      filpath = path
      len = getlen(co2ramp_path)
      co2ramp_name(1:len) = co2ramp_path
c-----Fire mode can be enabled only in onestep mode
      if (onestep) then
        fire_mode = fm
      else
        fire_mode = 0
      endif
c-----Climate zone C. Daly 2-12-98.
      zone = zn
	  frost_index = fi

c-----If initializing from an EQ run, get ending vveg type.  This
c-----will overwite the vveg passed to this routine from MAPSS,
c-----which is the MAPSS veg calculated for year 1 climate only.

c-----C. Daly 4-30-97.  Also, set mix index for deciduous forest
c-----from init file, rather than use first year's calculation.
c-----If not initializing, set global tree mix index to current
c-----value.

c-----Save original vveg type.
      ovveg = vveg

c-----Get file directory length.
      pathlen = getlen(filpath)

c-----Build up Century file names from Vemap veg type.
c-----For standard century, prefix is "CenStd".
c-----For one_step mapss_century, prefix is "MapssCen".
c-----Both standard and one-step need schedule and output binary names.
c-----Only one-step needs extension binary name.
c-----Fix.100 file needed for both versions.

      if(onestep) then
c-------New and old (extension) binary filenames
        nname = 'MapssCen'
        ename = 'MapssCen'
c-------Schedule file (biome-specific)
        sname(1:8) = 'MapssCen'
        if(vveg.ge.10) then
          write(sname(9:10),301) vveg
        else
          write(sname(9:9),302) vveg
        endif
      else
c-------New binary filename
        nname = 'CenStd'
c-------Schedule file (biome-specific)
        sname(1:8) = 'vvegType'
        if(vveg.ge.10) then
          write(sname(9:10),301) vveg
        else
          write(sname(9:9),302) vveg
        endif
      endif

c-----Fixed parameters filename
      call vetofix(vveg, fname)
      if (diags) write(*,*) 'fname = ',fname

c-----Initialize state vars (need this anymore?).
      runoff_x_intercept = state(188)
      runoff_slope = state(189)
      lait_low_lim = state(190)
      do 20 i = 191, 200
        if (state(i) .eq. 0.) then
          code_flags(i - 190) = .false.
        else
          code_flags(i - 190) = .true.
        endif
20    continue

      call stateinit(state)

c-----If no netCDF initialization, set outvars to zero.
      if (initflag .eqv. .false.) then
        call varsinit(outvars)
        call varsin(outvars)
      endif
          
c-----Transfer mix indexes to global var for use here, and for output.
cx    mx_index = mx
      c3c4_index = c3c4
      tmp_index = tmpi
      ppt_index = ppti
      if (diags) then
        write(*,*) 'Temp mix index is ',tmp_index,'.'  
        write(*,*) 'Ppt mix index is ',ppt_index,'.'  
        write(*,*) 'C3C4 mix index is ',c3c4_index,'.'
      endif
      
c-----Do century initialization.
10    call detiv()
        
c-----If netCDF init, pass in outvars from array, outvars from
c-----detiv will now be overwritten.
      if (initflag) call varsin(outvars)
        
c-----Make sure month is set to 1 for onestep runs.
      if (onestep) month = 1

c-----Report back to MAPSS with new outvars and state.
      call stateout(state)
      call varsout(outvars)

301   format(i2)
302   format(i1)
      end
c-------------------------------------------------------------------


      subroutine stand_step(state, outvars, eq_flag, final_growth,
     .                      diag_flag, os_flag, fi)

c-----Routine runs century in standard mode, once model has been
c-----initialized.  Called by MAPSS each month, returns state
c-----and output vars, and equilibrium flag. 

c-----C. Daly 5-6-96.  Split out from csa_main.

      include 'zztim.inc'
      include 'schvar.inc'
      include 'timvar.inc'
      include 't0par.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'plot4.inc'
      include 'param.inc'
      include 'wth.inc'
      include 'site.inc'
      include 'parfx.inc'

      include 'mapssvar.inc'
	  include 'vemap2.inc'

      real      state(200)
      real      outvars(600)
      logical   eq_flag
      real      final_growth
	  real fi
      integer   i
      logical   diag_flag, os_flag

c-----Set diags and onestep flags.
      diags = diag_flag
      onestep = os_flag
      frost_index = fi
      do 22 i = 191, 200
        if (state(i) .eq. 0.) then
          code_flags(i - 190) = .false.
        else
          code_flags(i - 190) = .true.
        endif
22    continue
      if (.not. code_flags(9)) then
        runoff_x_intercept = state(188)
        runoff_slope = state(189)
        lait_low_lim = state(190)
      endif
      
c        write(*,*) 'time, som2c(1) = ', time, som2c(1)
	  
      if(time .le. 1.0 + (0.5 * dt)) then
       if (diags) then 
        write(*,*) 'Initial values of selected variables:'
        write(*,*)
        write(*,*) 'precip = ', precip
        write(*,*) 'tmx2m = ', tmx2m
        write(*,*) 'tmn2m = ', tmn2m
        write(*,*) 'site latitude = ', sitlat
	write(*,*) 'site elevation= ', sitelev
        write(*,*) 'sand = ', sand
        write(*,*) 'silt = ', silt
        write(*,*) 'clay = ', clay
        write(*,*) 'bulk density = ', bulkd
        write(*,*) 'hist n deposition = ', hdp
        write(*,*) 'current n deposition = ', cdp
        write(*,*) 'number of soil layers = ', nlayer
        write(*,*) 'number rooted soil layers = ', nlaypg
        write(*,*)
        write(*,*) 'layer effective thicknesses: '
        write(*,*) '----------------------------'
        do 89 i = 1, nlayer
          write(*,*) i, ' thickness = ',adep(i)
89      continue
        write(*,*)
       endif
      endif

c.....Update month
      month = mod(month,12) + 1

c...If time is greater than the ending time for the current block,
c     read the next block
      if ((abs(time - blktnd) .lt. (0.5 * dt)) .and.
     $     (abs(time - tend)   .gt. (0.5 * dt))) then

         call readblk()
      endif

c...Perform annual tasks
      if (month .eq. 1) then
         call eachyr()
      endif

c...The main driver for the model; call decomposition, growth, etc.
      call simsom()

c...Update time
      time = time + dt
      if (time .ge. -1.0e-07 .and. time .le. 1.0e-07) then
         time = 0.0
      endif

c-----C. Daly - reactivate output interval dtpl, read from schedule file

c...Write out values
      if (tplt - time .lt. dt * 0.5) then
         if (gobin) call wrtbin(time)       
         tplt = time + dtpl
      endif

c-----Write output vars and state to MAPSS each time step.
      call stateout(state)
      call varsout(outvars) 

c-----C. Daly - 2/23/96.  For standard Century runs, run model only until
c-----slow pool organic matter (som2c) achieves relative equilibrium.

      call eq_test (eq_flag, 200, final_growth)
      if (eq_flag) then

c-------We're done. Equilibrium reached.
       if(diags) then
         write(*,*) 'time, tend, rptyrs = ', time, tend, rptyrs
         write(*,*)'Soil slow pool SOM has reached equilibrium,'
         write(*,*)' or time is within rptyrs of tend.'
         endif

      endif

      return 
      end
c---------------------------------------------------------------------

      /***********
      C. Daly 2-5-96.  We are now going to use a mapss-century schedule
      file, rather than the old EQ one, for transient runs.  Anything
      in this schedule file that is not overwritten by info from mapss
      will take effect as usual.  We are targeting grazing and fire
      as specific examples of events we might want to prescribe in
      the transient run.  
      ************/

      subroutine cen_step (
     .     vveg_in,
     $     wbegin, wend,		/* woody growing season */
     $     cbegin, cend, csene,		/* grass growing season */
     $     wfire, cfire,		/* fire names */
     $     treec, cropc,		/* results sent back to MAPSS */
     $     state, outvars,		/* I/O arrays */
     $	   mx,	mx_out,			/* I/O Decid/Evergreen mix index*/
     $	   c3c4, c3c4_out,		/* I/O C3/C4 mix index*/
     $     tname, cname,                /* Century tree and crop names */
     $     diag_flag, os_flag,		/* Flags for diags and onestep */
     $     burn_out,			/* Output fire flag */
     $     burn_count,			/*Years since last burn */
     $     tmpi, tmpi_out,		/* Tmp-based mix index */
     $     ppti, ppti_out,		/* Ppt-based mix index */
     $     zn,                          /* climate zone defined in mapss TempZone */
     $     warmstart_init_flag,		/* TRUE=>1st call of a warmstart */
     $     time_in,			/* warmstart: initial time */
     $     rptyrs_in,			/* warmstart: rptyrs */
     $     unlimitedNflag,
     $     fi)                  /* frost_index */

      character*5 treename, cropname, tname, cname, wfire, cfire
cx      character*5 new_tn1, new_tn2
      integer wbegin, wend, cbegin, cend, csene, vveg_in, zn
	  integer :: vveg = 0	  
      integer rptyrs_in
      real mx, c3c4, mx_out, c3c4_out, burn_count,tmpi, tmpi_out,
     .     ppti, ppti_out, time_in, fi
      logical vegtype_burn, burn_out, warmstart_init_flag
      logical unlimitedNflag
cx      logical disturb
	integer getlen

      /*************
      These are the things that are actully sent back to MAPSS.  They
      are the leaf carbon for trees and grasses.  The actual assignments
      occur at the very bottom of this routine.  The returned numbers
      come from rleavc and aglivc.
      **************/

      real treec, cropc

      include 't0par.inc'/* Keep a close eye on these, I may not */
      include 'timvar.inc'/* need all of them in both the init and */
      include 'zztim.inc'	/* step routines. */
      include 'pass.inc'
      include 'schvar.inc'
      include 'chrvar.inc'
      include 'const.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'plot4.inc'
      include 'param.inc'
      include 'wth.inc'
      include 'parfx.inc'

      include 'mapssvar.inc'
      include 'vemap2.inc'

      integer	tevent
      save tevent
      integer 	ctemp, ttemp
      integer   prevsys
      equivalence (prevsys, prevsys_cen_step) /* in /save_vars/ */
      integer   crtyr, i
      real      state(200)
      real      outvars(600)
      integer cen_burn_count

      logical   diag_flag, os_flag
      
      integer cfireLen, wfireLen,j

      cen_burn_count = INT(burn_count)

c-----Set diags and onestep flags.
      diags = diag_flag
      onestep = os_flag
      unl_N_flag = unlimitedNflag
      zone = zn
	  frost_index = fi
      
      if (warmstart_init_flag) then
        time = time_in
	rptyrs = rptyrs_in
	vveg = vveg_in
      endif

c-----Pass in state from MAPSS each time step.
      call statein(state)
	  
      if(time .le. 1.0 + (0.5 * dt)) then
       if (diags) then
        write(*,*) 'Initial values of selected variables:'
        write(*,*)
        write(*,*) 'precip (cm) = ', precip
        write(*,*) 'tmx2m (C) = ', tmx2m
        write(*,*) 'tmn2m (C) = ', tmn2m
        write(*,*) 'site elevation= ', sitelev
        write(*,*) 'site latitude = ', sitlat
        write(*,*) 'sand = ', sand
        write(*,*) 'silt = ', silt
        write(*,*) 'clay = ', clay
        write(*,*) 'bulk density = ', bulkd
        write(*,*) 'hist n deposition = ', hdp
        write(*,*) 'current n deposition = ', cdp 
        write(*,*) 'number of soil layers = ', nlayer
        write(*,*) 'number rooted soil layers = ', nlaypg
        write(*,*)
        write(*,*) 'layer effective thicknesses: '
        write(*,*) '----------------------------'
        do 89 i = 1, nlayer
          write(*,*) i, ' thickness = ',adep(i)
89      continue
        write(*,*)
       endif
      endif

c     Each year I will completely recreate the events in the cmdary().
c     It will contain the veg stuff, the growing season stuff, and
c     any fire stuff.
c     Each month I may add events into the cmdary(), fire stuff.

c-----C. Daly 4-29-96.  First time step was coming up with crtyr of
c-----0.  This caused a mismatch between events here (crtyr=0) and
c-----event search in schedl (crtyr=1).  Now, they match.
c-----Changed again 6-5-96.  moved (+ 1) outside modulus calc.
c-----Prevents getting crtyr = 0 during last year of simulation.
c-----Not sure if it matters that its different than in schedl.

      crtyr = aint(time + .001)
      crtyr = mod((crtyr - strtyr), rptyrs) + 1

cx      if (diags) write(*,*) __LINE__, ' cen_step called ', month, time


c*************Beginning of the year.  Do schedule stuff.**************

      if (month .eq. 1) then

c...        write(*,*) 'time = ',time,' current year = ',crtyr

c-------C. Daly 6-16-97. Set disturbance flag for hardwiring tree types.
cx	if (crtyr .le. 1 ) disturb = .false.
cx        if (vegtype_burn .eq. .true.) disturb = .true.

c-------C. Daly 10-29-96.  Move this year's scheduled events
c-------from storage arrays into the working arrays.

        call update_sched(crtyr)
        tevent = ttlind
        if (diags) write(*,*) "tevent = ", tevent
c-------Vary supertree and supergrass parameters in first year. 

        if ( (crtyr .le. 1) .and. .not. warmstart_init_flag ) then
          call veg_change(vveg_in, vveg, mx, c3c4, mx_out, 
     .                c3c4_out,tmpi,tmpi_out,ppti,ppti_out)
        endif 


        treename = "SUPRT"
        cropname = "SUPRG"

c        write(*,*) 'mx_index = ',mx_index,' percent evergreen.'
c        write(*,*) 'c3c4_index = ',c3c4_index,' percent C3.'

c-------C. Daly 6-16-97.For testing purposes, read an external file with 
c-------a tree and crop name to be inserted when a disturbance
c-------occurs, or if its the first year.  For years leading up
c-------to the first disturbance, file 'tree1' contains the
c-------treename.  For post-disturbance years, file 'tree2'
c-------contains the tree name.  If an error occurs 
c-------in the read (files don't exist, bad names), go on as normal.

cx        if (disturb .eq. .false.) then
cx          open(unit=2, file='tree1',status='OLD',err=95)
cx          read(2, 90, err=95) new_tn1
cx90        format(a5)
cx          close(2)
cx          treename = new_tn1
cx          write(*,*) 'Pre-disturbance treename hardwired to = ',treename
cx          go to 94

cx95        write(*,*)'Unsucessful read of tree1 file.Treetype unchanged.'
cx        endif

cx94      if (disturb .eq. .true.) then
cx          open(unit=2, file='tree2',status='OLD',err=99)
cx          read(2, 90, err=99) new_tn2
cx          close(2)
cx          treename = new_tn2
cx          write(*,*)'Post-disturbance treename hardwired to = ',treename
cx          go to 98

cx99        write(*,*)'Unsucessful read of tree2 file. ',
cx     .               'Treetype unchanged.'
cx        endif

98      tname = treename
        cname = cropname

        ctemp = 0
        ttemp = 0

cx        if (diags) write(*,*) __LINE__, ' >> new year ',
cx     $        treename(1:5), ' ', cropname(1:5), month, time

c-------C. Daly 2-5-96.  Search through cmdary and find last incidence
c-------of given event type.  Overwrite there.  If no occurrence,
c-------add to end of event list.  Ordering of events does not seem
c-------to be critical.

c-------This is for trees.

        if (getlen(treename) .ne. 0) then
           if(wbegin.gt.12.or.wbegin.lt.1) wbegin = 1
           if(wend.gt.12.or.wend.lt.1) wend = 12
           call update_event('TREE', treename, 0, crtyr, month,tevent)
           ttemp = FORSYS
        endif 

c-------This is for crops.

        if (getlen(cropname) .ne. 0) then
            if(cbegin.gt.12.or.cbegin.lt.1) cbegin = 1
            if(cend.gt.12.or.cend.lt.1) cend = 12
            call update_event('CROP', cropname, 0, crtyr, month,tevent)
            ctemp = CRPSYS
        endif

        /**************
        There are some questions about wether the cursys should ever
        be either of the "pure" types, FORSYS or CROPSYS.  The feeling
        is that MAPSS will have __some__ grass in all but the thickest
        forest.

        C. Daly (note) Since ctemp = 1 if grass is present and ttemp = 2
        if forest is present, adding them together will give a cursys of
        1 (CROPSYS) for grass only, 2 (FORSYS) for forest only, and 3 
        for (SAVSYS) both grass and forest.
        ***************/

        prevsys = cursys
        cursys = ctemp + ttemp

        /************
        The following chunk of code is copied from the bottom of
        readblk()

c       Check if decsys needs to be updated.  If it does, predec()
c       will also have to be called so new decomposition parameters
c       are initialized.
         *************/
c-------C. Daly (note) - adding woody compartment to the system...
        if (decsys .eq. CRPSYS .and. ttemp .eq. FORSYS) then
            decsys = FORSYS
            call predec(sand)
        endif

      endif	/* if (month .eq. 1) then */

c-----Fire events
      wfireLen = getlen(wfire)
      if (wfireLen .ne. 0) then
         if (diags) write(*,11) 'wfireLen, wfire = ', wfireLen, wfire
         /********
         Some sort of woody fire is scheduled for the current month.
         Due to the way the models interact, this really means that
         the previous month was ripe for a fire.
         *********/
         call update_event('TREM', wfire, 0, crtyr, month, tevent)
      endif
      cfireLen = getlen(cfire)
      if (cfireLen .ne. 0) then
         if (diags) write(*,11) 'cfireLen, cfire = ', cfireLen, cfire
         /********
         Some sort of crop fire is scheduled for the current month.
         Due to the way the models interact, this really means that
         the previous month was ripe for a fire.
         *********/
         call update_event('FIRE', cfire, 0, crtyr, month, tevent)
      endif
 11   format(A18, I4, A10)
c-----C. Daly 5-30-97.  At end of year,
c-----Search event array for either a grass
c-----burn (FIRE) or a tree removal event (TREM).  If either of
c-----these occur, green light is given for veg type to change
c-----the following January.

      if (month .eq. 12) then

c-------C. Daly 1-5-98.  Change mix indexes and veg type right
c-------here, instead of waiting until the following January.
c-------This is necessary in case transient run is seeded the
c-------very next year.

c-------C. Daly 2-8-98.  Changing veg type every year, now,
c-------so took it out of the if(fire) loop.

        call veg_change(vveg_in, vveg, mx, c3c4, mx_out, 
     .                c3c4_out,tmpi,tmpi_out,ppti,ppti_out)
	 
        do 20 i=1,tevent
          if(cmdary(i) .eq. 'TREM' .or. cmdary(i) .eq. 'FIRE') then
            if (diags) then
             write(*,*)
             write(*,*) 'Disturbance occurred in year ',crtyr,
     .                  ' month ', timary(i, 2), '.'
             write(*,*) 'i, tevent ='
             write(*,12) i, tevent
 12          format(2I4)
             write(*,*) 'Years since last burn: ',burn_count
             write(*,*)
             do 19 j = 1,tevent
              write(*,*) j, cmdary(j)
 19          continue            
            endif
            vegtype_burn = .true.
            burn_out = .true.
            cen_burn_count = 0.
            go to 21
          endif
20      continue
        vegtype_burn = .false.
        burn_out = .false.
        cen_burn_count = cen_burn_count + 1
21    continue  
      burn_count = REAL(cen_burn_count)          
      endif

c...Perform annual tasks
      if (month .eq. 1) then
         call eachyr()
      endif

c...The main driver for the model; call decomposition, growth, etc.
      if (cursys .ne. 3) then
      print *, '*** csa_main.F line 777 cursys = ', cursys
      endif
      call simsom()

      nppx = fcacc + agcacc + bgcacc

c.....annual accumulation of production and runoff
 
      if ( month .eq. 12) then
         fcaccx = fcacc
	 afcaccx = afcacc
         agcaccx = agcacc
         bgcaccx = bgcacc
         rspx = resp(1) + resp(2)
         nepx = nppx - rspx 
         soilnmx = soilnm(1) 
      else
         fcaccx = 0.
         agcaccx = 0.
         bgcaccx = 0.
         rspx = 0.
         nepx = 0.
         soilnmx = 0.
      endif
 
 
c.....calculating new monthly vars from monthly vars
 
      aetx = evap + tran
      frstcx = frstc
      crootcx = crootc
      frootcx = frootc
      aglivcx = aglivc
      bglivcx = bglivc
      vegc = frstc + aglivc + bglivc
      tslc = somtc + som1c(1)+ strucc(1) + metabc(1) + 
     *      stdedc + wood1c + wood2c + wood3c - stream(5)
      adeadc = wood1c + wood2c + stdedc + strucc(1) + metabc(1) 
      bdeadc = wood3c + strucc(2) + metabc(2) + som1c(2) + 
     *      som2c(1) + som3c(1) - stream(5)
     
c... August 24, 1999 DB
c... checking how minx has been calculated: we are missing a few terms

c      minx = metmnr(1,1) + metmnr(2,1) + w1mnr(1)
c      minx = minx + w2mnr(1) + w3mnr(1) + strmnr(1,1)

c... new definition identical to century variable tnetmn

      minx = metmnr(1,1) + metmnr(2,1) + 
     *       strmnr(1,1) + strmnr(2,1) +
     *       w1mnr(1) + w2mnr(1) + w3mnr(1) +
     *       s1mnr(1,1) + s1mnr(2,1) +
     *       s2mnr(1) + s3mnr(1)
 
      litf = lleaf + lfwood + llwood + lfroot + lcroot + lgrass
      clsl = stream(5)

c--------------------------------------------------------------
c...END OF VEMAP OUTPUT CALCULATIONS
c--------------------------------------------------------------

c...Update time
      time = time + dt
      if (time .ge. -1.0e-07 .and. time .le. 1.0e-07) then
         time = 0.0
      endif

c-----C. Daly - reactivate output interval dtpl, read from schedule file
c...Write out values
      if (tplt - time .lt. dt * 0.5) then
         if (gobin) call wrtbin(time)
         tplt = time + dtpl
      endif

c-----C. Daly 9-10-97.  Set burn flag for output.
c      burn = float(did_burn)
c...  new code from Ed 4-27-98

       if (did_burn) then
		burn = 1.0
       else 
		burn = 0.0
       endif

c-----Collect output vars for pass back to MAPSS.
      call varsout(outvars)
      if (diags) write(*,*) 'minerl(5,1) = ', minerl(5,1)
      call stateout(state)
c      /***************
c      These are the assignments that actually pass the tree leaf
c      carbon and the crop (grass) leaf carbon back to MAPSS.
c      Other assignments that pass data back into MAPSS should be
c      done here, like all the stuff required for Jim's super hot
c      fire model.
c      ****************/

      treec = rleavc
      cropc = aglivc

cx      if (diags) write(*,*) 'rleavc = ',rleavc
cx      if (diags) write(*,*) 'aglivc = ',aglivc

c-----C. Daly (note) This is where month is incremented
      month = mod(month,12) + 1

cx      if (diags) write(*,*) __LINE__, evtptr, ttlind


      end		/****	return to MAPSS		****/



c--------------------------------------------------------------------------
      subroutine cen_end

      /*********************
      Do a little clean up and end this Century session.
      **********************/

c...Close data files

c...Close the weather file
      close(unit=9)
c...Close the c14data file
      close(unit=10)
c...Close the schedule file
      close(unit=15)

c...Mark end of file
      endfile(unit=1)
 
c...Close binary file
      close(unit=1)

      end
      
      
	subroutine cen_zero()
c	Write zeroes into all the labeled commons.	

	include 'chrvar.inc'
	include 'comput.inc'
	include 'const.inc'  /* /save_vars/ */
	include 'dovars.inc'
	include 'fertil.inc'
	include 'flow.inc'
c                /zzzflo/	
	include 'forrem.inc'
	include 'isovar.inc'
	include 'ligvar.inc'
	include 'mapssvar.inc'
c	         /mapssvar/, /mapsschr/
	include 'param.inc'
	include 'parcp.inc'
	include 'parfs.inc'
	include 'parfx.inc'
	include 'pass.inc'
	include 'pheno.inc'	
	include 'plot1.inc'
	include 'plot2.inc'
	include 'plot3.inc'
	include 'plot4.inc'
	include 'potent.inc'
        include 'schvar.inc'
	include 'seq.inc'
	include 'site.inc'
	include 'supergrass.inc'
c                /suprg/	
	include 'supertree.inc'
c                /suprt/	
	include 't0par.inc'
	include 'timvar.inc'
	include 'vemap2.inc'
	include 'wth.inc'
	include 'zztim.inc'
	
	call write_zeroes(chrvarStart, chrvarEnd)
	call write_zeroes(computStart, computEnd)
	call write_zeroes(save_varsStart, save_varsEnd)
	call write_zeroes(dovarsStart, dovarsEnd)
	call write_zeroes(fertilStart, fertilEnd)
	call write_zeroes(zzzfloStart, zzzfloEnd)
	call write_zeroes(forremStart, forremEnd)
	call write_zeroes(isovarStart, isovarEnd)
	call write_zeroes(ligvarStart, ligvarEnd)
	call write_zeroes(mapssvarStart, mapssvarEnd)	
	call write_zeroes(mapsschrStart, mapsschrEnd)
	call write_zeroes(paramStart, paramEnd)
	call write_zeroes(parcpStart, parcpEnd)
	call write_zeroes(parfsStart, parfsEnd)
	call write_zeroes(parfxStart, parfxEnd)
	call write_zeroes(passStart, passEnd)
	call write_zeroes(phenoStart, phenoEnd)	
	call write_zeroes(plot1Start, plot1End)
	call write_zeroes(plot2Start, plot2End)
	call write_zeroes(plot3Start, plot3End)
	call write_zeroes(plot4Start, plot4End)
	call write_zeroes(potentStart, potentEnd)
        call write_zeroes(schvarStart, schvarEnd)
	call write_zeroes(seqStart, seqEnd)
	call write_zeroes(siteStart, siteEnd)
	call write_zeroes(suprgStart, suprgEnd)
	call write_zeroes(suprtStart, suprtEnd)
	call write_zeroes(t0parStart, t0parEnd)
	call write_zeroes(timvarStart, timvarEnd)
	call write_zeroes(vemap2Start, vemap2End)
	call write_zeroes(wthStart, wthEnd)
	call write_zeroes(zztimStart, zztimEnd)
	
	return
	end

