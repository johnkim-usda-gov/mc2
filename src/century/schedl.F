c*****************************************************************
c	$RCSfile: schedl.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: schedl.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/08/30 16:51:25  conklind
cImplement /save_vars/.
c
cRevision 1.1  2004/08/30 15:59:57  conklind
cInitial revision
c
c******************************************************************

c...  SCHEDL.F
c-----Eliminated references to pass1 -- Becky 4/7/96 (CD)
c-----C. Daly 5-3-96
c-----Got rid of ifdefs by deleting chunks of code and adding diags,
c-----a global flag that is 1 for write and 0 for no write to diags.

      subroutine schedl()

c...  Determine the next set of scheduling options from the 
c     schedule file

      include 'chrvar.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'schvar.inc'
      include 'seq.inc'
      include 'timvar.inc'
      include 'zztim.inc'
      
c-----C. Daly 5-3-96.  Added global flags.
      include 'mapssvar.inc'

c...  Local variables
      integer   crtyr, plntd
      real      savedfert
      character*5 curcult, curfert, curfire, curgraz, curharv, 
     $     curirri, curomad, curtrm
c Equivalence persistent local variables to /save_vars/.     
      equivalence (plntd, plntd_schedl), (curcult, curcult_schedl),
     - (curfert, curfert_schedl), (curfire, curfire_schedl), 
     - (curgraz, curgraz_schedl), (curharv, curharv_schedl),
     - (curirri, curirri_schedl), (curomad, curomad_schedl),
     - (curtrm, curtrm_schedl), (savedfert, savedfert_schedl)

 100  format(a5)


c...  Check crop dolast to reset crpgrw; done here so that crop grows
c     through the last month of growth
      if (dolast) then
         crpgrw = 0
         msplt = 0
         dolast = .false.
      endif

c...  Check forest doflst to reset forgrw; done here so that forest grows
c     through the last month of growth
      if (doflst) then
         forgrw = 0
         doflst = .false.
      endif

c...  Check if months since planting (msplt) needs to be updated
      if (plntd .eq. 1. and .stemp .ge. rtdtmp) then
         msplt = msplt + 1
      endif

c...  Reset do variables to false
      docult = .false.
      doerod = .false.
      dofert = .false.
      doflst = .false.
      dofone = .false.
      dofrst = .false.
      dograz = .false.
      dohrvt = .false.
      doirri = .false.
      dolast = .false.
      doomad = .false.
      doplnt = .false.
      dosene = .false.
      dotrem = .false.
      dofire(CRPSYS) = .false.
      dofire(FORSYS) = .false.
      dofire(SAVSYS) = .false.
      aufert = 0.0
      harmth = 0
c     evntyp = 3

c...  Convert time to integer year
      crtyr = aint(time + .001)
      crtyr = mod((crtyr - strtyr + 1), rptyrs)

c...  Working with real numbers - inexact so check on small number -rm
      if (crtyr .lt. 0.1) then
         crtyr = rptyrs
      endif

c...  Determine if evtptr needs to go back to 1
 10   continue

      if (ttlind .ne. 1 .and. evtptr .gt. ttlind) then
         evtptr = 1
      endif

c...  Look for events in timary that match the current time
c     If found, handle the event

      evtptr = 1

c-----Begin event loop---------------------------------------

 9010 continue

cx      if (evtptr .lt. ttlind) then
c-----J. Lenihan - 3-22-96.  Missing last event, so set to .le.
      if (evtptr .le. ttlind) then

         if ((timary(evtptr,1) .eq. crtyr) .and.
     $        timary(evtptr,2) .eq. month) then

            if (cmdary(evtptr) .eq. 'CROP') then

c--------------C. Daly 5-3-96 new write flag
cx               if( diags )  write(*,*) __LINE__, ' CROP called ', 
cx     .                      typary(evtptr)

               if ((typary(evtptr) .ne. 'F') .and. 
     $              (curcrp .ne. typary(evtptr))) then
                  /*************
                  It seems that pass1 is always .true., so this is just
                  when the crop type changes.  This means that every
                  time then the crop type changes, the crop.100 file is
                  opened and reread.  Furture versions of the linked
                  models should store the values of the crop.100 file in
                  an array/structure so that it is read only once.
                  **************/
                  call cropin(typary(evtptr))
                  call co2eff(time)
               endif
            endif
            if (cmdary(evtptr) .eq. 'PLTM') then
               doplnt = .true.
               seedl = 1
               plntd = 1
               msplt = 0
               crpgrw = 1
               falprc = 0
               prcfal = 0
            endif
            if (cmdary(evtptr) .eq. 'HARV') then
               dohrvt = .true.
               plntd = 0
               falprc = 1
               prcfal = 0
               harmth = 1
               if (curharv .ne. typary(evtptr)) then
c                  call harvin(typary(evtptr),curharv)
         stop 'Unexpected code for MC1'
               endif 
            endif
            if (cmdary(evtptr) .eq. 'FRST') then
cx               if( diags ) write(*,*) __LINE__, ' FRST called'
               dofrst = .true.
               crpgrw = 1
            endif
            if (cmdary(evtptr) .eq. 'LAST') then
c--------------C. Daly 5-3-96 new write flag
cx               if( diags ) write(*,*) __LINE__, ' LAST called'
               dolast = .true.
            endif
            if (cmdary(evtptr) .eq. 'SENM') then
c--------------C. Daly 5-3-96 new write flag
cx               if( diags ) write(*,*) __LINE__, ' SENM called'
               dosene = .true.
            endif
            if (cmdary(evtptr) .eq. 'FERT') then
               dofert = .true.
               aufert = savedfert
               if (curfert .ne. typary(evtptr)) then
         stop 'Unexpected code for MC1'
c                  call fertin(typary(evtptr),curfert,savedfert)
               endif
            endif
            if (cmdary(evtptr) .eq. 'CULT') then
               docult = .true.
               if (curcult .ne. typary(evtptr)) then
                   stop 'Unexpected code for MC1'            
c                  call cultin(typary(evtptr),curcult)
               endif
            endif
            if (cmdary(evtptr) .eq. 'OMAD') then
               doomad = .true.
               if (curomad .ne. typary(evtptr)) then
c                  call omadin(typary(evtptr),curomad)
         stop 'Unexpected code for MC1'
                  
               endif
            endif
            if (cmdary(evtptr).eq. 'IRRI') then
               doirri = .true.
               if (curirri .ne. typary(evtptr)) then
c                  call irrgin(typary(evtptr),curirri)
         stop 'Unexpected code for MC1'
               endif
            endif
            if (cmdary(evtptr) .eq. 'GRAZ') then
               dograz = .true.
               /*************
               It seems that pass1 is always .true., so this is just
               when the graze style changes.  This means that every time
               then the graze style changes, the graz.100 file is opened
               and reread.  Furture versions of the linked models should
               store the values of the graz.100 file in an array/structure
               so that it is read only once.
               **************/
c--------------C. Daly 5-3-96 new write flag
cx               if( diags ) write(*,*) __LINE__, ' GRAZ called'
               if (curgraz .ne. typary(evtptr)) then
c                  call grazin(typary(evtptr),curgraz)
         stop 'Unexpected code for MC1'
               endif
            endif
            if (cmdary(evtptr).eq. 'EROD') then
               doerod = .true.
               psloss = fltary(evtptr, 1)
            endif
            if (cmdary(evtptr) .eq. 'FIRE') then

cx               dofire(cursys) = .true.
               dofire(cursys) = .true.

               /*************
               It seems that pass1 is always .true., so this is just
               when the fire type changes.  This means that every time
               then the fire type changes, the fire.100 file is opened
               and reread.  Furture versions of the linked models should
               store the values of the fire.100 file in an array/structure
               so that it is read only once.
               **************/

c--------------C. Daly 5-3-96 new write flag
cx               if( diags ) write(*,*) __LINE__, ' FIRE called'
               if (curfire .ne. typary(evtptr)) then
                  call firein(typary(evtptr),curfire)
               endif
            endif
            if (cmdary(evtptr) .eq. 'TREE' .and. 
     $              curtre .ne. typary(evtptr)) then
c--------------C. Daly 5-3-96 new write flag
cx               if( diags ) write(*,*) __LINE__, ' TREE called ',
cx     .                       typary(evtptr)

               /*************
               It seems that pass1 is always .true., so this is just
               when the tree type changes.  This means that every time
               then the tree type changes, the tree.100 file is opened
               and reread.  Furture versions of the linked models should
               store the values of the tree.100 file in an array/structure
               so that it is read only once.
               **************/
               call treein(typary(evtptr))
c              write(*,*) ' treein() called by schedl.F at time =', time
               call co2eff(time)
            endif
            if (cmdary(evtptr) .eq. 'TREM') then
cx               dotrem = .true.
               dotrem = .true.
               /*************
               It seems that pass1 is always .true., so this is just
               when the trem type changes.  This means that every time
               then the trem type changes, the trem.100 file is opened
               and reread.  Furture versions of the linked models should
               store the values of the trem.100 file in an array/structure
               so that it is read only once.
               **************/
c--------------C. Daly 5-3-96 new write flag
cx               if( diags ) write(*,*) __LINE__, ' TREM called'
               if (curtrm .ne. typary(evtptr)) then
                  call tremin(typary(evtptr),curtrm)
               endif
            endif
            if (cmdary(evtptr) .eq. 'TFST') then
c--------------C. Daly 5-3-96 new write flag
cx               if( diags ) write(*,*) __LINE__, ' TFST called'
               dofone = .true.
               forgrw = 1
            endif
            if (cmdary(evtptr) .eq. 'TLST') then
c--------------C. Daly 5-3-96 new write flag
cx               if( diags ) write(*,*) __LINE__, ' TLST called'
               doflst = .true.
            endif
         endif

         evtptr = evtptr + 1
         goto 9010

      endif

c-----C. Daly 5-3-96 new write flag
cx      if( diags ) write(*,*) __LINE__,
cx     $     timary(evtptr,1), crtyr, timary(evtptr,2), month

      /********
      For some reason, this has to be here.  Otherwise you will get
      errors from cpp about a missing endif.
      *********/

 99   continue

      end
