
c*****************************************************************
c	$RCSfile: csa_detiv.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: csa_detiv.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:04:00  bachelet
cheader
c
c******************************************************************

c		Copyright 1993 Colorado State University
c			All Rights Reserved



c****Need to pass latitude and soils****

      subroutine detiv()

c...Determine name of schedule file, which contains the
c     name of the site file, values of timing variables,
c     and order of events
 
      include 'wth.inc'
      include 'chrvar.inc'
      include 'param.inc'
      include 'seq.inc'
      include 't0par.inc'
      include 'timvar.inc'

      include 'mapssvar.inc'

c...Local variables

      integer clen, getlen 
      character*132 newbin, oldbin, schnam
      logical ext, goahead
      character*132 filename

c      /***********
c      Now that Ive made these fixed length strings, I note that (Sun)
c      FORTRAN has a syntax to allow me to use variable length strings
c      as arguments to functions.  It should, since the caller has
c      the responsibility of allocation.  I should change to that
c      syntax.  See getlen() for an example.
c      ************/

      if (diags) write(*,*) 'years_to_run = ', years_to_run
      if (diags) write(*,*) 'ename = ', ename
      if (diags) write(*,*) 'nname = ', nname
      if (diags) write(*,*) 'sname = ', sname

c...Initialize weather labels - Becky 4-7-96 (CD)
      wlabel(1) = 'prec'
      wlabel(2) = 'tmin'
      wlabel(3) = 'tmax'

c...Initialize potential command line arguments
      ext = .false.
      schnam = ' '
      newbin = ' '
      oldbin = ' '

c...VAX NOTE: Need to get information interactively from user

c----Check for existence of interval file.  If none, don't do bin file.
      filename = ""
      filename(1:pathlen) = filpath
      filename(pathlen+1:pathlen+8) = 'interval'
c      inquire(file = filename, exist = gobin)

      if (gobin) then
         stop 'Unexpected code for MC1'
         
c-----Get extension bin file only in one-step mode.
cx      if (onestep) then

cx        clen = getlen(ename)
cx        oldbin = ename(1:clen)
cx        oldbin(clen + 1 : clen + 4) = '.bin'
cx        inquire(file = cx'/data/liveoak1/gclim_tools/models/century/model/'//oldbin, exist = goahead)

cx        clen = getlen(ename)
cx        oldbin(1:pathlen) = filpath
cx        oldbin(pathlen+1:pathlen+clen) = ename(1:clen)
cx        oldbin(pathlen+clen+1:pathlen+clen+4) = '.bin'

cx        inquire(file = oldbin, exist = goahead)
cx        write(*,*) oldbin(1:pathlen+clen+4)

cx        ext = .true.
cx        if (.not. goahead) then
cx           call message(' ')
cx           call message('   The old binary file could not be read.')
cx           call message(' ')
cx           STOP 'Execution error.'
cx        endif

cx      endif /* onestep */

cx        clen = getlen(nname)
cx        newbin = nname(1:clen)
cx        newbin(clen + 1 : clen + 4) = '.bin'
cx      inquire(file = cx'/data/liveoak1/gclim_tools/models/century/model/'//newbin, exist = goahead)

      clen = getlen(nname)
      newbin(1:pathlen) = filpath
      newbin(pathlen+1:pathlen+clen) = nname(1:clen)
      newbin(pathlen+clen+1:pathlen+clen+4) = '.bin'

c      inquire(file = newbin, exist = goahead)
      if (diags) write(*,*) newbin(1:pathlen+clen+4)

c-----C. Daly - 2/28/96.  I think Jesse had this backwards.
cx      if (.not. goahead) then
      if (goahead) then
         stop 'Unexpected code for MC1'
         
         call message(' ')
         call message('   The new binary file already exists.')
         call message(' ')
         call message('   So I will over-write it.')
      endif

      endif /* if (gobin) */

cx      clen = getlen(sname)
cx      schnam = sname(1:clen)
cx      schnam(clen + 1 : clen + 4) = '.sch'
cx      inquire(file = '/data/liveoak1/gclim_tools/models/century/model/'//schnam, exist = goahead)
cx      write(*,*) schnam(1:clen + 4)

      clen = getlen(sname)
      schnam(1:pathlen) = filpath
      schnam(pathlen+1:pathlen+clen) = sname(1:clen)
      schnam(pathlen+clen+1:pathlen+clen+4) = '.sch'

      inquire(file = schnam, exist = goahead)
      if (diags) write(*,*) schnam(1:pathlen+clen+4)

      if (.not. goahead) then
         call message(' ')
         call message('   The schedule file could not be read.')
         call message(schnam)
         call message(' ')
         STOP 'Execution error.'
      endif

c...Check that minimal information was entered
      if (schnam .eq. ' ') then
         call message(' ')
         call message('   No schedule file name was given.')
         call message(' ')
         STOP 'Execution error.'
      endif

      if (gobin) then

      if (newbin .eq. ' ') then
         if (ext) then
            newbin = oldbin
         else
            call message(' ')
            call message('   No binary output file name was given.')
            call message(' ')
            STOP 'Execution error.'
         endif
      endif

      endif

c...Open binary file to write to
c-----Only open if interval file exists.  Otherwise, assume no bin writing.

        if (gobin) then
          if (ext .and. newbin .eq. oldbin) then
            open(unit=1,file=newbin,form='UNFORMATTED',status='OLD')
          else
            open(unit=1,file=newbin,form='UNFORMATTED',status='UNKNOWN')
          endif
        endif

c...Open the schedule file and read the header lines

cx      open( unit = 15, file="/data/liveoak1/gclim_tools/models/century/model/"//schnam , status = 'OLD')

      open( unit = 15, file = schnam, status = 'OLD')

c  /* Starting year */
      read(15,915) strtyr
 915  format (I4)  
      
      read(15,*) tend /* Last year */
      tend = tend + 1

      read(15,*) sitnam /* Site file name */

c-----Read equilibrium check cycle here, instead of labeling type.
cx      read(15,*) labtyp /* Labeling type */
      labtyp = 0
      read(15,*) eqyears

      read(15,*) labyr  /* Labeling year */

      read(15,*) mctemp /* Microcosm */
      micosm = 0
      if (mctemp .ge. 0) then
         micosm = 1
      endif

      read(15,*) co2sys /* CO2 Systems */
      if (co2sys .gt. 0) then
         read(15,*) co2tm(1), co2tm(2)
      endif

      read(15,*) decsys /* Initial system */
      if (decsys .eq. 3) then
c        /**************
c        Whoa dawgies!!!  What is going on here???  I had quite a good
c        discussion with Dominique, Dennis, and Becky about this.  It
c        comes down to this being (in Beckys own words) a historical
c        idiosyncrasy.  The cursys is set from the combination of
c        the CROP and TREE entries in the current block.  The decsys
c        variable could (should) be set from that.
c
c        Though it would be better to use the symbolic names FORSYS and
c        SAVSYS for this strange assignment, those names are not
c        known at this point.
c        ***************/
        decsys = 2
      endif

      read(15,*) initcp /* Initial crop */
40    format(a5)
      if (initcp .eq. 'Initi') then
        initcp = ' '
      endif

      read(15,*) initre /* Initial tree */
      if (initre .eq. 'Initi') then
        initre = ' '
      endif

      read(15,*) /* nothing, a blank line */
      read(15,*) /* Year Month Option */

c...Read starting values from fixed parameter file
      call fixin()
 
c...Read starting values from site-specific file
c-----Becky 4-7-96 (CD)
c      open(unit=7,file=sitnam,status='OLD')

c-----Only do site file read if in standard mode.
c-----In onestep mode, site vars have already been read.
cx    if(onestep .eq. .false.) then
cx        open(unit=7,file='/data/liveoak1/gclim_tools/models/century/model/'//sitnam,status='OLD',err=1000)

        filename = ""
        filename(1:pathlen) = filpath
	clen = getlen(sitnam)
        filename(pathlen+1:pathlen+clen) = sitnam(1:clen)
        open(unit=7, file=filename,status='OLD',err=1000)
 
cx        open(unit=7,file=sitnam,status='OLD',err=1000)

        call sitein(ext)
cx      endif

c----Becky 4-7-96 (CD)
c...Moved the read calls for the initial tree and crop to inside the extend
c     if statement.  This is done to prevent a rather subtle bug that occurs
c     when the initial crop/tree do not match the final values in the 
c     original schedule.  In that case, the derived output values, (crpval ...)
c     do not match the current crop values.
c   The crop/tree reads must occur before the calciv call on a normal run.
c    7/20/95  K. Killian

c...Obtain the initial values for the crop or forest system

c...Determine initial values
      if (ext) then
        if (oldbin .ne. newbin) then
         stop 'Unexpected code for MC1'            
c        open(unit=3,file=oldbin,form='UNFORMATTED',status='OLD')
cx       call extend(3,strtyr,.TRUE.)
c          close(unit=3)
        else
cx          call extend(1,strtyr,.FALSE.)
        endif
c...Obtain the crop or forest system values (execute a CROP/TREE event)
        if (initcp .ne. ' ') call cropin(initcp)
        if (initre .ne. ' ') call treein(initre)
c        write (*,*) 'csa_detive.F called treein.'        

        call calciv
      else
c...Obtain the initial values for the crop or forest system
        if (initcp .ne. ' ') call cropin(initcp)
        if (initre .ne. ' ') call treein(initre)
c        write (*,*) 'csa_detive.F called treein.'
        call calciv
      endif

c...Sum up isotopes
c	  print *, 'calling sumcar from csa_detiv.F'
      call sumcar()

c...Do preliminary initializations and calculations
      call prelim()

c...Read the first block of events

      call readblk()

      if (onestep) then
        if (diags) call message('   From detiv: initialized.')
      else
        if (diags) call message(' ')
        if (diags) call message('   Model is running...')
      endif

      return
c----Becky 4-7-96 (CD)
1000    call message(' Fatal error: unknown site file :'//sitnam)
        stop ' Abnormal Termination'
      end



      /***********
      This fella is causing me problems finding the length of NULL
      terminated C strings.
      ************/
cx      integer function getlen(name)
cx      character*(*) name

cx      integer rlen, max

cx      max = len(name)
cx      do 10 rlen = 1, max
cx         if (name(rlen:rlen) .eq. 0 .or. 
cx     .       name(rlen:rlen) .eq. ' ')  goto 20
cx 10   continue
cx      rlen = 0
 
cx 20   if (rlen .eq. 0) then
cx         getlen = max
cx      else
cx         getlen = rlen - 1
cx      endif

cx      return
cx      end


cx      integer function getlen(name)
cx      character*(*) name
cx      integer j
C
C -----------------------------------------------------------------------------
C     this subroutine left justifies the file name and determines the length
C
C Variables
C      Input
C   name    character (*)  the input and processed file name
C
C  Modified by K. Killian 8/11/94
C              included the left justification on a subroutine coded by Laura
C
C -----------------------------------------------------------------------------
 
cx15    getlen = index(name,' ')-1

cx      if (getlen .eq. -1) then
cx        getlen = len(name)
cx      else if (getlen .eq. 0) then
cx        do 20 j= 1,len(name)
cx#if defined(MAPSS_Cen) || defined(MC_STANDARD)
cx          if (name(j:j) .ne. 0) then
cx#else /* MAPSS_Cen or MC_STANDARD */
cx          if (name(j:j) .ne. ' ') then
cx#endif /* MAPSS_Cen or MC_STANDARD */
cx            name = name(j:)
cx            goto 15
cx          endif
cx20      continue
cx        getlen = 0
cx      endif
 
cx      return
cx      end
