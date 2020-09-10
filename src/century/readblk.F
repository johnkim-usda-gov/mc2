c*****************************************************************
c	$RCSfile: readblk.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: readblk.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/08/30 16:50:08  conklind
cImplement /save_vars/.
c
cRevision 1.1  2004/08/30 15:56:33  conklind
cInitial revision
c
c******************************************************************

cx      subroutine readblk( years_to_run)
      subroutine readblk()

c...Reads the next block of events from the schedule file

c-----C. Daly 2-5-96 - passing years to run
c-----C. Got rid of ifdefs by using version flag (onestep)
c-----onestep = 0 for standard century, = 1 for one_step.

cx      integer years_to_run
cx      integer year_new, ttlind_orig


      include 'chrvar.inc'
      include 'const.inc' /* constants and /save_vars/ */
      include 'param.inc'
      include 'schvar.inc'
      include 'seq.inc'
      include 't0par.inc'
      include 'timvar.inc'
      include 'zztim.inc'
      include 'site.inc'

c-----C. Daly 5-3-96. Include new flags.
      include 'mapssvar.inc'


c...Local variables
      integer   blknum, ctemp, i, modt, ttemp, yrnum
      integer   prevsys
      equivalence (prevsys, prevsys_readblk) /* in /save_var/ */
      real      pltmo, strplt

c...Local variables for new input block
      integer ststrt, stend, readint, lch, istat
      real readreal
      character buffr*80, readstring*20
      character*132 filename

c...Set the starting year of the block
c...Using ANINT which rounds to the nearest whole number.  -rm
      strtyr = anint(time)

c...Read the block number, ending year, and # of years set up
      read(15,*) blknum /* Block # */
      read(15,*) blktnd /* Last year */
      read(15,*) rptyrs /* Repeats # years */
      read(15,*) strplt /* Output starting year */
      read(15,*) pltmo  /* Output month */
      read(15,10) dtpl  /* Output interval */
10    format(f6.3)

c******Resetting output interval****
c-----Read new interval only if file "interval" exists
c-----Otherwise, assume that binary file writing is not needed.

cx      if (blknum .eq. 1) then

cx        open(unit=2,file='/data/liveoak1/gclim_tools/models/century/model/'//'interval')
        filename = ""
        filename(1:pathlen) = filpath
        filename(pathlen+1:pathlen+8) = 'interval'
        inquire(file = filename, exist = gobin)

        if (gobin) then
          open( unit = 2, file = filename, status = 'OLD')
          read(2,*)  dtpl
          rewind(unit=2)
          if (diags) write(*,*) ' new output interval is ',dtpl
          close(unit=2)
        else
          dtpl = 100000.
cx          write(*,*) ' default output interval is ',dtpl
        endif

cx      endif

c-----C. Daly 2-5-96 Replace rptyrs with years to run, if one_step
      if (onestep) rptyrs = years_to_run

c...Add 1 more to blktnd for January of the ending year
      blktnd = blktnd + 1

c...Reset the dtpl value; the units are years, but it should
c     have a value representing an integral number of months as
c     accurately as possible.
      if (dtpl .lt. .01) then

c...Reset to monthly if lower than lowest limit
        dtpl = dt
      else

c...Rounding
        modt = int(12.*dtpl+.05)
        dtpl = real(modt)/12.
      endif

c...Set tplt, the next time to spit out output
      pltmo = pltmo - 1
      tplt = strplt + real(pltmo)/12.
      if (tplt-int(tplt) .eq. 0.0) then
        tplt = tplt + dtpl
      endif

c...Determine the weather type and call wthini
c     to initialize the weather for this block
      read(15,20) wthr /* Weather choice */
20    format(a1)
      if (wthr .eq. 'F') then
        read(15,30) wthnam
30      format(a20)
      endif

c      write(*,'(/a,i5)') ' reading block ', blknum
      call wthini(precip,prcstd,prcskw,mintmp,maxtmp)

c=======================================================================
c
c The following block reads the events from the schedule file.
c   Uses subroutines that allow free form number and string input.
c   It also simplifies the input format by relaxing the year input requirement
c   (it will assume the same year as the previous command)
c
c  History
c    Modified : 04/03/95 : Kendrick Killian
c              added error checking on Century events and dates
c    Written : 1/25/95 : Kendrick Killian 
c
c  New External References
c    character function READSTRING (UNIT, buffr, LCH)
c    real function      READREAL   (UNIT, buffr, LCH, ISTAT)
c    integer function   READINT    (UNIT, buffr, LCH, ISTAT)
c
c  Additional comments
c    1) input commands are converted to upper case 
c    2) input fields can be separated by commas, white space or endline 
c    3) UNIX like comments can be inserted in the input stream
c       - # is the comment character
c       - comments extend to the end of the line
c    4) Century events can be entered without the year.  The year is
c        pulled from the previous event.
c    5) Century commands with NO date arguments generate a fatal error
c
c========================== BEGIN MODIFIED CODE ========================

c...Read and save the events for this block
      ttlind = 0
c...set thebuffer pointer (lch) to beyond the buffer length to
c...   force the "get" routines to read an input line
      lch = 81
40    yrnum = readint(15,buffr,lch,istat)
c...check for the existance of the date arguments
      if (istat .le. 0) then
        write(*,*) ' fatal error: missing date information :',buffr
        stop
      endif
c
c...check for end of block
      if (yrnum .ne. -999) then
        ttlind = ttlind + 1
c...Fill in the date parameters

        timary(ttlind, 1) = yrnum
        timary(ttlind, 2) = readint(15,buffr,lch,istat)

c...check for the explicit month argument
        if (istat.eq.0) then
          if (ttlind.gt.1) then
            timary(ttlind, 1) = timary(ttlind-1, 1)
            timary(ttlind, 2) = yrnum
c...bad input
          else
            write(*,*) ' missing year designator :',buffr
            stop
          endif
        endif

c...check date input consistancy
c     use blknum to mark an error  (blknum < 0  bad input)
        if (timary(ttlind, 1) .le.0) then
          write(*,'(a,i5/a)') ' fatal error: year less than zero ',
     $          timary(ttlind, 1), buffr
          timary(ttlind, 1) = 1
          blknum = -abs(blknum)
        endif
c-------Kill off run only if std mode.  In onestep, repeat_event will
c-------catch it.
        if (timary(ttlind, 1) .gt.rptyrs) then
          if (onestep .eqv. .false.) then
            write(*,'(a,i5/a)') ' fatal error: year greater than '//
     $          ' rotation length ',timary(ttlind, 1), buffr
            timary(ttlind, 1) = rptyrs
            blknum = -abs(blknum)
          endif
        endif
        if ((timary(ttlind, 2).lt.0).or.(timary(ttlind, 2).gt.12)) then
          write(*,'(a,i5/a)') ' fatal error: illegal month ',
     $          timary(ttlind, 2), buffr
          timary(ttlind, 2) = 0
          blknum = -abs(blknum)
        endif

c...Read and record the command
        cmdary(ttlind)    = readstring(15,buffr,lch,ststrt,stend)
c
c...Fill in additional information   (if required)
c
        if (INDEX('CROP CULT FERT FIRE GRAZ HARV IRRI OMAD TREE TREM ',
     $      cmdary(ttlind)(1:4)) .ne. 0) then
c
c...get string options (types)
          typary(ttlind) = readstring(15,buffr,lch,ststrt,stend)
        elseif ('EROD' .eq. cmdary(ttlind)) then
c
c...get real number for erosion
          fltary(ttlind, 1) = readreal(15,buffr,lch,istat)
        else
c
c...check for valid commands
          if (INDEX('PLTM FRST LAST SENM TFST TLST ',
     $          cmdary(ttlind)(1:4)) .eq. 0) then
            write(*,'(2a/a)')' fatal error: unknown CENTURY command : ',
     $            cmdary(ttlind)(1:4), buffr
            cmdary(ttlind)(1:4) = 'err'
            blknum = -abs(blknum)
          endif
        endif


        if (onestep) then

c ***************************************************
        if (diags) write(*,'(a,i5,a,i5,i3,a,a4,3a,g12.5)') 'Read  
     $  event #',ttlind,'   ->', timary(ttlind, 1),timary(ttlind, 2),
     $  ' "',cmdary(ttlind),'" "',typary(ttlind),'"',
     $  fltary(ttlind, 1)
c ***************************************************

c---------C. Daly 2-7-95 - For one_step, use yrnum as a return interval.
c---------So, a yrnum of 4 means to repeat the event every 4th year, 
c---------starting with year 4; yrnum=1 means repeat event every year 
c---------starting with yr 1.

c---------C. Daly 5-21-97 - When Jim's fire mode is used, do not schedule
c---------fires from schedule file.

          if (fire_mode .ne. 0 .and. (cmdary(ttlind)(1:4) .eq. 
     .        'TREM' .or.
     .        cmdary(ttlind)(1:4) .eq. 'FIRE') ) then
            if (diags) write(*,*) 'External fire mod -Event ignored.'
          else
            call store_event()
          endif

        endif

        goto 40
      endif /* if (yrnum .ne. -999) */

c
c...abort on illegal event input
      if (blknum .lt. 0) stop
c============================ End Modified Code ========================

c...Reset evtptr, the array index into cmdary, timary, typary
      evtptr = 1

c...Set up cursys, the current system(s) in use if 
c     this is the first block read
c     cursys: CRPSYS = crop/grass; FORSYS = forest; SAVSYS = savanna
c     Reset system for each block to handle changing systems.
      
c...Store the value of the previous system
      if (blknum .gt. 1) then
         prevsys = cursys
      endif

c...Initialize the temporary variables
      ctemp = 0
      ttemp = 0
      do 80 i = 1, ttlind
         if (cmdary(i) .eq. 'CROP') then
            ctemp = CRPSYS
         else if (cmdary(i) .eq. 'TREE') then
            ttemp = FORSYS
         endif
 80   continue
      cursys = ctemp + ttemp

c...If no crop or forest designation was given in the
c     first block, use the initial system designated
      if (cursys .eq. 0) then
         if (blknum .eq. 1) then
            if (initcp .ne. ' ') then
               ctemp = CRPSYS
            endif
            if (initre .ne. ' ') then
               ttemp = FORSYS
            endif
            cursys = ctemp + ttemp
         else

c...If no crop or forest designation was given in this
c     block and it is not the first, use the previous
c     block's system designation.
            cursys = prevsys
         endif
      endif

c...Check if decsys needs to be updated.  If it does, predec()
c    will also have to be called so new decomposition parameters
c    are initialized    
      if (decsys .eq. CRPSYS .and. ttemp .eq. FORSYS) then
         decsys = FORSYS
         call predec(sand)
      endif

      return
      end


C 3456789012345678901234567890123456789012345678901234567890123456789012 4567890
      CHARACTER*20 FUNCTION readstring (UNIT,buffr,lch,ststrt,stend)
C
C=======================================================================
C
C  Description
C    Parses character data separated by white spaces or hard delimiters.
C    see comments below.  Null record results in zero parsed field
C
C  Arguments
C
C  - Input
C    buffr, character*(*), input string
C    LCH   , integer,       is the index of the first character to parse
C                           updated by the routine as parsing occurs
C    UNIT,   integer,       IO unit to read input lines from
C
C  - Output
C    ststrt, integer,       index of string start
C    stend,  integer,       index of string end
C
C  History
C    Modified : 03/31/95 : Kendrick Killian
C              Forced a soft terminator at the buffer end.  Corrects the
C              bug which causes right justified strings to be dropped
C    Modified : 01/25/95 : Kendrick Killian
C              Modified routine to locate and return a single field
C    Written : 3/15/91 : Kendrick Killian 
C
C  Error Conditions
C    None
C
C  External References
C       subroutine READCLIN   reads an input line
C
C  Additional comments
C    1) a data field is established/bounded by a pair of upstream and
C       downstream delimiters.
C    2) commas are considered as hard delimiters:
C       - hard delimiters delimits field with no exception.
C       - consequtive hard delimiters constitute a null field.
C    3) White spaces (blanks and tabs) are soft delimeters:
C       - soft delimiters that preceed or follow an established data 
C         field are ignored.
C       - consecutive soft delimeters do not constitute a null field.
C    4) a record is assumed to be bounded by hard delimiters.
C
C      Copyright 1992 - 95    Kendrick Killian     All Rights Reserved
C========================== BEGIN EXECUTABLE CODE ======================
C
      CHARACTER    buffr*(*), CHAR*1
      INTEGER      UNIT, ststrt, stend, lch

      INTEGER I2, I, LENLIN
      LOGICAL ISWHITE, ISHARD
C
C  Logical statement functions to defines delimeters
C
      ISWHITE(CHAR) = ((CHAR.EQ.' ') .OR. (CHAR.EQ.'	'))
      ISHARD(CHAR) = (CHAR.EQ.',')
C
C Set Initial Values
C
      readstring = " "
      I = lch - 1
C
C Set status (I2) flag to a hard delimiter condition
C    I2 condition codes are:
C        -1  Field terminated with a white space
C         0  Last Field terminated with a hard delimiter
C        >0  tracking through a token; location of last hard delimeter
C
      I2 = 0
C
C determine record specified length
C
      LENLIN = LEN(buffr)
C
C read a record if required
C

20    I= I+1
      IF (I.GT.LENLIN) THEN
        IF (I2.GT.0) THEN
          ststrt = I2
          stend = LENLIN
          GOTO 90
        ELSE
          CALL READCLIN(UNIT,buffr,lch)
          I = lch
        ENDIF
      ENDIF
C
C Check to see if we have encountered a hard delimeter
C
        IF (ISHARD(buffr(I:I))) THEN
C
C Check if the previous field was terminated with a hard terminator
C  If so this is a null field
C
          IF (I2 .EQ. 0 ) THEN
            ststrt = lch
            stend  = lch
C
C Check to see if we are stepping over a data field
C If so we need to terminate the field
C
          ELSE IF (I2 .GT. 0 ) THEN
            ststrt = I2
            stend = I-1
          ENDIF
C
C  Found a field
          GOTO 90
C
C Have we encountered a white character with an intervening field
C
        ELSE IF ((I2 .GT. 0) .AND. ISWHITE(buffr(I:I))) THEN
          ststrt = I2
          stend = I-1
          GOTO 90
C
C If the character is not white and we are starting a new field and
C we need to store its location  
C (NOTE: don't worry about hard terminators since they won't get here
C
        ELSE IF ((I2 .LE. 0) .AND. .NOT.ISWHITE(buffr(I:I))) THEN
          I2 = I
        END IF
C
C If we are here we are either stepping over white spaces after 
C terminating a field or we are stepping over a valid field
C  In either case    KEEP GOING
C
      GOTO 20
C
C Clean up and end processing
C
   90 CONTINUE
      lch = I+1
      readstring = buffr(ststrt:stend)
      RETURN
      END

C2345678901234567890123456789012345678901234567890123456789012345678901234567890
      REAL FUNCTION readreal(UNIT, buffr, LCH, ISTAT)
C
C=======================================================================
C
C  This is a general purpose routine for parsing a numeric character
C  string and converting it to a number. 
C
C  Arguments
C  - Input
C    buffr, character*(*), character string to be converted.
C    LCH  , INTEGER,       is the index of the first character to parse
C                          this is updated as parsing occurs
C    ISTAT, integer,       NOTE: if READREAL is called with ISTAT =-2 or -3,
C                                the previous real is returned.
C                                This allows access to the real decoded 
C                                in a failed integer read.
C  - Output
C    LCH  , INTEGER,    the index of the next token
C    ISTAT, integer,    conversion status
C                   < -1, Conversion error
C                      0, no numerical conversion (alpha field)
C                      1, valid integer (and real)
C                      2, valid floating point and integer approximation
C                      3, valid floating point and integer overflow
C
C  External References
C       subroutine READCLIN   reads an input line
C
C  History
C    Modified : 02/12/95 : Kendrick Killian
C              Substantial recoding to reduce the possibility of misreporting
C              unusual strings (strings starting with .,-,+,e)
C    Modified : 01/25/95 : Kendrick Killian
C              Modified arguments
C    Modified : 01/23/95 : Kendrick Killian
C              Removed FORTRAN77 non-compliant internal unformatted READ
C    Written  : 03/27/91 : Kendrick Killian
C
C      Copyright 1992 - 95    Kendrick Killian     All Rights Reserved
C========================== BEGIN EXECUTABLE CODE ======================
C
      CHARACTER buffr*(*),CHR*1
      INTEGER LCH, ISTAT, readint, UNIT
C
C  local variables
      LOGICAL NOSPAC, RDEC, INTEG, NEG, NEGE
      INTEGER I, J, LOCE, LENLIN, Cip, ND
      DOUBLE PRECISION Cfp
      SAVE Cfp
      INTEGER*2 Cex, Cdec
C
C  Overflow/underflow value for a 32 bit integer
      INTEGER ovrflw, Ndig
      PARAMETER (ovrflw = 2147483647, Ndig = 9)
C
C  Statement functions
      LOGICAL  ISUPC,ISLOWC,ISDIGIT,ISWHITE
C
C      ISUPC(CHR) = CHR.GE.'A' .AND. CHR.LE.'Z'
C      ISLOWC(CHR) = CHR.GE.'a' .AND. CHR.LE.'z'
      ISDIGIT(CHR) = CHR.GE.'0' .AND. CHR.LE.'9'
      ISWHITE(CHR)  = CHR.EQ.' ' .OR. CHR.EQ.'	'
C
      IF ((ISTAT .EQ. -2) .OR. (ISTAT .EQ. -3)) THEN
        readreal = Cfp
        ISTAT=0
        RETURN
      ENDIF
      readreal  = 0.
      RDEC = .TRUE.
      GOTO 10

      ENTRY readint(UNIT, buffr, LCH, ISTAT)
      RDEC = .FALSE.
      readint = 0
C
 10   LENLIN = LEN(buffr)
      IF (LCH.GT.LENLIN) CALL READCLIN(UNIT,buffr,LCH)
C
      I     = LCH-1
      J     = 0
      ND    = 0
      Cex   = 0
      Cip   = 0
      Cfp   = 0.
      Cdec  = 0
      LOCE  = 0
      ISTAT = 0
      NEG    = .FALSE.
      NEGE   = .FALSE.
      INTEG  = .TRUE.
      NOSPAC = .TRUE.
C
C check for numeric entries.  Increment I until a non numeric character
C
      I = LCH-1
20    I = I+1
      IF (I .LE. LENLIN) THEN
c     write (*,'(a,i3,3a,l2,f18.0,i10,3i4,3L3,2i4)')
c    1' character i=',I ,' /',buffr(I:I),'/',integ,
c    1   Cfp, Cip, ND, Cdec, Cex, NEG, NEGE, NOSPAC, J, loce
C   Accept a digit
        IF (ISDIGIT(buffr(I:I)) .AND. NOSPAC) THEN
          J = ICHAR(buffr(I:I)) - ICHAR('0')
          IF (LOCE.GT.0) THEN
C   record an exponent digit
            Cex = Cex * 10 + J
          ELSEIF (ND .GE. Ndig) THEN
C   Suppress excess precision  (record decimal point shift on integer)
            IF (INTEG) Cdec = Cdec+1
          ELSEIF (ND+J .EQ. 0)  THEN
C   Suppress leading zero's (record decimal shift if is a decimal
            ISTAT = 1
            IF (.NOT.INTEG) Cdec = Cdec-1
          ELSE
C   record a significant digit
            ISTAT = 1
            Cip = Cip * 10 + J
            ND = ND + 1
            IF (.NOT.INTEG) Cdec = Cdec-1
          ENDIF
          J=I
          NOSPAC = .TRUE.
C  Accept the first decimal point
        ELSE IF ((buffr(I:I).EQ.'.') .AND. INTEG .AND. NOSPAC) THEN
          INTEG = .FALSE.
          IF (J.EQ.0) J = I
C  Accept a plus sign
        ELSE IF ((buffr(I:I).EQ.'+') .AND. (J.LE.LOCE)) THEN
          IF (J.EQ.0) J = I
C  Accept a minus sign
        ELSE IF ((buffr(I:I).EQ.'-') .AND. (J.LE.LOCE)) THEN
          IF (J.EQ.0) J = I
          IF (LOCE.EQ.0) THEN
            NEG  = .TRUE.
          ELSE
            NEGE = .TRUE.
          ENDIF
C   Step over spaces
        ELSE IF (ISWHITE(buffr(I:I))) THEN
          NOSPAC = ((LOCE.ge.J).or.(J.eq.0))
C   Check for an E in an exponent field
        ELSE IF (LOCE.EQ.0 .AND. ISTAT.NE.0 .AND.
     1    ( (buffr(I:I).EQ.'E') .OR. (buffr(I:I).EQ.'e') .OR.
     2      (buffr(I:I).EQ.'D') .OR. (buffr(I:I).EQ.'d') ) ) THEN
          INTEG = .FALSE.
          NOSPAC = .TRUE.
          LOCE = I
        ELSE
c      WRITE (*,*) ' decoding',LCH,i,J,ND,LOCE,NOSPAC,ISTAT,buffr(I:I)
C         IF (I.GT.LENLIN .AND. J.EQ.0) GOTO 10
          GOTO 30
        ENDIF
C step to the next character
      GOTO 20
      ENDIF
C
C was anything accepted
30        LCH = I
          IF (ISTAT .EQ. 0) then
            IF (I.GT.LENLIN) GOTO 10
            IF(J.NE.0) LCH = J
            RETURN
          ENDIF
          LCH = I
C  update the buffer pointer for special cases
          IF (LOCE .GT. J) LCH = LOCE
          IF (buffr(I:I).EQ.',') LCH = I + 1
C
CC  if you have reached here start the conversion to a real number
C         start by combining the mantissa and the exponents
c      WRITE (*,*) ' parsed',Cfp,Cip,Cdec,Cex,ND,NEG,NEGE,ISTAT
      IF (NEGE) Cex = -Cex
      Cex = Cex+Cdec
C
c      Cfp =  DBLE(Cip) *10.d0**Cex
      Cfp =  DBLE(Cip) *10.d0**(Cex*1.0)
C
      IF (INTEG) THEN
         ISTAT = 1
      ELSE IF (Cfp.le.ovrflw) THEN
         Cip = Cfp
         ISTAT = 2
      ELSE
         Cip = ovrflw
         ISTAT = 3
      ENDIF
      IF (NEG) THEN
         Cfp = -Cfp
         Cip = -Cip
      ENDIF
C
      IF (RDEC) THEN
         readreal = Cfp
C        WRITE (*,*) ' --- Real Return --',readreal, LCH, ISTAT
      ELSE
         readint  = Cip
C        WRITE (*,*) ' --- Integer Return --',readint, LCH, ISTAT
      ENDIF
      RETURN
      END

      SUBROUTINE READCLIN(UNIT,buffr,LCH)
C
C=======================================================================
C
C  Description
C    Reads a character buffer from the specified unit.
C    The routine skips comment lines and upper cases alpha characters
C
C  - Input
C    UNIT,   integer,       IO unit to read input lines from
C
C  - Output
C    buffr, character*(*), input buffer
C    LCH   , integer,       is the index of the first character to parse
C
C  History
C    Modified : 03/14/95 : Kendrick Killian
C              added a 5 character over read to check for input overflow
C    Rewritten : 1/25/95 : Kendrick Killian 
C                          Recoded to remove special field processing
C
C  Error Conditions
C    None
C
C  External References
C    NONE
C
C  Additional comments
C    1) input commands are converted to upper case 
C    2) UNIX like comments can be inserted in the input stream
C       - # is the comment character
C       - comments extend to the end of the line
C
C      Copyright 1995   Colorado State University    All Rights Reserved
C=======================================================================
C
      CHARACTER buffr*(*), PAD*5, comment*1 
      INTEGER LCH, I, UNIT
      PARAMETER (comment='#')
C
10      READ(UNIT,'(a,a5)') buffr,PAD

C  Check for comments and buffer overflow
        I = INDEX(buffr,comment)

C  error if non comment data extends past the buffer dimension
        IF (I.EQ.0) THEN
          IF (PAD.NE.'     ') THEN
            WRITE(*,'(a/4a)') ' Fatal Error: Input line to long:',
     1                      buffr,"/",PAD,"/"
            STOP
          ENDIF
        ELSE

C  Remove comments
          IF(I .EQ. 1) GOTO 10
          buffr = buffr(:I-1)
        ENDIF

C  check for a blank line
        LCH = 0
 15     LCH = LCH +1
C  NOTE: the second string is a TAB  This does not comply with FORTRAN 77
C        standards.  If this is a problem: ASCII tab is ICHAR(9)
        IF(LCH .GT. LEN(buffr)) GOTO 10
        IF(buffr(LCH:LCH).EQ.' ' .OR. buffr(LCH:LCH).EQ.'  ') GOTO 15

C  upper case the input line
      DO 50 I=LCH,LEN(buffr)
        IF (buffr(I:I).GE.'a' .AND. buffr(I:I).LE.'z') buffr(I:I) = 
     1               CHAR(ICHAR(buffr(I:I)) + (ICHAR('A')-ICHAR('a')))
 50   CONTINUE
      RETURN
      END
