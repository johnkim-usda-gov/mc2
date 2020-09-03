      
c*****************************************************************
c	$RCSfile: ckdata.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: ckdata.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:01:16  bachelet
cheader
c
c******************************************************************
	subroutine ckdata(routin,expect,found)

      character*(*) routin, expect, found

c...Check a parameter name read by fixin or sitein.
c     routin is the name of the calling routine.
c     expect is the name that should have been read.
c     found is the name that was read from the data file.

c...Local variables
      integer i, loc
      character par*6, string*80


c...Extract junk from the parameter name

      loc = index(found,'(')
      if (loc .eq. 0) then
     	  loc = index(found,'|')
	  if (loc .eq. 0) then
	      loc = index(found,',')
          endif
      endif
      if (loc .eq. 0) then
      	par = found( :6)
      else
        par = found( :loc-1)
      endif

c...Convert to lower case if needed (MACHINE DEPENDANT)
      do 100 i = 1, 6
         if (par(i:i) .ge. 'A' .and. par(i:i) .le. 'Z') then
           par(i:i) = char( ichar(par(i:i)) + 32 )
         endif
100   continue

c...Test name for expected name
      if (expect .ne. par) then
	if (routin .eq. 'fixin') then
          call message('   There is an error in your fix.100 file.')
	else if (routin .eq. 'sitein') then
          call message('   There is an error in your <site>.100 file.')
	endif
        string = '   The data for ' // par // ' was read when the ' //
     +           'data for ' // expect // ' was expected.'
	call message(string)
	  STOP
      endif

      return
      end
