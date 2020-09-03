
c*****************************************************************
c	$RCSfile: wrtbin.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: wrtbin.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:55:08  bachelet
cheader
c
c******************************************************************

c		Copyright 1993 Colorado State University
c			All Rights Reserved



      subroutine wrtbin(time)
      real time

      include 'outval.inc'

c**** Local variables 
      integer :: ierr = 0

c**** Write all output values to the binary file
c*****Write full set of one_step vars, even if they are not used
c-----(in mc_standard, for example), to maintain consistency

cx	return
cx      write(unit=1,iostat=ierr) time, vals1, vals2, vals3, vals4
cx      write(unit=1,iostat=ierr) time, vals1, vals2, vals3

c**** Check ierr for an error on writing
      if (ierr .ne. 0) then
         call message('   ***Error on writing to binary file.')
         call message('      Is the disk flooded?')
         STOP 'Execution error.'
      endif

      return
      end
