      
c*****************************************************************
c	$RCSfile: faterr.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: faterr.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:24:21  bachelet
cheader
c
c******************************************************************
	subroutine faterr(itell1,itell2)                                          

      character*4 itell1, itell2

c...An error on the weather data file has been detected;                      
c     terminate execution.                                                      
                                                                               
c...Local varialbes
      character*80 string


      call message('    ')
      call message('   The weather file is out of sequence.')
      string = '   A record for ' // itell2 // ' was expected.'
      call message(string)
      string = '   The record read was labeled ' // itell1 // '.'
      call message(string)

      STOP                                                                      
      end                                                                       
