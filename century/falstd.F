      
c*****************************************************************
c	$RCSfile: falstd.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: falstd.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:23:53  bachelet
cheader
c
c******************************************************************
	subroutine falstd(pltlig)
      real pltlig(2)

c...Simulate fall of standing dead for the month.

      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot1.inc'

c...Local variables
      integer   iel
      real      fr14, fsdc, recres(MAXIEL)

      if (stdedc .gt. 0) then
        fsdc = stdedc * fallrt
        do 10 iel = 1, nelem
          recres(iel) = stdede(iel)/stdedc
 10     continue
        fr14 = stdcis(LABELD)/stdedc
        call partit(fsdc,recres,1,stdcis,stdede,pltlig(ABOVE),fr14)
      endif

      return
      end
