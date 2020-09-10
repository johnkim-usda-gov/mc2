
c*****************************************************************
c	$RCSfile: mnracc.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: mnracc.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:55:22  bachelet
cheader
c
c******************************************************************
C...MNRACC.F

      subroutine mnracc (mnrflo,gross,net)

      real      mnrflo, gross, net

c...Update mineralization accumulators.
c     written by vek 05/91

c     Input:
c         mnrflo = mineralization value returned by esched.
c                  A negative value indicates immobilization.
c      
c     Transput:
c         gross = gross mineralization (sums only mineralization, not
c                     immoblization)
c         net = net mineralization (mineralization - immobilization)
c     


c...Gross mineralization
      if (mnrflo .gt. 0.0) then
        gross = gross + mnrflo
      endif

c...Net mineralization
      net = net + mnrflo

      return
      end
