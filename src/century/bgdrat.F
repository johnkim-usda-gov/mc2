      
c*****************************************************************
c	$RCSfile: bgdrat.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: bgdrat.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 19:59:55  bachelet
cheader
c
c******************************************************************
	real function bgdrat (aminrl,varat,iel)

c...BelowGround Decomposition RATio computation.

      integer   iel
      real      aminrl(3), varat(3,3)



c...Determine C/E of new material entering 'Box B'.

c...Input:
c     aminrl
c     varat : varat(1,iel) = maximum C/E ratio for material
c                              entering 'Box B'
c             varat(2,iel) = minimum C/E ratio for material
c             varat(3,iel) = amount of E present when minimum
c                              ratio applies
c     iel = index to element for which ratio is being computed;
c             1 for N, 2 for P, 3 for S

c...Output:
c     bgdrat = C/E ratio of new material where E is N, P, or S
c                depending on the value of iel



c...Ratio depends on available E 

      if (aminrl(iel) .le. 0.) then

c...Set ratio to maximum allowed
	  bgdrat = varat(1,iel)

      elseif (aminrl(iel) .gt. varat(3,iel)) then

c...Set ratio to minimum allowed
          bgdrat = varat(2,iel)

      else

c... aminrl(iel) > 0 and <= varat(3,iel)
         bgdrat = (1.-aminrl(iel)/varat(3,iel))
     $         *(varat(1,iel)-varat(2,iel))+varat(2,iel)

c...where   varat(1,iel) = maximum C/E ratio for material
c                            entering 'Box B'
c           varat(2,iel) = minimum C/E ratio for material
c           varat(3,iel) = amount of E present when minimum
c                            ratio applies

      endif

      return
      end
