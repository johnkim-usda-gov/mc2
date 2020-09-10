
c*****************************************************************
c	$RCSfile: pschem.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: pschem.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:28:33  bachelet
cheader
c
c******************************************************************
C...PSCHEM.F

      subroutine pschem(dtm)

      real      dtm

c...This subroutine calculates the Phosophorus and Sulfur chemistry
c     - decomposition flows.  These calculations were removed from
c     the DECOMP subroutine, and slightly modified to include the
c     calculation for the fraction of mineral P in solution.  Also
c     removed the flow from secondary to parent material.
c
c     Date:        6/91
c     Coded by:        McKeown
c
c     Called From:        SIMSOM
c     Calls:                FLOW
c                             ATANF
c                        FSFUNC
c
c     Local Variables:        FMNSEC        - flow from mineral to secondary
c                        FPARNT        - flow from parent to mineral
c                        FSECND        - flow from secondary to mineral
c                        FSECOC        - flow from secondary to occluded P
c                        FSOL        - fraction of mineral P in solution
c
c    ********************************************************************

      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'site.inc'
      include 'zztim.inc'

c...Function declarations
      real      atanf, fsfunc
      external  atanf, fsfunc

c...Local variables
      integer   iel, lyr
      real      fmnsec, fparnt, fsecnd, fsecoc, fsol


c...For Phosphorus and Sulfur...
      do 20 iel = 2, nelem

c...Determine the fraction of mineral P in solution
        if (iel .eq. P) then
          fsol = fsfunc(minerl(SRFC,iel), pslsrb, sorpmx)
        else
          fsol = 1.0
        endif

c...Flow from parent material to mineral compartment.
c     Soil Texture may affect weathering of Phosophorus.

c...This calculation is actually done in prelim, it is
c     shown here for clarity.
c       if ((iel .eq. P) .and. (texepp(1) .eq. 1.0)) then
c         Include effect of texture
c         textur = clay + silt
c         Weathering factor should be per year
c         wfact = 12.0 * atanf(textur, texepp(2), texepp(3),
c    +                             texepp(4), texepp(5))
c       else
c         wfact = pparmn(iel)
c       endif

        fparnt = pparmn(iel) * parent(iel) * defac * dtm
        call flow(parent(iel), minerl(1,iel), time, fparnt)

c...Flow from secondary to mineral compartment.
c     Soil texture may affect mineralization of secondary P.

c...This calculation is actually done in prelim, it is
c     shown here for clarity.
c       if ((iel .eq. P) .and. (texesp(1) .eq. 1.0)) then
c         Include effect of texture
c         wfact = 12.0 * (texesp(2) + texesp(3) * sand)
c       else
c         wfact = psecmn(iel)
c       endif

        fsecnd = psecmn(iel) * secndy(iel) * defac * dtm
        call flow(secndy(iel), minerl(1,iel), time, fsecnd)

c...Flow from mineral to secondary
        do 10 lyr = 1, nlayer
          fmnsec = pmnsec(iel) * minerl(lyr,iel) * (1 - fsol)
     $                        * defac * dtm
          call flow(minerl(lyr,iel), secndy(iel), time, fmnsec)
10        continue
20    continue

c...Flow from secondary Phosophorus to occluded Phosophorus.
      fsecoc = psecoc * secndy(P) *defac * dtm
      call flow(secndy(P), occlud, time, fsecoc)
      return
      end
