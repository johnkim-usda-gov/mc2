
c*****************************************************************
c	$RCSfile: killrt.F,v $
c	$Revision: 1.2 $
c	$Date: 2004/09/24 20:50:10 $
c	$Locker:  $
c	
c$Log: killrt.F,v $
cRevision 1.2  2004/09/24 20:50:10  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

c...KILLRT
      subroutine killrt(accum)

c...Death of roots due to cutting or fire in a forest.

c...Called from:	frem

      include 'const.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'zztim.inc'

c...Argument Declarations
      real	accum(ISOS)

c...Local Variables
      integer	iel
      real	crd, dethe, frc14, frd, recres(MAXIEL)

c...LENIHAN 5/28/97 

      if ( did_burn ) then
       fd(1) = mf_fd(1)
       fd(2) = mf_fd(2)
      endif

c...Death of FINE ROOTS
c...Bachelet 11-5-97 

      frd = frootc * fd(1)
      if (frootc .ne. 0.0) then	
        do 10 iel = 1, nelem
          recres(iel) = froote(iel) / frootc
 10     continue
      else
        do 15 iel = 1, nelem
          recres(iel) = 0.0	
 15     continue
      endif
      frc14 = frtcis(LABELD) / frootc
      call partit(frd,recres,2,frtcis,froote,wdlig(FROOT),frc14)

c...Death of COARSE ROOTS
c...Bachelet 11-5-97

      crd = crootc * fd(2)
      do 20 iel = 1, nelem
	if(crootc.ne. 0.0) then
       		dethe = crd * (croote(iel)/crootc)
	else
		dethe = 0.0
	endif
       	call flow(croote(iel),wood3e(iel),time,dethe)
 20   continue
      call csched(crd,crtcis(LABELD),crootc,
     $            crtcis(UNLABL),wd3cis(UNLABL),
     $             crtcis(LABELD),wd3cis(LABELD),
     $             1.0,accum)

      return
      end
