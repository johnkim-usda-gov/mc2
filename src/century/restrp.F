
c*****************************************************************
c	$RCSfile: restrp.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: restrp.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:37:06  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

C...RESTRP.F
c...Restrict the actual production based on C/E ratios.  Calculate
c     minimum, and maximum whole plant nutrient concentrations.

      subroutine restrp(nelem, nparts, avefrc, cerat, cfrac, potenc,
     $                  rimpct, storage, snfxmx, cprodl, eprodl,
     $                  eup, uptake, elimit, nfix, relyld,
     $                  snfxac, svuptk)

      include 'const.inc'
      include 'fertil.inc'
c mapssvar.inc is here because it contains unl_N_flag
      include 'mapssvar.inc'
      include 'parfx.inc'
c...schvar.inc is included here because it has the code_flags
      include 'schvar.inc'
      include 'seq.inc'
      include 'zztim.inc'

      integer	nelem, nparts
      real elimit
      real	avefrc(MAXIEL), cerat(2,nparts,MAXIEL), cfrac(nparts),
     $		cprodl, eprodl(MAXIEL), eup(nparts,MAXIEL),
     $          nfix, potenc, relyld, rimpct, snfxac, snfxmx,
     $		storage(MAXIEL), svuptk(MAXIEL), uptake(4,MAXIEL)

c...Local variables
c   NOTE:  Local variables cannot have adjustable array size.  MINECI
c          and MAXECI are set to the largest array size which may occur.
      integer	iel, ipart
      real	afert(MAXIEL), eavail(MAXIEL), minec(MAXIEL), 
     $		maxec(MAXIEL), mineci(5,MAXIEL),
     $          maxeci(5,MAXIEL), ustorg(MAXIEL), ediff
      logical fastrtnflag

c...Reset variables to zero
      cprodl = 0.0
      do 30 iel = 1, nelem
	eprodl(iel) = 0.0
	afert(iel) = 0.0
	do 20 ipart = 1, nparts
	  eup(ipart,iel) = 0.0
 20     continue
 30   continue

c...There is no production if one of the mineral elements in not
c     available.
      if ((snfxmx .eq. 0.0) .and. (aufert .eq. 0.0)) then
        do 40 iel = 1, nelem
          fastrtnflag = (avefrc(iel) .le. 0.0) 
          if ((iel .eq. N) .and. fastrtnflag) then
            fastrtnflag = .not. unl_N_flag
     $                    .or. code_flags(9)
          endif
          if (fastrtnflag) return
 40     continue
      endif

c...Initialize cprodl
      cprodl = potenc

c...Calculate soil available nutrients, adding storage.
      do 60 iel = 1, nelem
	eavail(iel) = favail(iel) * rimpct * avefrc(iel)
     $                  + storage(iel)
 60   continue

c...Calculate average E/C of whole plant (crop, grass, or tree)
      do 100 iel = 1, nelem
	minec(iel) = 0.0
	maxec(iel) = 0.0
	do 90 ipart = 1, nparts
	  mineci(ipart,iel) = 1 / cerat(IMAX,ipart,iel)
	  maxeci(ipart,iel) = 1 / cerat(IMIN,ipart,iel)
	  minec(iel) = minec(iel) + cfrac(ipart) * mineci(ipart,iel)
	  maxec(iel) = maxec(iel) + cfrac(ipart) * maxeci(ipart,iel)
 90     continue
 100  continue

c...Compute the limitation
      call nutrlm(nelem, nparts, cprodl, eprodl, maxec,
     $            maxeci, mineci, cfrac, eavail, nfix,
     $            snfxmx, snfxac, elimit, eup)


c...Calculate uptakes from all sources: storage, soil, nfix
      do 200 iel = 1, nelem
	ustorg(iel) = min(storage(iel), eprodl(iel))

c...If storage pool contains all needed for uptake
      if (eprodl(iel) .le. ustorg(iel)) then
	uptake(ESTOR,iel) = eprodl(iel)
	uptake(ESOIL,iel) = 0.0

c...Otherwise, extra necessary from the soil pool
        elseif (eprodl(iel) .gt. ustorg(iel)) then
	  uptake(ESTOR,iel) = storage(iel)
	  uptake(ESOIL,iel) = eprodl(iel) - storage(iel)

        endif

c...Save soil uptakes for use by savanna
        svuptk(iel) = uptake(ESOIL,iel)

 200  continue

c...N fixation uptake was computed in the limitation routines
      uptake(ENFIX,N) = nfix

c...Check for consistency
c
      ediff = eprodl(N) - uptake(ESTOR,N) - uptake(ESOIL,N)
     $           - uptake(ENFIX,N)
      if (ediff .gt. 0.000001) then
       call message('RESTRP not balancing N')
c      STOP
      endif


      return
      end
