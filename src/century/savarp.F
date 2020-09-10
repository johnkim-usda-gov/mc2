
c*****************************************************************
c	$RCSfile: savarp.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: savarp.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:39:23  bachelet
cheader
c
c******************************************************************
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine savarp

c...Compute variables for printing or plotting

      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

c...Function declarations
      real      del13out, fsfunc
      external  del13out, fsfunc

c...Local variables
      integer  iel, iso, lyr, m, nll
      real     fsol, temp


c...Calculate total non-living C, minimum over the year
      temp = totc
      totc = som1c(SOIL) + som1c(SRFC) + som2c(1) + som3c(1) +
     $       strucc(SOIL) + strucc(SRFC) +
     $       metabc(SOIL) + metabc(SRFC)
      totc = min(totc, temp)

c...Compute soil organic matter sums for plotting
      somsc = som1c(SOIL) + som2c(1) + som3c(1)
      somtc = somsc + strucc(SOIL) + metabc(SOIL)
      woodc = wood1c + wood2c + wood3c
      frstc = rleavc + frootc + fbrchc + rlwodc + crootc
c      print *, 'rleavc, frootc, fbrchc, rlwodc, crootc', 
c     $ rleavc, frootc, fbrchc, rlwodc, crootc
      fsysc = somtc + woodc + frstc + strucc(SRFC) + metabc(SRFC) +
     $        som1c(SRFC)

      do 10 iel = 1, nelem
        somse(iel) = som1e(SOIL,iel) + som2e(iel) + som3e(iel)
        somte(iel) = somse(iel) + struce(SOIL,iel) + metabe(SOIL,iel)
	woode(iel) = wood1e(iel) + wood2e(iel) + wood3e(iel)
	frste(iel) = rleave(iel) + froote(iel) + fbrche(iel) +
     $               rlwode(iel) + croote(iel)
	fsyse(iel) = somte(iel) + woode(iel) + frste(iel) +
     $               struce(SRFC,iel) + metabe(SRFC,iel) +
     $               som1e(SRFC,iel)
10    continue

c...Compute soil organic matter sums by isotope
      do 20 iso = UNLABL, LABELD
          somsci(iso) = som1ci(SOIL,iso) + som2ci(iso) + som3ci(iso)
          somtci(iso) = somsci(iso) + strcis(SOIL,iso) +
     $                  metcis(SOIL,iso)

c...Add litter layer components, including surface som1
c     to get total organic matter including residue.  vek 08-91
          tomres(iso) = somtci(iso) + strcis(SRFC,iso) +
     $                  metcis(SRFC,iso) + som1ci(SRFC,iso)
20    continue

c...Sum the co2 loss and save
      totco2 = totco2 
     $         + mt1c2(UNLABL) + mt1c2(LABELD)
     $         + mt2c2(UNLABL) + mt2c2(LABELD)
     $         + s11c2(UNLABL) + s11c2(LABELD)
     $         + s2c2(UNLABL)  + s2c2(LABELD)
     $         + s3c2(UNLABL)  + s3c2(LABELD)
     $         + st1c2(UNLABL) + st1c2(LABELD)
     $         + s21c2(UNLABL) + s21c2(LABELD)
     $         + st2c2(UNLABL) + st2c2(LABELD)

c...Sum all state variables
c     Include som1c(SRFC) since it is not included in somtc.   vek 08$91
      totalc = somtc + strucc(SRFC) + metabc(SRFC) + som1c(SRFC)
     $         + aglivc +  stdedc + bglivc + csrsnk(UNLABL)
     $         + csrsnk(LABELD) + woodc + frstc + strm5u + strm5l
     $         + totco2

c...Calculate tminrl
      plabil = 0.0
      do 50 iel = 1, nelem
        tminrl(iel) = 0.0
        do 40 lyr = 1, nlayer
          if (minerl(lyr,iel).gt.0.0) then
            tminrl(iel) = tminrl(iel) + minerl(lyr,iel)
            if (iel .eq. P) then
              fsol = fsfunc(minerl(lyr,P), pslsrb, sorpmx)
              plabil = plabil + (minerl(lyr,iel) * fsol)
            endif
          endif
40      continue
50    continue

      nll = nlayer + 1
      do 60 iel = 1, nelem

c...Include som1e(1,iel) since it is not included in somte.  vek 08-91
        totale(iel) = tminrl(iel) + somte(iel) + struce(SRFC,iel) 
     $                + metabe(SRFC,iel) + som1e(SRFC,iel) + aglive(iel) 
     $                + stdede(iel) + bglive(iel) + esrsnk(iel) 
     $                + minerl(nll,iel) + parent(iel) + secndy(iel) 
     $                + woode(iel) + frste(iel) + stream(iel+5) 
     $                + crpstg(iel) + forstg(iel)
	if (iel .eq. P) then
          totale(iel) = totale(iel) + occlud
        endif
60    continue

c...Above and below ground live C/N ratio
      if (aglive(N) .gt. 0.0 .and. bglive(N) .gt. 0.0) then
         aglcn = aglivc/aglive(N)
         bglcn = bglivc/bglive(N)
      else
         aglcn = -999.0
         bglcn = -999.0
      endif

c...Overall c/n, c/p, and c/s ratios in soil organic matter
      do 70 iel = 1, nelem

c...	tcerat becomes a NaN when somte goes to zero. 
c...	added this to prevent the problem but demand should be high if no N
c...	so this may need some reworking
c...	however it looks like tcerat is never used anywhere else: BIZARRE
c...	D. Bachelet 10-8-97

	if(somte(iel) .ne. 0.0) then
        	tcerat(iel) = somtc/somte(iel)
	endif
70    continue

c...Average annual value of defac, the decomposition factor which
c     combines the effects of temperature and moisture
      adefac = 0.0
      do 80 m = 1, MONTHS
c...A negative value indicates that defac has not yet been calculated 
c     for month m
        if (defacm(m) .lt. 0.) then
          goto 90
        endif
        adefac = adefac + defacm(m)
80    continue
      m = 13
90    m = m-1
      if (m .gt. 0) then
        adefac = adefac/float(m)
      endif

c...Clittr output
      clittr(SRFC,UNLABL) = metcis(SRFC,UNLABL) + strcis(SRFC,UNLABL)
      clittr(SRFC,LABELD) = metcis(SRFC,LABELD) + strcis(SRFC,LABELD)
      clittr(SOIL,UNLABL) = metcis(SOIL,UNLABL) + strcis(SOIL,UNLABL)
      clittr(SOIL,LABELD) = metcis(SOIL,LABELD) + strcis(SOIL,LABELD)

c...Delta 13C output
      dsomsc = del13out(somsci(LABELD), somsci(UNLABL), dsomsc)

      dsomtc = del13out(somtci(LABELD), somtci(UNLABL), dsomtc)

      dsom1c(SRFC) = del13out(som1ci(SRFC,LABELD), som1ci(SRFC,UNLABL),
     $                        dsom1c(SRFC))
      dsom1c(SOIL) = del13out(som1ci(SOIL,LABELD), som1ci(SOIL,UNLABL),
     $                        dsom1c(SOIL))
      dsom2c = del13out(som2ci(LABELD), som2ci(UNLABL), dsom2c)
      dsom3c = del13out(som3ci(LABELD), som3ci(UNLABL), dsom3c)

      dslit = del13out(clittr(SRFC,LABELD), clittr(SRFC,UNLABL), dslit)
      dblit = del13out(clittr(SOIL,LABELD), clittr(SOIL,UNLABL), dblit)

      dstruc(SRFC) = del13out(strcis(SRFC,LABELD), strcis(SRFC,UNLABL),
     $                        dstruc(SRFC))
      dstruc(SOIL) = del13out(strcis(SOIL,LABELD), strcis(SOIL,UNLABL),
     $                        dstruc(SOIL))
      dmetc(SRFC) = del13out(metcis(SRFC,LABELD), metcis(SRFC,UNLABL),
     $                        dmetc(SRFC))
      dmetc(SOIL) = del13out(metcis(SOIL,LABELD), metcis(SOIL,UNLABL),
     $                        dmetc(SOIL))
 
      return
      end

      real function del13out(labeled, unlabeled, oldval)
      real labeled, unlabeled, oldval

      include 'const.inc'

      if (unlabeled .ne. 0.0) then
        del13out = ((labeled/unlabeled)/PEEDEE - 1) * 1000
      else
	del13out = oldval
      endif
     
      return
      end
