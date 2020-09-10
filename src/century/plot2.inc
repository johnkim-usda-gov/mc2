c*****************************************************************
c	$RCSfile: plot2.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: plot2.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 17:30:32  conklind
cInitial revision
c
c******************************************************************
c**** The variables in /plot2/ are written to the plot file in the UNIX
c**** version of the model, but not the Microsoft version.
c
      common /plot2/ accrst, adefac, agcisa(2), aglcn, 
     -    bgcisa(2), bglcn,
     -    cgracc, cisgra(2), cltfac(4), csrsnk(2),
     -    dblit, dmetc(2), dslit, dsom1c(2), dsom2c, dsom3c,
     -    dsomtc, dstruc(2),
     -    egracc(3), ereta(3),
     -    esrsnk(3),
     -    gromin(3),
     -    lhzcac, lhzeac(3),
     -    metmnr(2,3),
     -    prcfal, rmvsti(2), rnpml1,
     -    sdrmae(3), sdrmai(2), shrmai(2), 
     -    shrmae(3), somsci(2), somtci(2), somte(3), 
     -    strmnr(2,3), s1mnr(2,3), s2mnr(3), s3mnr(3),
     -    tcerat(3), tcnpro, tomres(2), totalc, totale(3),
     -    volex, volgm, volpl,
     -    wdfx, wdfxa, wdfxma, wdfxms, wdfxs,
     -    plot2End

      real accrst, adefac, agcisa, aglcn
     -    ,bgcisa, bglcn
     -    ,cgracc, cisgra, cltfac, csrsnk
     -    ,dblit, dmetc, dslit, dsom1c, dsom2c, dsom3c
     -    ,dsomtc, dstruc
     -    ,egracc, ereta
     -    ,esrsnk
     -    ,gromin
     -    ,lhzcac, lhzeac
     -    ,metmnr
     -    ,prcfal, rmvsti, rnpml1
     -    ,sdrmae, sdrmai, shrmai
     -    ,shrmae, somsci, somtci, somte
     -    ,strmnr, s1mnr, s2mnr, s3mnr
     -    ,tcerat, tcnpro, tomres, totalc, totale
     -    ,volex, volgm, volpl
     -    ,wdfx, wdfxa, wdfxma, wdfxms, wdfxs

      real plot2Start, plot2End
      equivalence (plot2Start, accrst)
      save /plot2/
