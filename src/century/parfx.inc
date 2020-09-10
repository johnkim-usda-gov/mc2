c*****************************************************************
c	$RCSfile: parfx.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: parfx.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 19:37:23  conklind
cInitial revision
c
c******************************************************************

c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c-----C. Daly 4-7-96.  Had to add tmax, topt, tshl, and tshr back on
c-----to list to make updated code work.  Update was missing these vars.

      common/parfx/adep(10), agppa, agppb, aneref(3), animpt, 
     $    awtl(10), bgppa, bgppb, co2ppm(2), co2rmp,
     $    damr(2,3), damrmn(3),
     $    dec1(2), dec2(2), dec3(2), dec4, dec5, deck5, dligdf,
     $    dresp, edepth, elitst, enrich, favail(6), fleach(5), 
     $    fwloss(4),
     $    fxmca, fxmcb, fxmxs, fxnpb, gremb, idef, lhzf(3),
     $    minlch, nsnfix, ntspm,
     $    omlech(3), 
     $    p1co2a(2), p1co2b(2), p2co2, p3co2, pabres, pcemic(3,3),
     $    peftxa, peftxb, phesp(4),
     $    pligst(2), pmco2(2), pmnsec(3),
     $    pmntmp, pmxbio, pmxtmp, pparmn(3), pprpts(3),
     $    ps1co2(2), ps1s3(2), ps2s3(2),
     $    psecmn(3), psecoc, rad1p(3,3), rcestr(3),
     $    rictrl, riint, rsplig, seed,
     $    spl(2), strmax(2),
     $    texepp(5), texesp(3), teff(3), tmelt(2), tmax, topt, tshl, 
     $    tshr, varat1(3,3), varat2(3,3), varat3(3,3), vlosse, vlossg,
     $    parfxEnd

      integer idef, nsnfix, ntspm, seed

      real adep, agppa, agppb, aneref, animpt
     $    ,awtl, bgppa, bgppb, co2ppm, co2rmp
     $    ,damr, damrmn
     $    ,dec1, dec2, dec3, dec4, dec5, deck5, dligdf
     $    ,dresp, edepth, elitst, enrich, favail, fleach
     $    ,fwloss
     $    ,fxmca, fxmcb, fxmxs, fxnpb, gremb, lhzf
     $    ,minlch, omlech
     $    ,p1co2a, p1co2b, p2co2, p3co2, pabres, pcemic
     $    ,peftxa, peftxb, phesp
     $    ,pligst, pmco2, pmnsec
     $    ,pmntmp, pmxbio, pmxtmp, pparmn, pprpts
     $    ,ps1co2, ps1s3, ps2s3
     $    ,psecmn, psecoc, rad1p, rcestr
     $    ,rictrl, riint, rsplig
     $    ,spl, strmax
     $    ,texepp, texesp, teff, tmelt, tmax, topt, tshl 
     $    ,tshr,varat1, varat2, varat3, vlosse, vlossg

      real parfxStart, parfxEnd
      equivalence (parfxStart, adep(1))
      save /parfx/
