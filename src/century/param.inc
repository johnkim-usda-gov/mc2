c*****************************************************************
c	$RCSfile: param.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: param.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.3  2004/05/19 17:28:25  conklind
cDBs reconciliation.
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/19 23:22:47  conklind
cInitial revision
c
c******************************************************************

c               Copyright 1993 Colorado State University
c                       All Rights Reserved


       common/param/afiel(10),amov(10),awilt(10),basef,bulkd,
     - co2ipr(2),co2ice(2,2,3),co2irs(2),co2itr(2),co2sys,co2tm(2),
     - drain,epnfa(2),epnfs(2),falprc,
     - hpttr(12),htran(12),ivauto,labtyp,labyr,
     - maxtmp(12),mctemp,micosm,mintmp(12),
     - nelem,nlayer,nlaypg,ph,ppdf(4,2),prcskw(12),prcstd(12),prdx(3),
     - precip(12),psloss,pslsrb,rcelit(2,3),rces1(2,3),
     - rces2(3),rces3(3),remwsd,
     - satmos(2),satmt,sirri,snfxmx(2),sorpmx,stormf,
     - strm5l,strm5u,svuptk(3),swflag,totco2,trbasl,
     - mf_lmc(3),did_burn,annet,hdp,cdp,
     - paramEnd

       integer falprc,ivauto,labtyp,labyr,micosm,nelem
     - ,nlayer,nlaypg,swflag

       logical did_burn

       real afiel,amov,awilt,basef,bulkd
     - ,co2ipr,co2ice,co2irs,co2itr,co2sys,co2tm
     - ,drain,epnfa,epnfs
     - ,hpttr,htran,maxtmp,mctemp,mintmp
     - ,ph,ppdf,prcskw,prcstd,prdx,precip,psloss
     - ,pslsrb,rcelit,rces1,rces2,rces3,remwsd
     - ,satmos,satmt,sirri,snfxmx,sorpmx,stormf
     - ,strm5l,strm5u,svuptk,totco2,trbasl
     - ,mf_lmc,annet,hdp,cdp

       real paramStart, paramEnd
       equivalence (paramStart, afiel(1))
       save /param/
