c*****************************************************************
c	$RCSfile: parcp.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: parcp.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:47:56  conklind
cInitial revision
c
c******************************************************************
c23456789012345678921234567893123456789412345678951234567896123456789712
       common/parcp/astrec(3),aglrem,astgc,astlbl,astlig,
     - auirri,awhc,basfc2,bglrem,bioflg,biok5,biomax,clteff(4),
     - crpgrw,crprtf(3),cultra(7),efrgrn(3),
     - fallrt,fawhc,fdfrem(2),fdgrem,
     - fecf(3),feclig,ffcret,flfrem,flghrv,flgrem,fligni(2,2),fnue(2),
     - fret(3),frtc(3),frtsh,fsdeth(4),fulcan,gfcret,gret(3),grwprc,
     - grzeff,hibg,himax,himon(2),hiwsf,irramt,irraut,msplt,pltmrf,
     - pramn(3,2),pramx(3,2),prbmn(3,2),prbmx(3,2),rdr,rmvstr,rtdtmp,
     - sdethc,seedl,twhc,vlossp,mf_flfrem,mf_fdfrem(2),
     - parcpEnd

       integer auirri,bioflg,crpgrw,flghrv,grzeff,himon,msplt,seedl

       real astrec,aglrem,astgc,astlbl
     - ,astlig,awhc,basfc2,bglrem,biok5,biomax,clteff,crprtf,cultra
     - ,efrgrn,fallrt,fawhc,fdfrem,fdgrem,fecf,feclig
     - ,ffcret,flfrem,flgrem,fligni,fnue,fret,frtc,frtsh
     - ,fsdeth,fulcan,gfcret,gret,grwprc,hibg,himax,hiwsf,irramt,irraut
     - ,pltmrf,pramn,pramx,prbmn,prbmx
     - ,rdr,rmvstr,rtdtmp,sdethc,twhc,vlossp,mf_flfrem,mf_fdfrem

      real parcpStart, parcpEnd
      equivalence (parcpStart, astrec(1))
      save /parcp/
