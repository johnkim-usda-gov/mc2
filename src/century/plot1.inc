c*****************************************************************
c	$RCSfile: plot1.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: plot1.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.3  2004/08/30 16:43:06  conklind
cMake the spelling of plot1Start and plot1End consistent with other files.
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/19 19:34:09  conklind
cInitial revision
c
c******************************************************************

c               Copyright 1993 Colorado State University
c                       All Rights Reserved
c Dimensions of minerl have been changed from (10,3) to (11,2) to
c accomodate leaching into nlayer + 1 when nlayer = 10
       common/plot1/ agcacc,aglcis(2),aglivc,aglive(3),aminrl(3),
     - amt1c2,amt2c2,anerb,
     - as11c2,as21c2,as2c2,as3c2,asmos(10),ast1c2,ast2c2,
     - avh2o(3),bgcacc,bglcis(2),
     - bglivc,bglive(3),cgrain,cinput,clittr(2,2),
     - co2cce(2,2,3),co2crs(2),co2cpr(2),co2ctr(2),cproda,cprodc,
     - cprodf,creta,crmvst,crpstg(3),crpval,defac,dsomsc,
     - egrain(3),elimit,eprodc(3),eprodf(3),ermvst(3),
     - eupacc(3),eupaga(3),eupbga(3),evap,fertot(3),
     - harmth,hi,irract,irrtot,metabc(2),metabe(2,3),metcis(2,2),
     - minerl(11,2),spare(8),mt1c2(2),mt2c2(2),nfix,nfixac,occlud,
     - parent(3),pet,petann,plabil,prcann,ptagc,ptbgc,
     - pttr,rain,resp(2),relyld,rwcf(10),
     - s11c2(2),s21c2(2),s2c2(2),s3c2(2),satmac,sclosa,scloss,
     - sdrema,secndy(3),shrema,sirrac,snfxac(2),snlq,snow,soilnm(3),
     - somsc,somse(3),somtc,som1c(2),som1ci(2,2),som1e(2,3),som2c(1),
     - som2ci(2),som2e(3),som3c(1),som3ci(2),som3e(3),stdcis(2),
     - st1c2(2),st2c2(2),stemp,strcis(2,2),stream(8),
     - strlig(2),strucc(2),struce(2,3),sumnrs(3),stdedc,stdede(3),tave,
     - tminrl(3),tnetmn(3),totc,tran,volgma,volexa,volpla,wdfxaa,wdfxas,
     - plot1End

       real agcacc,aglcis,aglivc,aglive,aminrl
     -,amt1c2,amt2c2,anerb
     -,as11c2,as21c2,as2c2,as3c2,asmos,ast1c2,ast2c2
     -,avh2o,bgcacc,bglcis
     -,bglivc,bglive,cgrain,cinput,clittr
     -,co2cce,co2crs,co2cpr,co2ctr,cproda,cprodc
     -,cprodf,creta,crmvst,crpstg,crpval,defac,dsomsc
     -,egrain,elimit,eprodc,eprodf,ermvst
     -,eupacc,eupaga,eupbga,evap,fertot
     -,harmth,hi,irract,irrtot,metabc,metabe,metcis
     -,minerl,mt1c2,mt2c2,nfix,nfixac,occlud
     -,parent,pet,petann,plabil,prcann,ptagc,ptbgc
     -,pttr,rain,resp,relyld,rwcf
     -,s11c2,s21c2,s2c2,s3c2,satmac,sclosa
     -,scloss,sdrema,secndy,shrema,sirrac,snfxac,snlq,snow,soilnm
     -,somsc,somse,somtc,som1c,som1ci,som1e,som2c
     -,som2ci,som2e,som3c,som3ci,som3e,spare,stdcis
     -,st1c2,st2c2,stemp,strcis,stream,strlig
     -,strucc,struce,sumnrs,stdedc,stdede,tave
     -,tminrl,tnetmn,totc,tran,volgma,volexa,volpla,wdfxaa,wdfxas
           
       real plot1Start, plot1End
       equivalence (plot1Start, agcacc)
       save /plot1/
