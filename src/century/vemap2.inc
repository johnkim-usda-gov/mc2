c*****************************************************************
c	$RCSfile: vemap2.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: vemap2.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.6  2005/05/25 16:23:11  lenihan
cchanges for output of production scalars
c
cRevision 1.5  2004/12/13 17:16:50  lenihan
cplaced tot_whc declaration before Vemap2End
c
cRevision 1.4  2004/12/03 18:32:51  lenihan
cadd tot_whc declaration
c
cRevision 1.3  2004/05/19 17:33:35  conklind
cDBs reconciliation.
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:49:36  conklind
cInitial revision
c
c******************************************************************

      common/vemap2/nppx,npp_prv,rspx,rsp_prv,nepx,
     - aetx,spare2,rnf_prv,vegc,tslc,minx,lleaf,lfroot,
     - lfwood,llwood,lcroot,lgrass,litf,clsl,evergw,
     - evergd, everg_index,needle_index,mapss_lai,
     - fcaccx, agcaccx, bgcaccx, frstcx, aglivcx, 
     - bglivcx, soilnmx, afcaccx, frootcx, crootcx,
     - adeadc, bdeadc, tot_whc, 
     - fprd_ppt, fprd_tmp, gprd_ppt, gprd_tmp,
     - mi_tree, mi_grass, bio_consume, ddecid, 
     - frost_index,
     - vemap2End
     
      real nppx,npp_prv,rspx,rsp_prv,nepx,aetx,spare2,
     - rnf_prv,vegc,tslc,minx,lleaf,lfroot,lfwood,
     - llwood,lcroot,lgrass,litf,clsl,evergw,evergd,
     - everg_index,needle_index,mapss_lai,
     - fcaccx, agcaccx, bgcaccx, frstcx, aglivcx,
     - bglivcx, soilnmx, afcaccx, frootcx, crootcx,
     - adeadc, bdeadc, tot_whc, fprd_ppt, fprd_tmp, 
     - gprd_ppt, gprd_tmp, mi_tree, mi_grass, 
     - bio_consume, ddecid, frost_index

      real vemap2Start, vemap2End
      equivalence (vemap2Start, nppx)
      save /vemap2/
