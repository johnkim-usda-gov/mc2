
c*****************************************************************
c	$RCSfile: respir.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: respir.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:36:39  bachelet
cheader
c
c******************************************************************
C...RESPIR.F

      subroutine respir(co2los,nlr,lyr,tcstva,cstatv,csrsnk,resp,
     $                  estatv,minerl,gromin,netmnr)

      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'zztim.inc'

      integer   nlr, lyr
      real      co2los, tcstva(nlr), cstatv(nlr,ISOS), 
     $          csrsnk(ISOS),resp(ISOS), estatv(nlr,MAXIEL), 
     $          minerl(MAXLYR,MAXIEL),gromin(MAXIEL), 
     $          netmnr(nlr,MAXIEL)

c...Compute flows associated with microbial respiration.

c     Input:
c         co2los = CO2 loss associated with decomposition
c         nlr = total number of layers modeled for Box A;
c                 (=2 for structural, metabolic, som1;
c                  =1 for som2, som3, wood compartments)
c         lyr = layer for which respiration is being computed
c         tcstva(nlr) = total C (unlabeled + labeled) in layer 'lyr' of 
c                      Box A.  For components with only 1 layer,
c                      tcstva will be dimensioned (1).
c         time           passed in /zztim/
c     
c     Transput:
c         cstatv = C state variables for Box A;
c                  cstatv(lyr,iso) represents C in layer 'lyr',
c                  isotope 'iso' where iso=1 for unlabeled C
c                  and iso=2 for labeled C.  For components with
c                  only 1 layer, the 1st dimension of cstatv will
c                  be 1.
c         csrsnk = C source/sink 
c         resp = output variable
c         estatv = N, P, and S variables for Box A.  In 'estatv(lyr,iel)',
c                  lyr represents layer and iel indicates N, P, or S.  For
c                  components with only 1 layer, the 1st dimension of
c                  estatv will be 1.
c         minerl(1,iel) = labile N, P, or S in layer 1.
c         gromin(iel) = gross mineralization -- stored in /comput/
c         netmnr(lyr,iel)= net mineralization for layer lyr (N, P, or S)
c                         For components with only 1 layer, the 1st
c                         dimension of netmnr will be 1.


c...Local variables
      integer   iel
      real      mnrflo


c...C flow from cstatv to CO2
      if (labtyp .eq. 2) then
         stop 'Unexpected code for MC1'
        call csched(co2los,cstatv(lyr,LABELD),tcstva(lyr),
     $              cstatv(lyr,UNLABL),csrsnk(UNLABL),
     $              cstatv(lyr,LABELD),csrsnk(LABELD),
     $              dresp,resp)
       else
        call csched(co2los,cstatv(lyr,LABELD),tcstva(lyr),
     $              cstatv(lyr,UNLABL),csrsnk(UNLABL),
     $              cstatv(lyr,LABELD),csrsnk(LABELD),
     $              1.0,resp)
      endif
 
c...Mineralization associated with respiration
      do 10 iel = 1, nelem
          mnrflo = co2los*estatv(lyr,iel)/tcstva(lyr)
          call flow(estatv(lyr,iel),minerl(SRFC,iel),time,mnrflo)

c...Update gross mineralization
          gromin(iel) = gromin(iel) + mnrflo

c...Update net mineralization
          netmnr(lyr,iel) = netmnr(lyr,iel) + mnrflo

10    continue
 
      return
      end
