
c*****************************************************************
c	$RCSfile: cropmix.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: cropmix.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:03:33  bachelet
cheader
c
c******************************************************************

      subroutine cropmix      

c-----Coded by C. Daly 8-1-97.  Routine is a companion to treemix.
c-----Calculates intermediate values of selected crop.100 parameters
c-----based on c3/c4_index.
 
c-----Here is present list of parameters changed:

c-----prdx1(1) and prdx2(1) - c3/c4 max production 
c-----ppdf1(1-4) and ppdf2(1-4) - c3/c4 temperature curve parameters
c-----pramn1 and pramn2 - min/max c/n, biomass = 0/max
c-----pramx1 and pramx2 - min/max c/n, biomass = 0/max
c-----prbmn1 and prbmn2 - intercept c/n as function of ppt
c-----prbmx1 and prbmx2 - slope c/n as function of ppt

      include 'chrvar.inc'
      include 'const.inc'
      include 'isovar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'site.inc'
      include 'plot4.inc'

      include 'mapssvar.inc'
      include 'supertree.inc'
      include 'supergrass.inc'

      external scale

c-----------------------------------------------------------------------
      integer switch
c-----If cropmix is to be switched off, set to 0; set to 1 to switch on.
c-----If switched off, default values of supergrass will be used.
      data switch /1/      

      real c3_fac
      integer i, j
c-----------------------------------------------------------------------

      if (switch.eq.1) then

              
        if (diags) write(*,*) 'Based on c3/c4 index of ',c3c4_index,',',  
     .           '  some crop.100 params changed. '

c-------Get mix index factor from original 0-100 range.
c-------c3_fac is level of c3-ness. 
        c3_fac = c3c4_index * 0.01

c-------Max production
          call scale (c3_fac, prdx1(1), prdx2(1), prdx(1))
cx          write(*,*) 'prdx(1) = ',prdx(1)

c-------Temp curves
        do 40 i = 1, 4
          call scale (c3_fac, ppdf1(i,1), ppdf2(i,1), ppdf(i,1))
cx          write(*,*) 'ppdf(',i,', 1) = ',ppdf(i,1)
40      continue

c-------C/N ratios
        do 60 i = 1, 2
          do 50 j = 1, MAXIEL
            call scale( c3_fac, pramn1(j,i), pramn2(j,i), pramn(j,i))
cx            write(*,*) 'pramn(',j,',',i,') = ',pramn(j,i)
50        continue
60      continue

        do 80 i = 1, 2
          do 70 j = 1, MAXIEL
            call scale( c3_fac, pramx1(j,i), pramx2(j,i), pramx(j,i))
cx            write(*,*) 'pramx(',j,',',i,') = ',pramx(j,i)
70        continue
80      continue

        do 100 i = 1, 2
          do 90 j = 1, MAXIEL
            call scale( c3_fac, prbmn1(j,i), prbmn2(j,i), prbmn(j,i))
cx            write(*,*) 'prbmn(',j,',',i,') = ',prbmn(j,i)
90        continue
100     continue

        do 120 i = 1, 2
          do 110 j = 1, MAXIEL
            call scale(c3_fac, prbmx1(j,i), prbmx2(j,i), prbmx(j,i))
cx            write(*,*) 'prbmx(',j,',',i,') = ',prbmx(j,i)
110       continue
120     continue


      endif
      
      return
      end
