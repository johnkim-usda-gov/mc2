      
c*****************************************************************
c	$RCSfile: veg_change.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: veg_change.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 23:53:15  bachelet
cheader
c
c******************************************************************
	subroutine veg_change(vveg_in, vveg, mx, c3c4, mx_out, 
     .                      c3c4_out,tmpi,tmpi_out,ppti,ppti_out) 

      real mx, c3c4, mx_out, c3c4_out,tmpi,tmpi_out,ppti,ppti_out
      integer vveg, vveg_in

c-----C. Daly 1-6-98.  Made subroutine, so that it can be called
c-----in Jaunary of first year or EQ-seeded runs, or in December
c-----of burn years.

      include 'mapssvar.inc'
      include 'plot4.inc'
      include 'zztim.inc'	

      if (diags)
     .   write(*,*) 'veg type changed from ',vveg,' to ',
     .                       vveg_in, ' at time = ',time
      vveg = vveg_in

c-----C. Daly 8-1-97.  Now using "supertree" and "supergrass,"
c-----so treetype is always SUPRT and croptype is always SUPRG.
c-----However, treemix and cropmix now vary the values of several
c-----parameters based on the mixed indexes.  The mixed indexes
c-----are updated only after a fire.

c-----C. Daly 12-9-96. Update decid/evergreen mix index.
c-----C. Daly 5-2-97.  Changed name from mix to mx_index.
c-----C. Daly 8-1-97.  Update C3/C4 mix index.
c-----C. Daly 2-18-98. Now using tmp_index and ppt_index.
cx    mx_index = mx
      c3c4_index = c3c4
      tmp_index = tmpi
      ppt_index = ppti
c-----Pass back out to mapss for biomass->lai translation.
      mx_out = mx_index
      c3c4_out = c3c4_index
      tmpi_out = tmp_index
      ppti_out = ppt_index

c-----C. Daly 5-30-97.  Check to see if new tree name
c-----should be changed due to mix index.
      call treemix

c-----C. Daly 8-1-97.  Check to see if new crop name
c-----should be changed due to mix index.
      call cropmix

      return
      end

