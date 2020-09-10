c*****************************************************************
c	$RCSfile: plot4.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: plot4.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 17:37:27  conklind
cInitial revision
c
c******************************************************************
       common/plot4/ mx_index, c3c4_index, burn, surface_runoff, 
     .               tmp_index, ppt_index, symb, afcacc, 
     . 				       bgrema, bgrmai(2), bgrmae(3), 
     .               asmos_subsoil, plot4End
 
       real mx_index, c3c4_index, burn, surface_runoff, 
     .               tmp_index, ppt_index, symb, afcacc,
     .               bgrema, bgrmai, bgrmae,
     .               asmos_subsoil
           
       real plot4Start, plot4End
       equivalence (plot4Start, mx_index)
       save /plot4/
