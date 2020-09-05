c*****************************************************************
c	$RCSfile: mapssvar.inc,v $
c	$Revision: 1.1 $
c	$Date: 2006/10/13 18:25:50 $
c	$Locker: wellsj $
c	
c$Log: mapssvar.inc,v $
cRevision 1.1  2006/10/13 18:25:50  wellsj
cInitial revision
c
cRevision 1.3  2003/11/26 17:19:03  conklind
cPut co2val at the end and define mapssvarPart to facilitate saving only 306 values from
cthe beginning of the 5000 value co2val array.
c
cRevision 1.2  2003/11/24 15:39:06  conklind
cAdd symbols for determining length at compile time.
c
cRevision 1.1  2003/11/20 21:47:47  conklind
cInitial revision
c
c******************************************************************
c      This is the MAPSS set of global variables.
c
c
c-----C. Daly 5-2-97 Removed mix from common block.  Replaced by
c-----mx_index in plot4.inc.
c-----C. Daly 5-21-97 Added which_soils, so that model will default to
c-----site file soils values if no soils available.
c-----Also added fire mode for scheduling purposes.


c Note that the initial part of /mapssvar/, thru the first 307 values of co2val,
c are saved and restored by cen_save, cen_restore, for warmstarting.  Variables at 
c the end of /mapssvar/ after co2val are not saved for warmstart.  In particular, 
c unl_n_flag defaults to FALSE and must be set to TRUE with the -N 1 command
c line switch when mc1 is invoked.

      common /mapssvar/ ovveg, onestep, diags, ppt(12), tmaxi(12), 
     .             tmin(12), years_to_run, eqyears, pathlen, gobin,
     .             timstor(50,2), fltstor(50,2), ttlind_stor,
     .             which_soils, fire_mode, zone, co2val(5000),
     .             unl_n_flag, co2ramp_name, 
     .             mapssvarEnd     

      integer ovveg, years_to_run, eqyears, pathlen, timstor, 
     .        ttlind_stor, which_soils, fire_mode, zone
      logical onestep, diags, gobin, unl_n_flag
      real ppt, tmaxi, tmin, fltstor, co2val
      character*100 co2ramp_name

      real mapssvarStart, mapssvarEnd, mapssvarPart
      equivalence (mapssvarStart, ovveg)
      equivalence (mapssvarPart, co2val(307))
      save /mapssvar/

c-----ppt - monthly ppt 
c-----tmaxi - monthly tmaxi
c-----tmin - monthly tmin
c-----years_to_run - simulation years (one-step only)
c-----pathlen - length of century input and output file pathname

c-----ovveg - ??
c-----onestep - flag for onestep version.  If false, standard version.
c-----diags - flag for writing to diagnostics file.

c-----eqyears -- interval of years over which to check for som2c 
c-----EQ (std only)

c-----timstor -- one-step storage of timary from schedule file.
c-----fltstor -- one-step storage of fltary from schedule file.
c-----ttlind_stor -- one-step storage of ttlind from schedule file.

cx-----mix -- index for mixed tree types.

c-----LENIHAN 11/6/97 co2val -- yearly array of co2 ppm values for co2 ramp

c-----Daly 2-12-98.  zone - mapss temperature zone



      common /mapsschr/ ename, nname, sname, fname, filpath, 
     .                  cmdstor(50),typstor(50), rname,
     .                  wthrSave,
     .                  mapsschrPad3, mapsschrEnd

      character*15 ename, nname, sname, fname, typstor
      character*7 rname
      character*20 cmdstor
      character*100 filpath
      character*20 wthrSave

      real mapsschrStart, mapsschrEnd
      character*3 mapsschrPad3
      equivalence (mapsschrStart, ename)
      save /mapsschr/
 
c-----ename - extension binary filename (one-step only)
c-----nname - new binary filename
c-----sname - schedule filename
c-----fname - fix.100 filename
c-----rname - co2 ramp data filename
c-----filpath - pathname of century input and output files
c-----cmdstor -- one-step storage of cmdstor from schedule file.
c-----typstor -- one-step storage of typstor from schedule file.

