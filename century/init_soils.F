
c*****************************************************************
c	$RCSfile: init_soils.F,v $
c	$Revision: 1.1 $
c	$Date: 2007/01/04 22:26:56 $
c	$Locker:  $
c	
c$Log: init_soils.F,v $
cRevision 1.1  2007/01/04 22:26:56  wellsj
cInitial revision
c
cRevision 1.1  2006/10/13 18:14:20  wellsj
cInitial revision
c
cRevision 1.2  2004/09/24 20:47:12  bachelet
cheader
c
c******************************************************************

      subroutine cen_init_soils(bd,snd,cly,depth,rock,ws,ndep1,ndep2)      

      include 'site.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'mapssvar.inc'

      real bd,ndep1,ndep2
      real snd, cly, depth, rock(3), depth_cm
      integer ws, n

      which_soils = ws

c-----If gridded soils data not available, just set bulk density
c-----and get out.

      bulkd = bd
      hdp = ndep1 * 1000.
      cdp = ndep2 * 1000.

      if (which_soils .eq. 1 .or. which_soils .eq. 2) then

        sand = snd * 0.01
        clay = cly * 0.01
        silt = 1.0 - (sand + clay)

c-------C. Daly 4-28-97 Calculating soil layering structure from depth.
c-------First 4 layers are 15 cm deep (to 60 cm), 30 cm beyond that.
c-------Max soil layers allowed is 9 (210 cm)

        depth_cm = depth * 0.1
        if (depth_cm .lt. 15.) then 
           depth_cm = 15.
        endif
        if (depth_cm .gt. 240. ) depth_cm = 240.
        if (depth_cm .le. 60. ) then
          nlayer = nint (depth_cm / 15.)
        else
          nlayer = 4 + nint ((depth_cm - 60.) / 30.)
        endif
        if (nlayer .gt. 10) nlayer = 10
		if (nlayer .lt. 1) nlayer = 1

c-------Rock fragments reduce effective soil depth of each layer.
c-------Assume data are for 3 layers 1: 0-50 cm, 2: 50-150 cm, 3: 150+
c-------This translates into century layers 1-4 use layer 1 rocks,
c-------century layers 5-7 use layer 2 rocks, layers 8-10 use layer 3 rocks.

        do 10 n=1,nlayer
          if (n .le. 4) then
            adep(n) = 15. * (1.0 - rock(1)*.01)
          else if (n .le. 7) then
            adep(n) = 30. * (1.0 - rock(2)*.01)
          else
            adep(n) = 30. * (1.0 - rock(3)*.01)
          endif
10      continue

      endif
 
      return
      end
