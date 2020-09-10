/*
 *  gridspecs.h
 *  MC1
 *
 *
 */
#ifndef GRIDSPECS_H
#define	GRIDSPECS_H

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef Boolean
#define Boolean int
#endif
// #define TRUE 1
// #define FALSE 0 

/* Grid specs for the grids.  Row numbers increase from north to south.  Column numbers increase from west to east. */
/* Coordinates given here are of the northwest corner of the northwest grid cell. */
#define NOT_LAT_LON -360

#define US30asGRID_NAME "US30asToon5"
#define US30asLATITUDE0 49.0208333333333333333
#define US30asLONGITUDE0 -124.5208333
#define US30asCELL_SPACING 0.0416666666667 /* 150 arc-seconds, expressed as decimal degrees */
#define US30asNROWS_GRID 576
#define US30asNCOLS_GRID 1381

#define US4k_Toon6GRID_NAME "US4k_Toon6"
#define US4k_Toon6LATITUDE0 49.125  /* 49:07:30 */
#define US4k_Toon6LONGITUDE0 -124.875 /* 124:52:30 */
#define US4k_Toon6CELL_SPACING 0.25 /* 15' in decimal degrees */
#define US4k_Toon6NROWS_GRID 96
#define US4k_Toon6NCOLS_GRID 232

// GlobalToon5 - 2.5 degree by 2.5 degree (72 rows by 144 columns)
#define GlobalToon5GRID_NAME "GlobalToon5"
#define GlobalToon5LATITUDE0 90. /* 90Â° 00' N */
#define GlobalToon5LONGITUDE0 -180. /* 180Â° 00' W */
#define GlobalToon5CELL_SPACING 2.5 /* 2.5 decimal degrees */
#define GlobalToon5NROWS_GRID  72
#define GlobalToon5NCOLS_GRID 144

// GlobalToon3 - 1.5 degree by 1.5 degree (120 rows by 240 columns)
#define GlobalToon3GRID_NAME "GlobalToon3"
#define GlobalToon3LATITUDE0 90. /* 90° 00' N */
#define GlobalToon3LONGITUDE0 -180. /* 180° 00' W */
#define GlobalToon3CELL_SPACING 1.5 /* 1.5 decimal degrees */
#define GlobalToon3NROWS_GRID 120
#define GlobalToon3NCOLS_GRID 240

// Global - Half degree by half degree (360 rows by 720 columns)
#define GlobalGRID_NAME "Global"
#define GlobalLATITUDE0 90. /* 90° 00' N */
#define GlobalLONGITUDE0 -180. /* 180° 00' W */
#define GlobalCELL_SPACING 0.5 /* 30 arc-minutes, expressed as decimal degrees */
#define GlobalNROWS_GRID 360
#define GlobalNCOLS_GRID 720
// 360 x 720 = 259,200 cells total
// rows 0-12 are over water, and rows 292-359 are over the Southern Ocean and Antarctica
// Excluding Antarctica, there are ~61,000 live cells
// the subgrid consisting of rows 13-291 and cols 0-719 has 200,880 cells
// the northern hemisphere subgrid consisting of rows 13-179 has 50,228 active cells
// the 0,0 cell of the VEMAP grid is at row 82, col 111 in the Global grid (VEMAP grid is -r 82,119, -c 111,225)

// VEMAP - Conterminous United States. 48 rows, 115 cols, 5520 gridpoints, 3261 populated points
#define VEMAPGRID_NAME "VEMAP"
#define VEMAPLATITUDE0 49.0 /* 49° 00' N */
#define VEMAPLONGITUDE0 -124.5 /* 124° 30' W */
#define VEMAPCELL_SPACING 0.5 /* 30 arc-minutes, expressed as decimal degrees */
#define VEMAPNROWS_GRID 48
#define VEMAPNCOLS_GRID 115

// US12km - Conterminous United States. 192 rows, 460 cols, 88,320 gridpoints, ~52,000 populated points
// This definition is preliminary on 11/22/08.
#define US12kmGRID_NAME "US12km"
#define US12kmLATITUDE0 49.0 /* 49° 00' N */
#define US12kmLONGITUDE0 -124.5 /* 124° 30' W */
#define US12kmCELL_SPACING 0.125 /* 7.5 arc-minutes, expressed as decimal degrees */
#define US12kmNROWS_GRID 192
#define US12kmNCOLS_GRID 460

// NA8km - The Lynx grid - North America at 5 arc-minutes: 720 rows, 1392 cols, 1,002,240 gridpoints, 392,132 populated points
#define NA8kmGRID_NAME "NA8km"
#define NA8kmLATITUDE0 85.0 /* 85° 00' N */
#define NA8kmLONGITUDE0 -168.0 /* 168° 00' W */
#define NA8kmCELL_SPACING 0.08333333333 /* 5 arc-minutes, expressed as decimal degrees */
#define NA8kmNROWS_GRID 720
#define NA8kmNCOLS_GRID 1392
// bounding rectangles
// OR-WA HUC5 -r 431,523 -c 519,621 (ARRA project) 93 rows x 103 cols  9,579 cells total  49deg 05' 00" N, 124deg 45' 00" W

// US4km - Conterminous United States at 2.5 arc-minute resolution
// This grid is collinear with, and has the same bounding rectangle as, the US800m grid, 
// at 1/5 the resolution, and hence 1/25 as many gridcells.
#define US4kmGRID_NAME "US4km"
#define US4kmLATITUDE0 49.9375 /* 49° 56’ 15” N */
#define US4kmLONGITUDE0 -125.0208333 /* 125° 01’ 15” W */
#define US4kmCELL_SPACING 0.041666667 /* 2’ 30”, expressed as decimal degrees */
#define US4kmNROWS_GRID 621
#define US4kmNCOLS_GRID 1405
// southern boundary = 24° 03’ 45” N
// eastern boundary = 66° 28’ 45” W
// bounding rectangles
// OR-WA HUC5 -r 20,205 -c 7,209 186 rows x 203 cols = 37,758 cells, of which more or less 33,569 are active
//    north edge 49deg 06' 15" N
//    south edge 41deg 21' 15" N
//    west edge 124deg 43' 45" W
//    east edge 116deg 16' 15" W
// AZ-NM HUC5 -r 299,447 -c 242,544  149 rows x 303 cols = 45,147 cells
// SkyIslands -r 402,447 -c 320,396  46 rows x 77 cols = 3542 cells, of which ~2063 are active
// Klamath-Siskiyou (aka "OSW") -r 153,193 -c 10,74  41 rows x 65 cols = 2665 cells
// Cascade Head/HJA/Metolius (CHM) study area -r 117,137 -c 20,83  21 rows x 64 cols  1344 cells total
// WW2100 study area (Willamette River basin) -r 95,157 -c 30,80  63 rows x 51 cols  3213 cells total
//   with Sulzman's 3 cell buffer -r 92,160 -c 27,83  69 rows x 57 cols  3933 cells total
//   north edge = 46deg 06' 15"
//   south edge = 43deg 13' 45"
//   west edge = -123deg 53' 45"
//   east edge = -121deg 31 15"

// USeast 4km - Eastern United States
#define USeast4kmGRID_NAME "USeast4km"
#define USeast4kmLATITUDE0 50.47916667 /* 50° 28’ 45” N */
#define USeast4kmLONGITUDE0 -106.1458333 /* 106° 08’ 45” W */
#define USeast4kmCELL_SPACING 0.041666667 /* 2’ 30”, expressed as decimal degrees */
#define USeast4kmNROWS_GRID 651
#define USeast4kmNCOLS_GRID 975
// southern boundary = 23° 21’ 15” N
// eastern boundary = 65° 31’ 15” W

// US 800m - Conterminous United States
#define US800mGRID_NAME "US800m"
#define US800mLATITUDE0 49.9375 /* 49° 56’ 15” N */
#define US800mLONGITUDE0 -125.0208333 /* 125° 01’ 15” W */
#define US800mCELL_SPACING 0.008333333333 /* 30 arc-seconds, expressed as decimal degrees */
#define US800mNROWS_GRID 3105
#define US800mNCOLS_GRID 7025
// 3105 x 7025 = 21,812,625 cells total
// southern boundary = 24° 03’ 45” N
// eastern boundary = 66° 28’ 45” W
// bounding rectangles
// AZ-NM HUC5 -r 1498,2234 -c 1211,2718  737 rows x 1508 columns  1,111,396 cells total
// OR-WA HUC5 -r 104,1022 -c 35,1043
// Eastern Oregon study area (WWETAC VDDT project) -r 730,845 -c 380,519  116 rows x 140 cols
// Apache-Sitgreaves study area (WWETAC VDDT project) -r 1891,2007 -c 1848,1926  117 rows x 79 cols  9243 cells total
//    5311 live cells in 6HUC5 study area
// Cascade Head/HJA/Metolius (CHM) study area -r 585,689 -c 100,419  105 rows x 320 cols  33,600 cells total
// Washington Coast Range (WCR) project -r 185,460 -c 33,303   276 rows x 271 cols 74,796 cells in rectangle, about 40,000 active
// For just the Olympic ecoregion, use only rows -r 185,360, which gives about 24,000 active cells.
//    N edge 48deg 23' 45" N
//    S edge 46deg 05' 45" N
//    W edge 124deg 44' 45" W
//    E edge 122deg 29' 15" W
// Blue Mtns Ecoregion bounding box -r 429,769 -c 417,1068  341 rows x 652 cols = 222,332 cells
//    N edge 46deg 21' 45" N
//    S edge 43deg 31' 15" N
//    W edge 121deg 32' 45" W
//    E edge 116deg 06' 45" W

// US West 800m - Conterminous United States west of longitude 104° W
#define USwest800mGRID_NAME "USwest800m"
#define USwest800mLATITUDE0 49.9375 /* 49° 56’ 15” N */
#define USwest800mLONGITUDE0 -125.0208333 /* 125° 01’ 15” W */
#define USwest800mCELL_SPACING 0.008333333333 /* 30 arc-seconds, expressed as decimal degrees */
#define USwest800mNROWS_GRID 3105
#define USwest800mNCOLS_GRID 2522

// BearSoil 800m - Conterminous United States
#define BearSoil800mGRID_NAME "BearSoil800m"
#define BearSoil800mLATITUDE0 49.37916667 /* 49° 22’ 45” N */
#define BearSoil800mLONGITUDE0 -124.7375 /* 124° 44’ 15” W */
#define BearSoil800mCELL_SPACING 0.008333333333 /* 30 arc-seconds, expressed as decimal degrees */
#define BearSoil800mNROWS_GRID 2980
#define BearSoil800mNCOLS_GRID 6934
//
// BearSoilToon6 - used for 1:36 samples of 30 arc-second (800m) grids
#define BearSoilToon6GRID_NAME "BearSoilToon6"
#define BearSoilToon6LATITUDE0 49.025 /* 49° 01' 30” N */
#define BearSoilToon6LONGITUDE0 -124.525 /* 124° 31’ 30” W */
#define BearSoilToon6CELL_SPACING 0.05 /* 3 arc-minutes, expressed as decimal degrees */
#define BearSoilToon6NROWS_GRID 490
#define BearSoilToon6NCOLS_GRID 1152


// ATtest800m - test grid for Appalachian Trail project
#define ATtest800mGRID_NAME "ATtest800m"
#define ATtest800mLATITUDE0 38.77083333 /* 38° 46’ 15” N */
#define ATtest800mLONGITUDE0 -79.25416667 /* 79° 15’ 15” W */
#define ATtest800mCELL_SPACING 0.008333333333 /* 30 arc-seconds, expressed as decimal degrees */
#define ATtest800mNROWS_GRID 50
#define ATtest800mNCOLS_GRID 20

// California 800m - Used to be only from the Oregon border down to about the latitude of Mt. Whitney.
// Now the whole state.  CA800m 0,0 = US800m 75,941 i.e. row 941, col 75 
#define CA800mGRID_NAME "CA800m"
#define CA800mLATITUDE0 42.09583333 /* 42° 05’ 45” N */
#define CA800mLONGITUDE0 -124.3958333 /* 124° 23’ 45” W */
#define CA800mCELL_SPACING 0.008333333333 /* 30 arc-seconds, expressed as decimal degrees */
#define CA800mNROWS_GRID 1148 /* was 719; 1148 reaches the southern border at 32°32'S */
#define CA800mNCOLS_GRID 1243 /* was 743; 1243 reaches the eastern border at 114°08'W */ 

// YNP800m - Yosemite National Park
// This is a subgrid of CA800m defined by -r 468,552 -c 540,623
// This is also a subgrid of US800m defined by -r 1409,1493 -c 615,698
#define YNP800mGRID_NAME "YNP800m"
#define YNP800mLATITUDE0 38.19583333 /* 38° 11’ 45” N */
#define YNP800mLONGITUDE0 -119.8958333 /* 119° 53’ 45” W */
#define YNP800mCELL_SPACING 0.008333333333 /* 30 arc-seconds, expressed as decimal degrees */
#define YNP800mNROWS_GRID 85
#define YNP800mNCOLS_GRID 84

// Pacific Northwest 800m
/* boundaries of "PNW grid" (Brendan's grid)
north: 49deg 00' 15" N
south: 41deg 54' 15" N
east: 118deg 45' 45" W
west: 124deg 44' 15" W
852 rows X 717 columns = 610,884 cells */
/* The 0,0 cell of the PNW800m grid is the 34,112 cell of the US800m grid, i.e. row 112, col 34. */
#define PNW800mGRID_NAME "PNW800m"
#define PNW800mLATITUDE0 49.00416667
#define PNW800mLONGITUDE0 -124.7375
#define PNW800mCELL_SPACING 0.008333333333 /* 30 arc-seconds, expressed as decimal degrees */
#define PNW800mNROWS_GRID 847 /* originally 852 */
#define PNW800mNCOLS_GRID 717

// Yellowstone 800m
/* boundaries of "Yellowstone grid" (Lauren's grid)
north: 47deg 29' 45" N
south: 40deg 24' 15" N
east: 108deg 18' 15" W
west: 112deg 41' 45" W
851 rows X 527 columns = 448,477 cells */
#define Yellowstone800mGRID_NAME "Yellowstone800m"
#define Yellowstone800mLATITUDE0 47.49583333
#define Yellowstone800mLONGITUDE0 -112.6958333
#define Yellowstone800mCELL_SPACING 0.008333333333 /* 30 arc-seconds, expressed as decimal degrees */
#define Yellowstone800mNROWS_GRID 851
#define Yellowstone800mNCOLS_GRID 527

// California 12km - the entire state of California.
#define CA12kmGRID_NAME "CA12km"
#define CA12kmLATITUDE0 42.0 /* 42° 00’ 00” N */
#define CA12kmLONGITUDE0 -124.5 /* 124° 30' 00” W */
#define CA12kmCELL_SPACING 0.125 /* 7' 30", expressed as decimal degrees */
#define CA12kmNROWS_GRID 76
#define CA12kmNCOLS_GRID 84

// California 10km - entire state
// an Albers projection, used in the California Scenarios 2006 study that Jim Lenihan did
/*
   10-192-138-127:Scenarios2006 daveconklin$ cat Header
north: 2480000
south: 1200000
east: -1620000
west: -2370000
rows: 128
cols: 75
10-192-138-127:Scenarios2006 daveconklin$ cat map_specs.asc
PROJECTION ALBERS
UNITS METERS
PARAMETERS
29 30 00  1st standard parallel
45 30 00  2nd standard parallel
-96 00 00  Central meridian
23 00 00    Latittude of projection's origin
0.0         False easting
0.0         False northing
 */
#define CA10kmAlbersGRID_NAME "CA10kmAlbers"
#define CA10kmAlbersLATITUDE0 2480000 
#define CA10kmAlbersLONGITUDE0 -2370000
#define CA10kmAlbersCELL_SPACING 10000
#define CA10kmAlbersNROWS_GRID 128
#define CA10kmAlbersNCOLS_GRID 75

// US 10km 
/*
north: 3200000
south: 250000
east: 2300000
west: -2400000
rows: 295
cols: 470
PROJECTION ALBERS
UNITS METERS
PARAMETERS
29 30 00  1st standard parallel
45 30 00  2nd standard parallel
-96 00 00  Central meridian
23 00 00    Latitude of projection's origin
0.0         False easting
0.0        / False northing
 */
#define US10kmAlbersGRID_NAME "US10kmAlbers"
#define US10kmAlbersLATITUDE0 3200000
#define US10kmAlbersLONGITUDE0 -2400000
#define US10kmAlbersCELL_SPACING 10000
#define US10kmAlbersNROWS_GRID 295
#define US10kmAlbersNCOLS_GRID 470

// BLM_PSW4kmAlbers 
// See e-mail from Ken Ferschweiler to Dave Conklin "Projection info for BLM projection" 8/24/11 1747 
/*
north: 2 801 063.802 303 31
south:   881 063.802 303 31
east:   -587 735.297 190 01
west: -2 502 735.297 190 01
rows: 480
cols: 480
116,472 active cells in a grid of 230,400 cells
PROJECTION ALBERS
UNITS METERS
PARAMETERS
29 30 00  1st standard parallel
45 30 00  2nd standard parallel
-96 00 00  Central meridian
23 00 00    Latitude of projection's origin
0.0         False easting
0.0         False northing
 */
#define BLM_PSW4kmAlbersGRID_NAME "BLM_PSW4kmAlbers"
#define BLM_PSW4kmAlbersLATITUDE0 2801063.80230331
#define BLM_PSW4kmAlbersLONGITUDE0 -2502735.29719001
#define BLM_PSW4kmAlbersCELL_SPACING 4000
#define BLM_PSW4kmAlbersNROWS_GRID 480
#define BLM_PSW4kmAlbersNCOLS_GRID 480


/* Time period definitions */

#define YEARS_1895_1909_NAME "1895_1909"
#define YEARS_1895_1909_FIRST_YEAR 1895
#define YEARS_1895_1909_LAST_YEAR 1909  

#define YEARS_1895_1993_NAME "1895_1993"
#define YEARS_1895_1993_FIRST_YEAR 1895
#define YEARS_1895_1993_LAST_YEAR 1993  

#define MONTHS_1895_1993_NAME "MONTHS_1895_1993"
#define MONTHS_1895_1993_FIRST_YEAR 1895
#define MONTHS_1895_1993_LAST_YEAR 1993  

#define YEARS_1895_2006_NAME "1895_2006"
#define YEARS_1895_2006_FIRST_YEAR 1895
#define YEARS_1895_2006_LAST_YEAR 2006  

#define MONTHS_1895_2006_NAME "MONTHS_1895_2006"
#define MONTHS_1895_2006_FIRST_YEAR 1895
#define MONTHS_1895_2006_LAST_YEAR 2006  

#define MONTHS_1900_1914_NAME "MONTHS_1900_1914"
#define MONTHS_1900_1914_FIRST_YEAR 1900
#define MONTHS_1900_1914_LAST_YEAR 1914  

#define MONTHS_1950_2099_NAME "MONTHS_1950_2099"
#define MONTHS_1950_2099_FIRST_YEAR 1950
#define MONTHS_1950_2099_LAST_YEAR 2099 

#define MONTHS_1950_2100_NAME "MONTHS_1950_2100"
#define MONTHS_1950_2100_FIRST_YEAR 1950
#define MONTHS_1950_2100_LAST_YEAR 2100  

#define MONTHS_1961_1990_NAME "MONTHS_1961_1990"
#define MONTHS_1961_1990_FIRST_YEAR 1961
#define MONTHS_1961_1990_LAST_YEAR 1990  

#define MONTHS_1900_2000_NAME "MONTHS_1900_2000"
#define MONTHS_1900_2000_FIRST_YEAR 1900
#define MONTHS_1900_2000_LAST_YEAR 2000  

#define YEARS_1895_2007_NAME "1895_2007"
#define YEARS_1895_2007_FIRST_YEAR 1895
#define YEARS_1895_2007_LAST_YEAR 2007  

#define MONTHS_1895_2008_NAME "MONTHS_1895_2008"
#define MONTHS_1895_2008_FIRST_YEAR 1895
#define MONTHS_1895_2008_LAST_YEAR 2008  

#define YEARS_1895_2008_NAME "1895_2008"
#define YEARS_1895_2008_FIRST_YEAR 1895
#define YEARS_1895_2008_LAST_YEAR 2008  

#define YEARS_1895_2009_NAME "1895_2009"
#define YEARS_1895_2009_FIRST_YEAR 1895
#define YEARS_1895_2009_LAST_YEAR 2009  

#define MONTHS_1895_2009_NAME "MONTHS_1895_2009"
#define MONTHS_1895_2009_FIRST_YEAR 1895
#define MONTHS_1895_2009_LAST_YEAR 2009  

#define YEARS_1895_2010_NAME "1895_2010"
#define YEARS_1895_2010_FIRST_YEAR 1895
#define YEARS_1895_2010_LAST_YEAR 2010  

#define MONTHS_1895_2010_NAME "MONTHS_1895_2010"
#define MONTHS_1895_2010_FIRST_YEAR 1895
#define MONTHS_1895_2010_LAST_YEAR 2010  

#define MONTHS_1895_2099_NAME "MONTHS_1895_2099"
#define MONTHS_1895_2099_FIRST_YEAR 1895
#define MONTHS_1895_2099_LAST_YEAR 2099  

#define MONTHS_1895_2100_NAME "MONTHS_1895_2100"
#define MONTHS_1895_2100_FIRST_YEAR 1895
#define MONTHS_1895_2100_LAST_YEAR 2100  

#define YEARS_1895_1896_NAME "1895_1896"
#define YEARS_1895_1896_FIRST_YEAR 1895
#define YEARS_1895_1896_LAST_YEAR 1896  

#define YEARS_1895_1950_NAME "1895_1950"
#define YEARS_1895_1950_FIRST_YEAR 1895
#define YEARS_1895_1950_LAST_YEAR 1950  

#define YEARS_1900_1999_NAME "1900_1999"
#define YEARS_1900_1999_FIRST_YEAR 1900
#define YEARS_1900_1999_LAST_YEAR 1999  

#define MONTHS_1900_1999_NAME "MONTHS_1900_1999"
#define MONTHS_1900_1999_FIRST_YEAR 1900
#define MONTHS_1900_1999_LAST_YEAR 1999  

#define YEARS_1901_2000_NAME "1901_2000"
#define YEARS_1901_2000_FIRST_YEAR 1901
#define YEARS_1901_2000_LAST_YEAR 2000  

#define MONTHS_1901_2000_NAME "MONTHS_1901_2000"
#define MONTHS_1901_2000_FIRST_YEAR 1901
#define MONTHS_1901_2000_LAST_YEAR 2000  

#define YEARS_1961_1990_NAME "1961_1990"
#define YEARS_1961_1990_FIRST_YEAR 1961
#define YEARS_1961_1990_LAST_YEAR 1990 

#define YEARS_1971_2000_NAME "1971_2000"
#define YEARS_1971_2000_FIRST_YEAR 1971
#define YEARS_1971_2000_LAST_YEAR 2000  

#define YEARS_1998_2007_NAME "1998_2007"
#define YEARS_1998_2007_FIRST_YEAR 1998
#define YEARS_1998_2007_LAST_YEAR 2007  

#define YEARS_2000_2099_NAME "2000_2099"
#define YEARS_2000_2099_FIRST_YEAR 2001
#define YEARS_2000_2099_LAST_YEAR 2100  

#define MONTHS_2000_2099_NAME "MONTHS_2000_2099"
#define MONTHS_2000_2099_FIRST_YEAR 2001
#define MONTHS_2000_2099_LAST_YEAR 2100  

#define YEARS_2001_2007_NAME "2001_2007"
#define YEARS_2001_2007_FIRST_YEAR 2001
#define YEARS_2001_2007_LAST_YEAR 2007  

#define YEARS_2001_2100_NAME "2001_2100"
#define YEARS_2001_2100_FIRST_YEAR 2001
#define YEARS_2001_2100_LAST_YEAR 2100  

#define MONTHS_2001_2100_NAME "MONTHS_2001_2100"
#define MONTHS_2001_2100_FIRST_YEAR 2001
#define MONTHS_2001_2100_LAST_YEAR 2100  

#define MONTHS_2006_2100_NAME "MONTHS_2006_2100"
#define MONTHS_2006_2100_FIRST_YEAR 2006
#define MONTHS_2006_2100_LAST_YEAR 2100  

#define YEARS_2007_2099_NAME "2007_2099"
#define YEARS_2007_2099_FIRST_YEAR 2007
#define YEARS_2007_2099_LAST_YEAR 2099  

#define YEAR_2009_2009_NAME "2009_2009"
#define YEAR_2009_2009_FIRST_YEAR 2009
#define YEAR_2009_2009_LAST_YEAR 2009  

#define YEARS_2009_2100_NAME "2009_2100"
#define YEARS_2009_2100_FIRST_YEAR 2009
#define YEARS_2009_2100_LAST_YEAR 2100  

#define YEARS_2011_2099_NAME "2011_2099"
#define YEARS_2011_2099_FIRST_YEAR 2011
#define YEARS_2011_2099_LAST_YEAR 2099  

#define YEARS_2011_2100_NAME "2011_2100"
#define YEARS_2011_2100_FIRST_YEAR 2011
#define YEARS_2011_2100_LAST_YEAR 2100  

#define MONTHS_2007_2099_NAME "MONTHS_2007_2099"
#define MONTHS_2007_2099_FIRST_YEAR 2007
#define MONTHS_2007_2099_LAST_YEAR 2099  

#define MONTHS_12_NAME "12_MONTHS"
#define MONTHS_12_FIRST_MONTH 0
#define MONTHS_12_LAST_MONTH 11

#define TIMELESS_NAME "TIMELESS"
#define TIMELESS_FIRST_YEAR 0
#define TIMELESS_LAST_YEAR 0  

#define TIMELESS_2D_NAME "TIMELESS_2D"
#define TIMELESS_2D_FIRST_YEAR 0
#define TIMELESS_2D_LAST_YEAR 0  


void instructions();


void getGrid(char * grid_name, double * latitude0P, double * longitude0P, double * cell_spacingP,
		int * nrows_gridP, int * ncols_gridP)
{
	if (strcmp(grid_name, US30asGRID_NAME)==0)
	{
		*latitude0P = US30asLATITUDE0;
		*longitude0P = US30asLONGITUDE0;
		*cell_spacingP = US30asCELL_SPACING;
		*nrows_gridP = US30asNROWS_GRID;
		*ncols_gridP = US30asNCOLS_GRID;
	}
	else if (strcmp(grid_name, US800mGRID_NAME)==0)
	{
		*latitude0P = US800mLATITUDE0;
		*longitude0P = US800mLONGITUDE0;
		*cell_spacingP = US800mCELL_SPACING;
		*nrows_gridP = US800mNROWS_GRID;
		*ncols_gridP = US800mNCOLS_GRID;
	}
	else if (strcmp(grid_name, USwest800mGRID_NAME)==0)
	{
		*latitude0P = USwest800mLATITUDE0;
		*longitude0P = USwest800mLONGITUDE0;
		*cell_spacingP = USwest800mCELL_SPACING;
		*nrows_gridP = USwest800mNROWS_GRID;
		*ncols_gridP = USwest800mNCOLS_GRID;
	}

	else if (strcmp(grid_name, US4k_Toon6GRID_NAME)==0)
	{
		*latitude0P = US4k_Toon6LATITUDE0;
		*longitude0P = US4k_Toon6LONGITUDE0;
		*cell_spacingP = US4k_Toon6CELL_SPACING;
		*nrows_gridP = US4k_Toon6NROWS_GRID;
		*ncols_gridP = US4k_Toon6NCOLS_GRID;
	}  

	else if (strcmp(grid_name, US4kmGRID_NAME)==0)
	{
		*latitude0P = US4kmLATITUDE0;
		*longitude0P = US4kmLONGITUDE0;
		*cell_spacingP = US4kmCELL_SPACING;
		*nrows_gridP = US4kmNROWS_GRID;
		*ncols_gridP = US4kmNCOLS_GRID;
	}
	else if (strcmp(grid_name, CA800mGRID_NAME)==0)
	{
		*latitude0P = CA800mLATITUDE0;
		*longitude0P = CA800mLONGITUDE0;
		*cell_spacingP = CA800mCELL_SPACING;
		*nrows_gridP = CA800mNROWS_GRID;
		*ncols_gridP = CA800mNCOLS_GRID;
	}
	else if (strcmp(grid_name, YNP800mGRID_NAME)==0)
	{
		*latitude0P = YNP800mLATITUDE0;
		*longitude0P = YNP800mLONGITUDE0;
		*cell_spacingP = YNP800mCELL_SPACING;
		*nrows_gridP = YNP800mNROWS_GRID;
		*ncols_gridP = YNP800mNCOLS_GRID;
	}
	else if (strcmp(grid_name, PNW800mGRID_NAME)==0)
	{
		*latitude0P = PNW800mLATITUDE0;
		*longitude0P = PNW800mLONGITUDE0;
		*cell_spacingP = PNW800mCELL_SPACING;
		*nrows_gridP = PNW800mNROWS_GRID;
		*ncols_gridP = PNW800mNCOLS_GRID;
	}
	else if (strcmp(grid_name, Yellowstone800mGRID_NAME)==0)
	{
		*latitude0P = Yellowstone800mLATITUDE0;
		*longitude0P = Yellowstone800mLONGITUDE0;
		*cell_spacingP = Yellowstone800mCELL_SPACING;
		*nrows_gridP = Yellowstone800mNROWS_GRID;
		*ncols_gridP = Yellowstone800mNCOLS_GRID;
	}
	else if (strcmp(grid_name, CA12kmGRID_NAME)==0)
	{
		*latitude0P = CA12kmLATITUDE0;
		*longitude0P = CA12kmLONGITUDE0;
		*cell_spacingP = CA12kmCELL_SPACING;
		*nrows_gridP = CA12kmNROWS_GRID;
		*ncols_gridP = CA12kmNCOLS_GRID;
	}  
	else if (strcmp(grid_name, CA10kmAlbersGRID_NAME)==0)
	{
		*latitude0P = CA10kmAlbersLATITUDE0;
		*longitude0P = CA10kmAlbersLONGITUDE0;
		*cell_spacingP = CA10kmAlbersCELL_SPACING;
		*nrows_gridP = CA10kmAlbersNROWS_GRID;
		*ncols_gridP = CA10kmAlbersNCOLS_GRID;
	}  
	else if (strcmp(grid_name, US10kmAlbersGRID_NAME)==0)
	{
		*latitude0P = US10kmAlbersLATITUDE0;
		*longitude0P = US10kmAlbersLONGITUDE0;
		*cell_spacingP = US10kmAlbersCELL_SPACING;
		*nrows_gridP = US10kmAlbersNROWS_GRID;
		*ncols_gridP = US10kmAlbersNCOLS_GRID;
	}  
	else if (strcmp(grid_name, BLM_PSW4kmAlbersGRID_NAME)==0)
	{
		*latitude0P = BLM_PSW4kmAlbersLATITUDE0;
		*longitude0P = BLM_PSW4kmAlbersLONGITUDE0;
		*cell_spacingP = BLM_PSW4kmAlbersCELL_SPACING;
		*nrows_gridP = BLM_PSW4kmAlbersNROWS_GRID;
		*ncols_gridP = BLM_PSW4kmAlbersNCOLS_GRID;
	}  
	else if (strcmp(grid_name, ATtest800mGRID_NAME)==0)
	{
		*latitude0P = ATtest800mLATITUDE0;
		*longitude0P = ATtest800mLONGITUDE0;
		*cell_spacingP = ATtest800mCELL_SPACING;
		*nrows_gridP = ATtest800mNROWS_GRID;
		*ncols_gridP = ATtest800mNCOLS_GRID;
	}  
	else if (strcmp(grid_name, BearSoil800mGRID_NAME)==0)
	{
		*latitude0P = BearSoil800mLATITUDE0;
		*longitude0P = BearSoil800mLONGITUDE0;
		*cell_spacingP = BearSoil800mCELL_SPACING;
		*nrows_gridP = BearSoil800mNROWS_GRID;
		*ncols_gridP = BearSoil800mNCOLS_GRID;
	}  
	else if (strcmp(grid_name, BearSoilToon6GRID_NAME)==0)
	{
		*latitude0P = BearSoilToon6LATITUDE0;
		*longitude0P = BearSoilToon6LONGITUDE0;
		*cell_spacingP = BearSoilToon6CELL_SPACING;
		*nrows_gridP = BearSoilToon6NROWS_GRID;
		*ncols_gridP = BearSoilToon6NCOLS_GRID;
	}  
	else if (strcmp(grid_name, NA8kmGRID_NAME)==0)
	{
		*latitude0P = NA8kmLATITUDE0;
		*longitude0P = NA8kmLONGITUDE0;
		*cell_spacingP = NA8kmCELL_SPACING;
		*nrows_gridP = NA8kmNROWS_GRID;
		*ncols_gridP = NA8kmNCOLS_GRID;
	}  
	else if (strcmp(grid_name, VEMAPGRID_NAME)==0)
	{
		*latitude0P = VEMAPLATITUDE0;
		*longitude0P = VEMAPLONGITUDE0;
		*cell_spacingP = VEMAPCELL_SPACING;
		*nrows_gridP = VEMAPNROWS_GRID;
		*ncols_gridP = VEMAPNCOLS_GRID;
	}  
	else if (strcmp(grid_name, GlobalGRID_NAME)==0)
	{
		*latitude0P = GlobalLATITUDE0;
		*longitude0P = GlobalLONGITUDE0;
		*cell_spacingP = GlobalCELL_SPACING;
		*nrows_gridP = GlobalNROWS_GRID;
		*ncols_gridP = GlobalNCOLS_GRID;
	}  
	else if (strcmp(grid_name, GlobalToon5GRID_NAME)==0)
	{
		*latitude0P = GlobalToon5LATITUDE0;
		*longitude0P = GlobalToon5LONGITUDE0;
		*cell_spacingP = GlobalToon5CELL_SPACING;
		*nrows_gridP = GlobalToon5NROWS_GRID;
		*ncols_gridP = GlobalToon5NCOLS_GRID;
	}  
	else if (strcmp(grid_name, GlobalToon3GRID_NAME)==0)
	{
		*latitude0P = GlobalToon3LATITUDE0;
		*longitude0P = GlobalToon3LONGITUDE0;
		*cell_spacingP = GlobalToon3CELL_SPACING;
		*nrows_gridP = GlobalToon3NROWS_GRID;
		*ncols_gridP = GlobalToon3NCOLS_GRID;
	}  
	else 
	{
		printf("*** Did not recognize grid_name %s.\n", grid_name);
		instructions();
	}
} // end of getGrid()


void getTime(char * time_name, int * first_yearP, int * last_yearP, 
		Boolean * timeless_flagP, Boolean * month_flagP, Boolean * twoD_flagP)  
{ // Get the time frame and initialize the time variables.
	*month_flagP = *twoD_flagP = FALSE;
	if (strcmp(time_name, YEARS_1895_2006_NAME)==0)
	{
		*first_yearP = YEARS_1895_2006_FIRST_YEAR;
		*last_yearP = YEARS_1895_2006_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, MONTHS_1895_2006_NAME)==0)
	{
		*first_yearP = MONTHS_1895_2006_FIRST_YEAR;
		*last_yearP = MONTHS_1895_2006_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, YEARS_1895_1993_NAME)==0)
	{
		*first_yearP = YEARS_1895_1993_FIRST_YEAR;
		*last_yearP = YEARS_1895_1993_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, MONTHS_1895_1993_NAME)==0)
	{
		*first_yearP = MONTHS_1895_1993_FIRST_YEAR;
		*last_yearP = MONTHS_1895_1993_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, MONTHS_1900_1914_NAME)==0)
	{
		*first_yearP = MONTHS_1900_1914_FIRST_YEAR;
		*last_yearP = MONTHS_1900_1914_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, MONTHS_1950_2099_NAME)==0)
	{
		*first_yearP = MONTHS_1950_2099_FIRST_YEAR;
		*last_yearP = MONTHS_1950_2099_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, MONTHS_1950_2100_NAME)==0)
	{
		*first_yearP = MONTHS_1950_2100_FIRST_YEAR;
		*last_yearP = MONTHS_1950_2100_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, MONTHS_1961_1990_NAME)==0)
	{
		*first_yearP = MONTHS_1961_1990_FIRST_YEAR;
		*last_yearP = MONTHS_1961_1990_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, MONTHS_1900_2000_NAME)==0)
	{
		*first_yearP = MONTHS_1900_2000_FIRST_YEAR;
		*last_yearP = MONTHS_1900_2000_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, YEARS_1895_2008_NAME)==0)
	{
		*first_yearP = YEARS_1895_2008_FIRST_YEAR;
		*last_yearP = YEARS_1895_2008_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_1895_2009_NAME)==0)
	{
		*first_yearP = YEARS_1895_2009_FIRST_YEAR;
		*last_yearP = YEARS_1895_2009_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_1895_2010_NAME)==0)
	{
		*first_yearP = YEARS_1895_2010_FIRST_YEAR;
		*last_yearP = YEARS_1895_2010_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, MONTHS_1895_2010_NAME)==0)
	{
		*first_yearP = MONTHS_1895_2010_FIRST_YEAR;
		*last_yearP = MONTHS_1895_2010_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, MONTHS_1895_2008_NAME)==0)
	{
		*first_yearP = MONTHS_1895_2008_FIRST_YEAR;
		*last_yearP = MONTHS_1895_2008_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, MONTHS_1895_2009_NAME)==0)
	{
		*first_yearP = MONTHS_1895_2009_FIRST_YEAR;
		*last_yearP = MONTHS_1895_2009_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, MONTHS_1895_2099_NAME)==0)
	{
		*first_yearP = MONTHS_1895_2099_FIRST_YEAR;
		*last_yearP = MONTHS_1895_2099_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, MONTHS_1895_2100_NAME)==0)
	{
		*first_yearP = MONTHS_1895_2100_FIRST_YEAR;
		*last_yearP = MONTHS_1895_2100_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, YEARS_1895_1909_NAME)==0)
	{
		*first_yearP = YEARS_1895_1909_FIRST_YEAR;
		*last_yearP = YEARS_1895_1909_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_1895_2007_NAME)==0)
	{
		*first_yearP = YEARS_1895_2007_FIRST_YEAR;
		*last_yearP = YEARS_1895_2007_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_1895_1896_NAME)==0)
	{
		*first_yearP = YEARS_1895_1896_FIRST_YEAR;
		*last_yearP = YEARS_1895_1896_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_1895_1950_NAME)==0)
	{
		*first_yearP = YEARS_1895_1950_FIRST_YEAR;
		*last_yearP = YEARS_1895_1950_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_1900_1999_NAME)==0)
	{
		*first_yearP = YEARS_1900_1999_FIRST_YEAR;
		*last_yearP = YEARS_1900_1999_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, MONTHS_1900_1999_NAME)==0)
	{
		*first_yearP = MONTHS_1900_1999_FIRST_YEAR;
		*last_yearP = MONTHS_1900_1999_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, MONTHS_1901_2000_NAME)==0)
	{
		*first_yearP = MONTHS_1901_2000_FIRST_YEAR;
		*last_yearP = MONTHS_1901_2000_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, YEARS_1961_1990_NAME)==0)
	{
		*first_yearP = YEARS_1961_1990_FIRST_YEAR;
		*last_yearP = YEARS_1961_1990_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_2001_2007_NAME)==0)
	{
		*first_yearP = YEARS_2001_2007_FIRST_YEAR;
		*last_yearP = YEARS_2001_2007_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_2000_2099_NAME)==0)
	{
		*first_yearP = YEARS_2000_2099_FIRST_YEAR;
		*last_yearP = YEARS_2000_2099_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, MONTHS_2000_2099_NAME)==0)
	{
		*first_yearP = MONTHS_2000_2099_FIRST_YEAR;
		*last_yearP = MONTHS_2000_2099_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, YEARS_2001_2100_NAME)==0)
	{
		*first_yearP = YEARS_2001_2100_FIRST_YEAR;
		*last_yearP = YEARS_2001_2100_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, MONTHS_2001_2100_NAME)==0)
	{
		*first_yearP = MONTHS_2001_2100_FIRST_YEAR;
		*last_yearP = MONTHS_2001_2100_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, MONTHS_2006_2100_NAME)==0)
	{
		*first_yearP = MONTHS_2006_2100_FIRST_YEAR;
		*last_yearP = MONTHS_2006_2100_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, YEARS_1998_2007_NAME)==0)
	{
		*first_yearP = YEARS_1998_2007_FIRST_YEAR;
		*last_yearP = YEARS_1998_2007_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_1901_2000_NAME)==0)
	{
		*first_yearP = YEARS_1901_2000_FIRST_YEAR;
		*last_yearP = YEARS_1901_2000_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_1971_2000_NAME)==0)
	{
		*first_yearP = YEARS_1971_2000_FIRST_YEAR;
		*last_yearP = YEARS_1971_2000_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_1900_1999_NAME)==0)
	{
		*first_yearP = YEARS_1900_1999_FIRST_YEAR;
		*last_yearP = YEARS_1900_1999_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_2007_2099_NAME)==0)
	{
		*first_yearP = YEARS_2007_2099_FIRST_YEAR;
		*last_yearP = YEARS_2007_2099_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_2011_2099_NAME)==0)
	{
		*first_yearP = YEARS_2011_2099_FIRST_YEAR;
		*last_yearP = YEARS_2011_2099_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_2011_2100_NAME)==0)
	{
		*first_yearP = YEARS_2011_2100_FIRST_YEAR;
		*last_yearP = YEARS_2011_2100_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEAR_2009_2009_NAME)==0)
	{
		*first_yearP = YEAR_2009_2009_FIRST_YEAR;
		*last_yearP = YEAR_2009_2009_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, YEARS_2009_2100_NAME)==0)
	{
		*first_yearP = YEARS_2009_2100_FIRST_YEAR;
		*last_yearP = YEARS_2009_2100_LAST_YEAR;
		*timeless_flagP = FALSE;
	}
	else if (strcmp(time_name, MONTHS_2007_2099_NAME)==0)
	{
		*first_yearP = MONTHS_2007_2099_FIRST_YEAR;
		*last_yearP = MONTHS_2007_2099_LAST_YEAR;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, MONTHS_12_NAME)==0)
	{
		*first_yearP = MONTHS_12_FIRST_MONTH;
		*last_yearP = MONTHS_12_LAST_MONTH;
		*timeless_flagP = FALSE;
		*month_flagP = TRUE;
	}
	else if (strcmp(time_name, TIMELESS_NAME)==0)
	{
		*first_yearP = TIMELESS_FIRST_YEAR;
		*last_yearP = TIMELESS_LAST_YEAR;
		*timeless_flagP = TRUE;
	}
	else if (strcmp(time_name, TIMELESS_2D_NAME)==0)
	{
		*first_yearP = TIMELESS_2D_FIRST_YEAR;
		*last_yearP = TIMELESS_2D_LAST_YEAR;
		*timeless_flagP = TRUE;
		*twoD_flagP = TRUE;
	}
	else 
	{
		printf("*** getTime(): Unrecognizable time_name = %s.\n", time_name);
		instructions();
	}

} // end of getTime()


#endif	/* GRIDSPECS_H */


