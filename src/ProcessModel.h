/*
 *  ProcessModel.h
 */

extern ScienceFcns sciFn;

class Simulation; // forward declaration

// Indices for discretionary output variables in VarDict[]
#define PPT 0
#define TMP PPT+1 
#define VPR TMP+1 
#define FPRD_PPT VPR+1 
#define FPRD_TMP FPRD_PPT+1 
#define GPRD_PPT FPRD_TMP+1 
#define GPRD_TMP GPRD_PPT+1 
#define RH GPRD_TMP+1 
#define N_VOLATIL RH+1 
#define C_HARVEST N_VOLATIL+1 
#define TMAX C_HARVEST+1 
#define TMIN TMAX+1 
#define AET TMIN+1 
#define PET AET+1
#define SFC_RUNOFF PET+1
#define H2O_STREAM_FLOW SFC_RUNOFF+1
#define NPP H2O_STREAM_FLOW+1 
#define TDMEAN NPP+1 
#define PPT_SMOOTHED TDMEAN+1 
#define TMAX_SMOOTHED PPT_SMOOTHED+1 
#define FIRE_MARGIN TMAX_SMOOTHED+1 

#define LAST_MONTHLY_VAR FIRE_MARGIN

#define VTYPE LAST_MONTHLY_VAR+1 
#define PHYSIOGNOMIC_CLASS VTYPE+1 
#define BIOME PHYSIOGNOMIC_CLASS+1 
#define PART_BURN BIOME+1 
#define CONSUMED PART_BURN+1 
#define C_VEG CONSUMED+1 
#define C_FOREST C_VEG+1 
#define C_MAX_LIVE_GRASS_ABOVEGR C_FOREST+1 
#define FIRE_KILLED C_MAX_LIVE_GRASS_ABOVEGR+1 
#define GFRAC FIRE_KILLED+1 
#define NEP GFRAC+1 
#define NBP NEP+1 
#define RSP NBP+1 
#define BIO_CONSUME_CENTURY RSP+1 
#define CONSUMED_LIVE BIO_CONSUME_CENTURY+1 
#define CONSUMED_DEAD CONSUMED_LIVE+1 
#define FFMC_ANN_MAX CONSUMED_DEAD+1 
#define BUI_ANN_MAX FFMC_ANN_MAX+1 
#define NPP_TREE BUI_ANN_MAX+1 
#define NPP_GRASS NPP_TREE+1 
#define MONTH_OF_FIRE NPP_GRASS+1 
#define MINERL_5 MONTH_OF_FIRE+1 
#define C_ECOSYS MINERL_5+1 
#define C_ECOSYS_DEC C_ECOSYS+1
#define C_LIVE_ABOVEGR C_ECOSYS_DEC+1
#define C_LIVE_BELOWGR C_LIVE_ABOVEGR+1
#define C_DEAD_ABOVEGR C_LIVE_BELOWGR+1
#define C_DEAD_BELOWGR C_DEAD_ABOVEGR+1
#define C_MAX_FOREST_LEAF C_DEAD_BELOWGR+1
#define C_FINE_BRANCH C_MAX_FOREST_LEAF+1 
#define C_BOLE C_FINE_BRANCH+1 
#define C_MAX_COARSE_ROOT C_BOLE+1 
#define C_MAX_FINE_ROOT C_MAX_COARSE_ROOT+1 
#define C_MAX_LIVE_GRASS_BELOWGR C_MAX_FINE_ROOT+1 
#define C_LITTER C_MAX_LIVE_GRASS_BELOWGR+1 
#define C_LITTER_METAB C_LITTER+1 
#define C_LITTER_STRUC C_LITTER_METAB+1 
#define C_DEAD_WOOD C_LITTER_STRUC+1 
#define C_MAX_STANDING_DEAD C_DEAD_WOOD+1 
#define C_SOIL_AND_LITTER C_MAX_STANDING_DEAD+1 
#define C_SOM_X_STRUC_METAB C_SOIL_AND_LITTER+1 
#define C_SOM C_SOM_X_STRUC_METAB+1 
#define D1HR C_SOM+1 
#define D10HR D1HR+1 
#define D100HR D10HR+1
#define D1000HR D100HR+1 
#define EM_CO D1000HR+1 
#define EM_CO2 EM_CO+1
#define EM_CH4 EM_CO2+1 
#define EM_NMHC EM_CH4+1 
#define EM_PM EM_NMHC+1 
#define C_GRAIN EM_PM+1 
#define TREE_HT C_GRAIN+1  
#define C_MAX_LIVE_GRASS TREE_HT+1 
#define C3_PCT_PROD C_MAX_LIVE_GRASS+1 
#define TREE_TYPE C3_PCT_PROD+1 
#define MIN_SMOOTHED_TMP TREE_TYPE+1 
#define MC_CLIMATE_ZONE MIN_SMOOTHED_TMP+1 
#define TMP_INDEX MC_CLIMATE_ZONE+1 
#define PPT_INDEX TMP_INDEX+1 
#define NEEDLE_INDEX PPT_INDEX+1 
#define EVERGREEN_INDEX NEEDLE_INDEX+1 
#define PSL EVERGREEN_INDEX+1 
#define SSZ_UB PSL+1 
#define PSFZ_LB SSZ_UB+1 
#define SAFZ_LB PSFZ_LB+1 
#define PKLZ_LB SAFZ_LB+1 
#define MATSL PKLZ_LB+1 
#define MAX_SFC_RUNOFF MATSL+1
#define MAX_GRASS_LAI MAX_SFC_RUNOFF+1
#define MAX_TREE_LAI MAX_GRASS_LAI+1
#define SNOWPACK_DAY91 MAX_TREE_LAI+1 
#define CONTINENTALITY SNOWPACK_DAY91+1 
#define GDD CONTINENTALITY+1 
#define SOIL_TMP_MAX GDD+1 
#define FIRE SOIL_TMP_MAX+1 
#define FIRE_UNSUPPRESSED FIRE+1 
#define FIRE_FLI FIRE_UNSUPPRESSED+1
#define FIRE_ROS FIRE_FLI+1 
#define FIRE_ERC FIRE_ROS+1 
#define GROUSE_HABITAT FIRE_ERC+1 
#define GROUSE_SMRPRE GROUSE_HABITAT+1 
#define GROUSE_AUGMAXT GROUSE_SMRPRE+1 
#define GROUSE_ANNTMP GROUSE_AUGMAXT+1 
#define FIRE_MARGIN_ANN_MIN GROUSE_ANNTMP+1 

#define LAST_YEARLY_VAR FIRE_MARGIN_ANN_MIN

#define LAST_MULTIYR_VAR LAST_YEARLY_VAR

#define ELEV LAST_MULTIYR_VAR+1 
#define SWHC_TOP ELEV+1 
#define SWHC_MID SWHC_TOP+1 
#define SWHC_DEEP SWHC_MID+1 
#define SWHC_ALL SWHC_DEEP+1 
#define SOIL_DEPTH SWHC_ALL+1 
#define LAST_SINGLE_VAR SOIL_DEPTH

#define NUM_OPTIONAL_OUTVARS LAST_SINGLE_VAR+1  /* Number of discretionary output variable values  */

// indices in mapssOutvars[]: "#define MAPSS<output variable name> <index>"
#define MAPSSmclass 0
#define MAPSSvclass MAPSSmclass+1
#define MAPSScanopy MAPSSvclass+1
#define MAPSSzone MAPSScanopy+1
#define NUM_MAPSS_OUTVARS (MAPSSzone+1)

// indices in cenOutvars[]: "#define CEN<output variable name> <index>"
#define CENcgrain 39 
#define CENcrmvst 67
#define CENmetabc_1 104
#define CENpet 156
#define CENsomsc 201
#define CENsomtc 205
#define CENstream_1 241
#define CENstrucc_1 251
#define CENstdedc 262
#define CENvolgma 275
#define CENvolexa 276
#define CENvolpla 277
#define CENcrootc 397
#define CENfbrchc 419
#define CENfrootc 430
#define CENrleavc 445
#define CENrlwodc 452
#define CENwood2c 483
#define CENmx_index 498
#define CENc3c4_index 499
#define CENsurface_runoff 501
#define CENtmp_index 502
#define CENppt_index 503
#define CENsymb 504
#define CENafcacc 505
#define CENbgrema 506
#define CENbgrmai_1 507
#define CENbgrmai_2 508
#define CENbgrmae_1 509
#define CENbgrmae_2 510
#define CENbgrmae_3 511
#define CENasmos_subsoil 512
// spares 513-521
#define CENtmp_prev_0 522 /* the first of 12 */
#define CENppt_prev_0 534 /* the first of 12 */
#define CENinitial_class_mapss 546
#define CENinitial_vclass_mapss 547
#define CENefold_t 548
#define CENfinal_year 550
#define CENmapss_canopy 551
#define CENmapss_mix_index 552
#define CENmapss_tmp_index 553
#define CENmapss_ppt_index 554
#define CENmapss_zone 555
#define CENfrost_index_this_year 556
#define CENfire_last_month CENfrost_index_this_year+1 
#define CENequil_time CENfire_last_month+1 
#define CENvclass_mapss CENequil_time+1 
#define CENclass_mapss CENvclass_mapss+1 

#define LAST_CEN_OUTVAR_INDEX CENclass_mapss
#define NUM_EQ_OUTVARS (LAST_CEN_OUTVAR_INDEX+1) 

// indices in FSoutputList[] - fire model state variables
#define FSsum_ann_ppt LAST_CEN_OUTVAR_INDEX+1 
#define FS_Pprev FSsum_ann_ppt+1
#define FS_Kprev FS_Pprev+1
#define FS_Ksum FS_Kprev+1
#define FSprev_l1hr FS_Ksum+1
#define FSprev_dstand FSprev_l1hr+1
#define FSprev_litter FSprev_dstand+1
#define FSprev_1hr FSprev_litter+1
#define FSprev_10hr FSprev_1hr+1
#define FSprev_100hr FSprev_10hr+1
#define FSprev_1000hr FSprev_100hr+1
#define FSprev_day_mc_100hr FSprev_1000hr+1 
#define FSprev_day_mc_1000hr FSprev_day_mc_100hr+1
#define FSprev_mxd FSprev_day_mc_1000hr+1
#define FSprev_mc_grass FSprev_mxd+1
#define FSprev_mc_tree FSprev_mc_grass+1
#define FSprev_depth FSprev_mc_tree+1
#define FSclai_in_midseason FSprev_depth+1
#define FSprev_snow FSclai_in_midseason+1
#define FS_F1 FSprev_snow+1
#define FS_F2 FS_F1+1
#define FSrand_seed_upper FS_F2+1
#define FSrand_seed_lower FSrand_seed_upper+1
#define FSyrs_in_sum_ann_ppt FSrand_seed_lower+1
#define FSyrs_since_fire FSyrs_in_sum_ann_ppt+1 
#define FSffmc_prev FSyrs_since_fire+1
#define FSdmc_prev FSffmc_prev+1
#define FSdc_prev FSdmc_prev+1 

#define LAST_FS_OUTVAR_INDEX FSdc_prev
#define NUM_FS_OUTVARS (LAST_FS_OUTVAR_INDEX-LAST_CEN_OUTVAR_INDEX)

// indices in BSoutputList - biogeography model state variables
#define BSvtype LAST_FS_OUTVAR_INDEX+1 
#define BSsmoothed_max_tree_lai BSvtype+1 
#define BSsmoothed_max_grass_lai BSsmoothed_max_tree_lai+1 

#define LAST_BS_OUTVAR_INDEX BSsmoothed_max_grass_lai 
#define NUM_BS_OUTVARS (LAST_BS_OUTVAR_INDEX-LAST_FS_OUTVAR_INDEX)

#define CENtmax_prev_0 LAST_BS_OUTVAR_INDEX+1 /* the first of 12 */
#define CENtmin_prev_0 CENtmax_prev_0+12 /* the first of 12 */
#define CENtdmean_prev_0 CENtmin_prev_0+12 /* the first of 12 */

#define LAST_WS_OUTVAR_INDEX CENtdmean_prev_0+11
#define NUM_WS_OUTVARS (LAST_WS_OUTVAR_INDEX+1) 

#define NUMBER_OF_ZONES	(mapssTROPICAL + 1)
#define NZONES NUMBER_OF_ZONES

// #define NLAYER_LE_9_FLAG 0  obsoleted in version 2B107
/* TRUE => limit nlayer to 9 for backward compatibility with some versions prior to 2B67 */
#define REGRESSION_TEST_2B101_FLAG 0 /* TRUE => changes for backward compatibility to versions prior to 2B101 */
#define FAST_WARMSTART_IN_FLAG 1 /* TRUE => read warmstart data as a single array rather than
				    as individual variables */
#define FAST_WARMSTART_OUT_FLAG 2 /* TRUE => write warmstart data as a single array rather than
				     as individual variables */
#define ALT_FUEL_LOAD_CODE_FLAG 3 /* TRUE => use alternate fuel load code */
#define ALT_TREE_ALLOMETRY_CODE_FLAG 4 /* TRUE => use alternate tree allometry code */
#define MAPSS_IN_SPINUP_CODE_FLAG 5 /* TRUE => call MAPSS at the beginning of spinup */
/* #define MASK_GE_997000_CODE_FLAG 6  obsoleted in versions after 2B89 */
/* Interpret mask values >= 997000 as cells to be skipped
   These are USFS PVT values for wetland, barren, water or ice, and not modeled. */
#define REGRESSION_TEST_2B89_FLAG 6 /* TRUE => changes for backward compatibility to versions prior to 2B89 */
#define FIRE_CODE_FLAG 7 /* TRUE => use m_vegc_mo<=105 instead of c_all_abovegr_mo<60 in fire threshold */
#define TREE_N_FLAG 8 /* TRUE => skip over the logic for unlimited N in restrp.F 
			 for backward compatibility with versions prior to 2B65 */
#define PRESCRIBED_NLAYER_FLAG 9 /* TRUE => use NLAYER from vvegTypeNN.100 instead of from 
				    soil depth in soil netCDF file */ 

#define CEN_STATE_SIZE 200
#define CEN_LENGTH 600
#define NUM_CODE_FLAGS 10
#define NCL 3 /* number of canopy layers */

#define NO_FIRE -1


typedef enum {unknownRunMode, MAPSS_EQ, CENTURY_EQ, SPINUP, TRANSIENT} RunModeEnum;				
typedef enum {unknownBaseCalibration = 0, 
	/* mc2VEMAP, mc2NA8km, mc2CA08, mc2YOSE, mc2WIES_YOSE, mc2VINCERA, mc2US10kmAlbers, */
	mc2GLOBAL, 
	mc2ConUS, // Conterminous US
	mc2ConUS_LC, // Conterminous US, Land Carbon project, uses a biogeography in the Lynx-WWETAC line
	mc2W_WA, // western Washington
	mc2California,
	mc2BlueMtns,
	mc2END} BaseCalibrationEnum;
typedef enum {NoSoilsData = '\0', FaoSoilsData, ScsSoilsData} SoilsType;
typedef enum 
{
	unknownGridName,
	Global, 
	VEMAP, 
	US12km, 
	NA8km, 
	US4km, 
	USeast4km, 
	US800m, 
	USwest800m, 
	BearSoil800m, 
	ATtest800m, 
	CA800m, 
	YNP800m, 
	PNW800m, 
	Yellowstone800m, 
	CA12km, 
	CA10kmAlbers, 
	US10kmAlbers, 
	BLM_PSW4kmAlbers	
} GridNameEnum;
enum MAPSSparameterSetName {defaultParams = 0, ORNLparams, US10kmFAOparams, US10kmSCSparams, WiesYoseParams };
typedef enum {
	C3C4Ignore = 0,
	C3Dominance,
	C3C4Mixed,
	C4Dominance
} C3C4Dominance;
typedef enum {UNKNOWNtree_typ=0, EN_TREES, EN_DB_TREES, DB_TREES, DB_EB_TREES, EN_EB_TREES, EB_TREES, 
	//                               1         2            3         4            5            6         
	DN_TREES, DN_EN_TREES} TreeType;
	//    7         8
	typedef enum {UNKNOWNagg_vclass=0, CONIFER_FORESTS, WINTER_DECIDUOUS_FORESTS, MIXED_FORESTS, 
		BROADLEAF_EVERGREENandDROUGHT_DECIDUOUS_FORESTS, SAVANNASandWOODLANDS, GRASSandSHRUBLANDS,  
		DESERTS} AggregatedVegClass;
#define MAX_AGG_VCLASS DESERTS

/*
   typedef enum { MAPSSmodel=1, CENTURYmodel, FIREmodel, BIOGEOGmodel } ProcessModel_ID;
#define NUM_PROCESS_MODELS BIOGEOGmodel
*/

typedef enum {
	UNKNOWNbiome=0, 
	DESERTbiome, 
	SHRUBLANDbiome, 
	GRASSLAND_SAVANNAbiome, 
	WOODLANDbiome, 
	FORESTbiome
} BiomeType;
#define MAX_BIOME FORESTbiome

typedef enum {
	UNKNOWNpclass=0, 
	UNVEGETATEDpclass, 
	SEMIDESERT_SHRUBLANDpclass, // as of B101 "shrubland without a significant grass component"
	SEMIDESERT_GRASSLANDpclass,
	SHRUB_STEPPEpclass, // as of B101 "shrub-steppe, meaning shrubland with grass"; was SHRUBLANDpclass before that
	GRASSLANDpclass, 
	SAVANNApclass,
	EVERG_WOODLANDpclass,
	DECID_WOODLANDpclass,
	MIXED_WOODLANDpclass,
	EVERG_FORESTpclass,
	DECID_FORESTpclass,
	MIXED_FORESTpclass
} PhysiognomicClass;
#define MAX_PCLASS MIXED_FORESTpclass

typedef enum {
	UNKNOWNveg=0,
	COLD_BARRENveg, // 1
	TUNDRAveg, // 2
	TAIGA_TUNDRAveg, // 3
	BOREAL_NEEDLELEAF_FORESTveg, // 4
	BOREAL_WOODLANDveg, // 5
	SUBALPINE_FORESTveg, // 6
	MARITIME_EN_FORESTveg, // 7
	TEMPERATE_NEEDLELEAF_FORESTveg, // 8
	TEMPERATE_DB_FORESTveg, // 9
	COOL_MIXED_FORESTveg, // 10 birch-maple
	TEMPERATE_WARM_MIXED_FORESTveg, // 11 oak-hickory
	TEMPERATE_EN_WOODLANDveg, // 12
	TEMPERATE_DB_WOODLANDveg, // 13
	TEMPERATE_COOL_MIXED_WOODLANDveg, // 14
	TEMPERATE_WARM_MIXED_WOODLANDveg, // 15
	SHRUB_STEPPEveg, // 16 prior to 2B101, was C3SHRUBveg
	C3GRASSveg, // 17
	TEMPERATE_DESERTveg, // 18
	SUBTROPICAL_EN_FORESTveg, // 19
	SUBTROPICAL_DB_FORESTveg, // 20
	WARM_EB_FORESTveg, // 21
	SUBTROPICAL_MIXED_FORESTveg, // 22
	SUBTROPICAL_EN_WOODLANDveg, // 23 
	SUBTROPICAL_DB_WOODLANDveg, // 24 
	SUBTROPICAL_EB_WOODLANDveg, // 25
	SUBTROPICAL_MIXED_WOODLANDveg, // 26
	DRY_SHRUB_STEPPEveg, // 27 prior to 2B101, was C4SHRUBveg
	C4GRASSveg, // 28
	SUBTROPICAL_DESERTveg, // 29
	TROPICAL_EB_FORESTveg, // 30
	TROPICAL_DECIDUOUS_WOODLANDveg, // 31
	TROPICAL_SAVANNAveg, // 32
	TROPICAL_SHRUBLANDveg, // 33 spare
	VTYPE34veg, // 34 spare
	TROPICAL_DESERTveg, // 35
	MOIST_TEMPERATE_NEEDLELEAF_FORESTveg, // 36
	VTYPE37veg, // 37 spare
	SUBALPINE_MEADOWveg, // 38
	WATERveg, // 39
	NATURAL_BARRENveg, // 40
	DEVELOPEDveg, // 41
	LARCH_FORESTveg, // 42
	SSZveg, // 43 Sitka spruce zone
	WHZveg, // 44 western hemlock zone
	PSFZveg, // 45 Pacific silver fir zone
	MHZveg, // 46 mountain hemlock zone
	SAFZveg, // 47 subalpine fir zone
	PKLZveg, // 48 subalpine parkland zone
	DRY_TEMPERATE_NEEDLELEAF_FORESTveg, // 49 cool dry needleleaf forest
	BOREAL_SHRUBLANDveg, // 50 boreal shrubland 
	SEMIDESERT_SHRUBLANDveg, // 51 temperate and subtropical shrubland without a significant grass component 
	LPPZveg, // 52 Lodgepole pine zone
	JPZveg, // 53 Jeffrey pine zone
	WWPZveg, // 54 Western white pine zone
	DFZ2veg, // 55 Douglas-fir zone 2
	POCZveg, // 56 Port Orford-cedar zone
	GFZveg, // 57 Grand fir zone
	WFZveg, // 58 White fir zone
	SRFZveg, // 59 Shasta red fir zone
	PPZveg, // 60 Ponderosa pine zone
} MC2VegType;
#define MAX_VTYPE PPZveg

// The #defines below represent equivalence in coded value, not semantic equivalence.
// Veg classes match up pretty well for codes 0 thru 36, but not for codes >36.
#define Lynx_Undefined UNKNOWNveg // 0
#define Lynx_Barren COLD_BARRENveg // 1
#define Lynx_Tundra TUNDRAveg // 2
#define Lynx_Taiga_Tundra TAIGA_TUNDRAveg // 3
#define Lynx_Boreal_Evergreen_Needleleaf_Forest BOREAL_NEEDLELEAF_FORESTveg // 4
#define Lynx_Boreal_Mixed_Woodland BOREAL_WOODLANDveg // 5
#define Lynx_Subalpine SUBALPINE_FORESTveg // 6
#define Lynx_Maritime_Evergreen_Needleleaf_Forest MARITIME_EN_FORESTveg // 7
#define Lynx_Temperate_Evergreen_Needleleaf_Forest TEMPERATE_NEEDLELEAF_FORESTveg // 8
#define Lynx_Temperate_Deciduous_Broadleaf_Forest TEMPERATE_DB_FORESTveg // 9
#define Lynx_Temperate_Cool_Mixed_Forest COOL_MIXED_FORESTveg // 10
#define Lynx_Temperate_Warm_Mixed_Forest TEMPERATE_WARM_MIXED_FORESTveg // 11
#define Lynx_Temperate_Evergreen_Needleleaf_Woodland TEMPERATE_EN_WOODLANDveg // 12
#define Lynx_Temperate_Deciduous_Broadleaf_Woodland TEMPERATE_DB_WOODLANDveg // 13 
#define Lynx_Temperate_Cool_Mixed_Woodland TEMPERATE_COOL_MIXED_WOODLANDveg // 14
#define Lynx_Temperate_Warm_Mixed_Woodland TEMPERATE_WARM_MIXED_WOODLANDveg // 15
#define Lynx_Temperate_Shrubland SHRUB_STEPPEveg // 16
#define Lynx_Temperate_Grassland C3GRASSveg // 17
#define Lynx_Temperate_Desert TEMPERATE_DESERTveg // 18
#define Lynx_Subtropical_Evergreen_Needleleaf_Forest SUBTROPICAL_EN_FORESTveg // 19
#define Lynx_Subtropical_Deciduous_Broadleaf_Forest SUBTROPICAL_DB_FORESTveg // 20
#define Lynx_Subtropical_Evergreen_Broadleaf_Forest WARM_EB_FORESTveg // 21
#define Lynx_Subtropical_Mixed_Forest SUBTROPICAL_MIXED_FORESTveg // 22
#define Lynx_Subtropical_Evergreen_Needleleaf_Woodland SUBTROPICAL_EN_WOODLANDveg // 23
#define Lynx_Subtropical_Deciduous_Broadleaf_Woodland SUBTROPICAL_DB_WOODLANDveg // 24
#define Lynx_Subtropical_Evergreen_Broadleaf_Woodland SUBTROPICAL_EB_WOODLANDveg // 25
#define Lynx_Subtropical_Mixed_Woodland SUBTROPICAL_MIXED_WOODLANDveg // 26
#define Lynx_Subtropical_Shrubland DRY_SHRUB_STEPPEveg // 27
#define Lynx_Subtropical_Grassland C4GRASSveg // 28
#define Lynx_Subtropical_Desert SUBTROPICAL_DESERTveg // 29
#define Lynx_Tropical_Evergreen_Broadleaf_Forest TROPICAL_EB_FORESTveg // 30
#define Lynx_Tropical_Deciduous_Woodland TROPICAL_DECIDUOUS_WOODLANDveg // 31
#define Lynx_Tropical_Savanna TROPICAL_SAVANNAveg // 32
#define Lynx_Tropical_Shrubland TROPICAL_SHRUBLANDveg // 33
#define Lynx_Tropical_Grassland VTYPE34veg // 34
#define Lynx_Tropical_Desert TROPICAL_DESERTveg // 35
#define Lynx_Cool_Needleleaf_Forest MOIST_TEMPERATE_NEEDLELEAF_FORESTveg // 36
#define Lynx_AgricultureGrazing VTYPE37veg // 37
#define Lynx_Developed SUBALPINE_MEADOWveg // 38
#define Lynx_Mining WATERveg // 39

typedef enum {
	LULC_Undefined = 0,			// use normal vtype
	LULC_Default,			// use normal vtype
	LULC_Mechanical_Disturb,	// Execute mechanical disturbance code
	LULC_Agriculture,		// Treat as grass, execute harvest code
	LULC_Developed,			// Treat as barren
	LULC_Mining				// Treat as barren
} LULCType;

typedef enum {
	UNKNOWNpnv=0,
	UNNAMED1pnv,
	UNNAMED2pnv,
	SDZpnv, // 3 Salt desert zone
	GRZpnv, // 4,GRZ,Grassland Zone
	STZpnv, // 5,STZ,Steppe Zone
	WJZpnv, // 6,WJZ,Western Juniper Zone
	UNNAMED7pnv, // 7
	LPPZpnv, // 8,LPPZ,Lodgepole Pine Zone
	SSZpnv, // 9,SSZ,Sitka Spruce Zone
	PPZpnv, // 10,PPZ,Ponderosa Pine Zone
	OWOZpnv, // 11,OWOZ,Oregon White Oak Zone
	JPZpnv, // 12,JPZ,Jeffrey Pine Zone
	POCZpnv, // 13,POCZ,Port Orford-Cedar Zone
	DFZpnv, // 14,DFZ,Douglas-fir Zone
	TOZpnv, // 15,TOZ,Tan Oak Zone
	GFZpnv, // 16,GFZ,Grand Fir Zone
	WWPZpnv, // 17,WWPZ,Western White Pine Zone
	DFZ2pnv, // 18,DFZ2,Douglas-fir 2 Zone
	WHZpnv, // 19,WHZ,Western Hemlock Zone
	WFZpnv, // 20,WFZ,White Fir Zone
	SRFZpnv, // 21,SRFZ,Shasta Red Fir Zone
	PSFZpnv, // 22,SFZ,Pacific Silver Fir Zone
	MHZpnv, // 23,MHZ,Mountain Hemlock Zone
	UNNAMED24pnv, // 24
	SAFZpnv, // 25,SAFZ,Subalpine Fir Zone
	UNNAMED26pnv, // 26
	UNNAMED27pnv, // 27
	UNNAMED28pnv, // 28
	UNNAMED29pnv, // 29
	UNNAMED30pnv, // 30
	UNNAMED31pnv, // 31
	PKLZpnv, // 32,PKLZ,Parkland Zone
	ALPZpnv // 33,ALPZ,Alpine Zone  
} PNVcode;


typedef struct {
	float soils_record[10];
	float bd;
} SoilDataStruct;


typedef struct 
{
	float   sum_ann_ppt;
	float		Pprev;
	float		Kprev;
	float		Ksum;

	// amount of dead fuel in each class in previous month
	float	prev_l1hr;
	float	prev_dstand;
	float	prev_litter;
	float	prev_1hr; 
	float	prev_10hr;
	float prev_100hr;
	float prev_1000hr;

	// moisture content of large dead fuel classes on previous day
	float prev_day_mc_100hr;
	float prev_day_mc_1000hr;

	float		prev_mxd;
	float		prev_mc_grass;
	float		prev_mc_tree;
	float		prev_depth;
	float		clai_in_midseason;
	float		prev_snow;	
	int		F1;
	int 	F2;
	unsigned long int		rand_seed;
	int		yrs_in_sum_ann_ppt;
	float ffmc_prev;
	float dmc_prev;
	float dc_prev;
	int yrs_since_fire;
} FireState;


typedef struct {
	float lai[MONTHS][NCL];
	float conductance[MONTHS][NCL];
	float lai_values[2*NCL];
	float at_woody_ppt;
	float at_woody_ppt_norm;
	float old_tree_lai;
	float new_tree_lai;
	float woody_alone_at;
	float k_factor;
	C3C4Dominance canopy_type;
	float h2o_capacity;
	int mclass;
	const char * msgP;
	int zone;
	float c3c4_ratio;
	float mapssOutvars[NUM_MAPSS_OUTVARS];
} MAPSSdataOutputs;


class RunParamsClass
{  
	public:
		RunParamsClass();
		~RunParamsClass();
		RunModeEnum interpretRunMode(char * run_mode);
		GridNameEnum interpretGridName(char * grid_name);

		// parameters as read from command file
		int first_calendar_year_of_run;
		int years_to_run;
		int years_offset_into_input_data;
		int col_begin;
		int col_end;
		int row_begin;
		int row_end;
		char * grid_name;
		char * run_mode;
		int multiyr_start;
		int multiyr_len;
		char * fire_model_switch;
		char * fire_suppression_switch;
		int fire_suppression_first_year;
		float suppressed_fire_cell_fraction;
		float fire_suppression_fli_threshold;
		float fire_suppression_ros_threshold;
		float fire_suppression_erc_threshold;
		int fire_set_jday;
		int fire_set_interval;
		char * climate_data_directory;
		char * earth_data_directory;
		char * soil_data_file;
		char * soil_bulk_density_file;
		char * CO2_file;
		char * mask_file;
		char * warmstart_file;
		char * output_variable_file;
		char * output_file_prefix;
		// char * model_parameter_set;
		char * monthly_output_switch;
		char * yearly_output_switch;
		char * multiyr_output_switch;
		char * time_invariant_output_switch;
		char * warmstart_output_switch;
		char * dummy_wind_switch;

		// parameters after interpretation
		RunModeEnum runMode;
		GridNameEnum gridName;
		double latitude0, longitude0, cell_spacing;
		bool latlonFlag; // true => coordinates are in degrees N, degrees E
		bool fireModelFlag;
		bool fireSuppressionFlag;
		BaseCalibrationEnum baseCalibration;
		bool diags;
		int years_of_climate_data;
		bool spaceBeforeTimeFlag;
		bool singleCellFlag;
		int num_multiyr_intervals;
		int actual_years_to_run;
		int last_calendar_year_of_run;
		int nominal_multiyr_start;
		bool monthlyOutputFlag; int monthlyOutputSwitchCount;
		bool yearlyOutputFlag; int yearlyOutputSwitchCount;
		bool multiyrOutputFlag; int multiyrOutputSwitchCount;
		bool timeInvariantOutputFlag; int timeInvariantOutputSwitchCount;
		bool warmstartOutputFlag; int warmstartOutputSwitchCount;
		bool dummyWindFlag; // true => use dummy wind instead of reading wnd.nc file

}; // end of class RunParamsClass


class ModelParamsClass
{
	public:
		ModelParamsClass();
		~ModelParamsClass();
		MAPSSparameterSetName interpretMAPSSparameterSet(char * name);

		// parameters as read from command file
		char * MAPSS_parameter_set;
		char * century_path;
		float m_century_runoff_x_intercept;
		float m_century_runoff_slope;
		float m_lait_lower_limit;
		float m_forest_thres_C; /* LOWER LIMIT OF TOTAL WOODY C FOR FOREST, g C m-2 */
		float m_woodl_thres_C;   /* LOWER LIMIT OF TOTAL WOODY C FOR WOODLAND, g C m-2 */
		float tmmin_threshold;
		float maritime_threshold;
		float subalpine_threshold;
		float moist_temperate_threshold; // mmH2O (25" = 635 mm)
		float dry_temperate_threshold; // mmH2O (17" = 432 mm)
		float MAPSS_subalpine_threshold;
		float p_hi_mult; 

		// The next several model parameters affect the fire model.
		float fire_min; // minimum fire effect fraction, i.e. min crown kill fraction in a crown fire
		float fire_ros_min; // minimum rate-of-spread to simulate fire, ft/min
		bool part_burnFlag; // enable or disable part_burn mechanism
		// When part_burnFlag is false, interpret vtype_fire_return_interval_range minimum
		// value as the minimum number of years between fires, and disregard the maximum value.      
		char * alt_fuel_load_switch;

		char * alt_tree_allometry_switch;
		char * unlimited_N_switch;
		float efold_t_max;
		char * southern_hemisphere_switch;
		float desert_grass_C_max;
		float c3_threshold;
		float grassfrac_thres;
		float crown_fire_mortality_pct; // % mortality from crown fire
		float bz_thres; // boreal zone temperature threshold deg C
		float double_CO2_grass_AT_mult; // multiplier that represents the effect of doubled CO2 on transpiration rate of herbaceous plants
		float double_CO2_grass_NPP_mult; // multiplier that represents the effect of doubled CO2 on NPP of herbaceious plants
		float double_CO2_tree_AT_mult; // multiplier that represents the effect of doubled CO2 on transpiration rate of woody plants
		float double_CO2_tree_NPP_mult; // multiplier that represents the effect of doubled CO2 on NPP of woody plants
		float max_LAI[5]; // max LAI of SUPRT	-DN-	-EN-	-DB-	-EB-, projected leaf area m2 m-2
		float max_grass_NPP[3]; // max monthly NPP of SUPRG	-C3-	-C4-, g C m-2 month-1
		float max_tree_NPP[5]; // max monthly NPP of SUPRT	-DN-	-EN-	-DB-	-EB-, g C m-2 month-1
		float needl_hi_threshold; // high threshold for needle index in ProcessModel::treeType()
		float ppt_hi_threshold; // high threshold for precipitation in ProcessModel::treeType()
		float shrub_steppe_precip_threshold; // divides moist shrub-steppe from dry shrub-steppe
		float c4grass_min_summer_precip; // mm H2O total for JJA (N hemi) or 92/90 DJF (S hemi)

		// Next 2 are parameters in the call to pprdwc().  Comments are from pprdwc.F
		// wc is the soil water holding capacity in the top soil layer, normalized by the depth of the layer
		// [0] A: The maximum ratio of available water to pet which would completely limit production assuming wc=0.
		// [1] B:  The effect of wc on the intercept, allows the user to increase the value of the intercept and thereby increase the slope of the line.
		// [2] C:  The lowest ratio of available water to PET at which there is no restriction on production.
		float pprdwc_grass[3]; 
		float pprdwc_tree[3]; 

		// Next 4 are coefficient sets for vegetation zone calculations in the W_WA base calibration.
		// Coefficient subscripts 0-8 correspond to A,B,...,H and P in Table 3, p. 25 of PNW-GTR-841
		float SSZ_ub_coeffs[9]; // Sitka spruce zone upper bound
		float PSFZ_lb_coeffs[9]; // Pacific silver fir zone lower bound
		float SAFZ_lb_coeffs[9]; // subalpine fir zone lower bound
		float PKLZ_lb_coeffs[9]; // subalpine parkland zone lower bound

		// Next 3 are for the grouse habitat model.
		float grouse_anntmp_threshold;
		float grouse_augmaxt_threshold;
		float grouse_smrpre_threshold;

		int regrowth_years; // number of years to wait after a fire before calling the biogeography method

		// parameters after interpretation
		MAPSSparameterSetName MAPSSparameterSet;
		bool altFuelLoadFlag;
		bool altTreeAllometryFlag;
		bool unlimitedNflag;
		bool code_flags[NUM_CODE_FLAGS]; 
		bool southernHemisphereFlag;   
		float z0[NZONES][NCL]; // roughness length 
		float vveg2mfri[MAX_VTYPE+1][3]; // vtype, min_fri, max_fri
		float ffmc_threshold[6][9]; // (1+#climatezones)x(1+#treetypes)
		float bui_threshold[6][9]; // (1+#climatezones)x(1+#treetypes)
		float ffmc_threshold_by_vtype[MAX_VTYPE+1];
		float bui_threshold_by_vtype[MAX_VTYPE+1];

}; // end of class ModelParamsClass


class InputDataClass
{
	public:
		InputDataClass();
		~InputDataClass(); 

		bool maskFlag;
		float elev; // m ASL

		// next 6 are for W_WA calibration
		float fog; // fog factor
		float topomoist; // topographic moisture factor
		float deltaTsl; // difference between mean air temperature for 1950-2000 and mean air temperature at sea level, deg F
		float aspect; // aspect
		float sw; // shortwave
		float cad; // cold air drainage

		double lat, lon;
		double northing, easting;
		SoilDataStruct soilData;
		int num_months;
		float * pptP;
		float * tminP;
		float * tmaxP;
		float * tmpP;
		float * vprP;
		float * tdmeanP;
		float * wndP;
		float * inFromMAPSS;
		float * inFromEQ;
		float * inFromWS;

}; // end of class InputDataClass

/*
// indices for floatOutvars[NUM_FLOAT_OUTVARS]  
#define OUTddecid 0
#define NUM_FLOAT_OUTVARS OUTddecid+1

// indices for intOutvars[NUM_INT_OUTVARS]
#define OUTvtype 0
#define OUTagg_vtype OUTvtype+1
#define OUTagg_vclass OUTagg_vtype // OUTagg_vtype and OUTagg_vclass are synonyms. 
#define OUTzone OUTagg_vclass+1
#define OUTtree_typ OUTzone+1
#define OUTgrass_typ OUTtree_typ+1
#define NUM_INT_OUTVARS OUTgrass_typ+1

class OutputDataClass
{
public:
OutputDataClass() { };
~OutputDataClass() { };

bool validFlag;
float floatOutvars[NUM_FLOAT_OUTVARS];
int intOutvars[NUM_INT_OUTVARS];
}; // end of class OutputDataClass
*/

class MC2centuryEQdataOutputs
{
	public:
		MC2centuryEQdataOutputs() {};
		~MC2centuryEQdataOutputs() {};

	public:
		float cenOutvars[NUM_WS_OUTVARS];

}; // end of class MC2centuryEQdataOutputs


/*
   class BiogeogState
   {
   public:
   bool deciduous;
   bool broadleaf;
   } // end of class BiogeogState
   */

class StateVariables
{
	public:
		void copyFStoCENoutvars();

	public:
		bool deciduous;
		bool broadleaf;
		MC2centuryEQdataOutputs eqState;
		FireState fireState;
		// BiogeogState biogeogState;

}; // end of StateVariables


class ProcessModel
{
	public:
		ProcessModel() { }
		virtual ~ProcessModel() = 0; // " = 0" makes it pure virtual, which makes ProcessModel an abstract class.

#ifdef UNIX_MAPSS
		virtual bool runModelAndCollectData(const int year_index, const int row_index, const int col_index) = 0;
		void ProcessModelInit(Simulation * pRun); // Called by each process model to
		// set up pointers to run parameters and model parameters, etc.
#else
		MAPSSvegClass simulateOnePointMAPSS();   
#endif

		float calcPET(float tmp, float vpr, float vp_sat, float wnd, float elevation, float z, float z0, int days_in_month);   
		VCategory convertMAPSStoVEMAP(MAPSSvegClass mclass); // Called in the process models
		// to convert MAPSS veg classes to VEMAP veg classes
		int convertVEMAPtoVtype(VCategory vclass, BaseCalibrationEnum baseCalib);
		void treeType(float everg, float needl, float ppt_smoothed_yr, BaseCalibrationEnum baseCalib, 
				TreeType * tree_typP, bool * deciduousP, bool * broadleafP);

		// machine-independent pseudo-random number generator
		int portable_rand(void);
		void portable_srand(unsigned long int seed);
		unsigned long int m_portable_next;

	public:
		Simulation * pS;
		ModelParamsClass * modelParamsP;
		RunParamsClass * runParamsP;
		InputDataClass * inputDataP;
		// OutputDataClass outputData;
		MAPSSdataOutputs mapssOutput;
		MC2centuryEQdataOutputs EQdata;
		StateVariables stateData;
		BiomeType m_biome;
		PhysiognomicClass m_physiognomic_class;

	public: // constants
		float c_frost;
		float c_n_decid_bound;
		float c_s_decid_bound;
		float c_az_thres; // upper GDD limit for arctic aka alpine zone, growing degree-days ref to 0 deg C
		// float c_bz_thres; // upper min_smoothed_tmp limit for boreal zone, deg C
		float c_tz_thres; // upper min_smoothed_tmp limit for temperate zone, deg C
		float c_stz_thres; // upper min_smoothed tmp limit for subtropical zone, deg C
		float c_z_all; // z for the trbxfr calculations for all zones and lifeforms

}; // end of class ProcessModel


/* Declarations for the interface with the CENTURY FORTRAN library */
#define logical int /* FORTRAN "logical" is represented by an "int" in C++ */
extern "C" { void     cen_init_climate_(float * ppt, float * tmax, float * tmin); }
extern "C" { void     cen_init_lat_(float * lat, float * elev); }
extern "C" { void     cen_end_(); }
extern "C" { void cen_init_(int * yr, VCategory * vveg, int * diag_flag, int * os_flag, float * state, float * outvars, char * path, 
		float * mx, float * c3c4, int * initflag, int * fm, int * zn, float * tmpi, float * ppti, int * unlimitedNflag, 
		float * fi, char * co2ramp_path); }
extern "C" { void cen_init_soils_(float * bd, float * snd, float * cly, float * depth, float * rock, SoilsType * ws, float * ndep1, float * ndep2); }
extern "C" { void	cen_step_(
		VCategory * prev_vclass, 
		int * wbegin,
		int * wend,
		int * cbegin,
		int * cend,
		int * csene,
		char * wfire,
		char * cfire,			
		float * treec,
		float * cropc,
		float * cen_state,
		float * cen_outvars,
		float * mx,
		float * mx_out,
		float * c3c4,
		float * c3c4_out,
		char * tname,
		char * cname,
		logical * diag_flag,
		logical * os_flag,
		logical * burn_out,
		float * burn_count,
		float * tmpi,
		float * tmpi_out,
		float * ppti,
		float * ppti_out,
		int * zn,
		logical * warmstart_init_flag,
		float * time_in,
		int * rptyrs_in, 
		logical * unlimited_N_flag,
		float * frost_index); }
		extern "C" { void cen_zero_(); }
		extern "C" { int getCenturyParameter(CenParamStruct paramArray[], const char * paramNameIn, int valNdx, float * paramValP); }
		extern "C" { int read_co2ramp_C(char * co2ramp_name, int first_calendar_year, int min_years, int years_offset_into_input_data); }
		extern "C" { int setCenturyParameter(CenParamStruct paramArray[], const char * paramNameIn, int valNdx, float paramVal); }
		extern "C" { void set_constant_CO2(float CO2conc); }
		extern "C" { void stand_step_(float * state, float * outvars, int * eq_flag, float * final_growth, int * diag_flag, int * os_flag, float * fi); }


