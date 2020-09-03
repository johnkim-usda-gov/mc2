#include <iostream>
#include <vector>
#include <map>
#include <iterator>
#include <sys/stat.h>

// #pragma hdrstop

using namespace std;

#define NO_DATA_TOKEN -9999.0f

#ifndef UNIX_MAPSS
#define diagMsg(_dmsg) \
{ \
   CString msg; \
   msg.Format(_dmsg); \
   Report::InfoMsg(msg); \
}
#define GRAIN_SUBDIVISIONS 1 /* number of pieces into which the height and width is divided when calculating cell area fractions */
#define GRAIN_SUBDIVISIONS_STRING "1"
#endif

#define FailMsg(FailMsg) \
	{msg = failureID + (_T(FailMsg)); \
		Report::ErrorMsg(msg); }

#define NC_FILL_FLOAT	(9.9692099683868690e+36f) /* near 15 * 2^119 */ /* from netcdf.h */
#define NC_FILL_INT (-2147483647L)  /* from netcdf.h */

#define INVALID_SOIL_DATA_MSG "*** mapss_model: invalid soil data."

#define EXCESS_H2O				0.05f	/* excess water at maximum lai */
#define INADEQ_H2O				0.025f	/* inadequate water at maximum lai */
# define INCREMENT				0.1f

# define LIFEFORM(_i) \
	((_i == GRASS) ? GRASS : (m_pet_adj == tree_pet_factor) ? TREE : SHRUB)


# define START_MONTH OCT
# define BASE_ZERO 0.0

#define SURFACE		0
#define INTERMEDIATE	1
#define DEEP			2
#define NUMBER_SOIL_LAYERS	(DEEP + 1)
#define NSL NUMBER_SOIL_LAYERS

#define SATURATED		0
#define UNSATURATED	1

/* h2o equilibrium parameters */
# define MAX_DIFF_EQUILIB		0.01f	/* initial pct diff */
# define MAX_DIFF_INCREMENT		0.005f	/* increment pct diff */
# define MAX_CYCLES				100
# define MINIMAL_LAI			0.05f	/* lower bound of lai iteration */
# define MINIMAL_LAI_REDUCTION	-0.01f	/* lower bound of lai reduction step */

# define INIT_GRASS_LAI_MIN		0.0f
# define INIT_WOODY_LAI_MIN		0.1f

#define MONTHS 12

/******
	Sometimes trees and shrubs are different, sometimes they are not.
	When thay are the same, they are called WOODY.
*******/
#define GRASS			0
#define WOODY			1	/*******					********/
#define TREE			1	/******* Note same as WOODY ********/
#define SHRUB			2
#define GRASS_SUM		3
#define TREE_SUM		4
#define SHRUB_SUM		5

# define FROST_ANY_MONTH	13
# define NEVER_FREEZES		14

# define PPT_CUTOFF			6.0

# define ICE			4
# define TUNDRA			5
# define TAIGA_TUNDRA     	6
# define NUMBER_OF_LAI_ZONES	(TAIGA_TUNDRA + 1)

# define EVERGREEN		0
# define DECIDUOUS		1

# define C3GRASS		0
# define C4GRASS		1

# define NEEDLELEAF		0
# define BROADLEAF		1

# define BINARY_MODE	                1
# define CANNOT_EXCEED_EXCESS_H2O	0
# define INCREMENT_MODE CANNOT_EXCEED_EXCESS_H2O

# ifndef MININT
#  define MININT			((int) ~(~(unsigned) 0 >> 1))
# endif /* MININT */

# define UNDIRECTED		0
# define SEEK_FLOOR		1
# define SEEK_CEILING           2
# define SEEK_FINAL		3


typedef struct {
	float			transpire[2],
					evap,
					place_holder, // evap_woody, in MC1
					infiltrate,	/* potential infiltration after evap, etc. */
					soil_h2o[NUMBER_SOIL_LAYERS],
					snow,		/* This is the amount of snow pack. */
					melt,
					surf_run,
          base_run;
	double			swp[2],
					pct_soil_h2o[2];
} State;
  

typedef struct {
	bool	broadleaf,	/* vs. needle, TRUE or FALSE */
      deciduous,	/* vs. evergreen, TRUE or FALSE */
      grassfire,
      mixed;
	MAPSSvegClass		mclass, initial_class;
  VCategory vclass, initial_vclass;
} Rulecat, *RulecatP;


typedef struct {
	float			events,
					melt_rate, // amount of snow which would melt in this month, in snow water equivalent.
					snow;		/* This is monthly snow. */
} Month;


typedef struct {
	float			C3[MONTHS],
					C4[MONTHS],
					C3_C4[MONTHS],
					C4_75[MONTHS],
					C3_75[MONTHS],
					soil_tmp[MONTHS],
					weight[MONTHS];
	float			ratio;
	int				ps_type[MONTHS];
	C3C4Dominance	canopy;
} PsPath;


typedef struct {
	float			beg,
					end;
} Interval;


class MAPSS_model_parameters
{
  public:
    MAPSS_model_parameters();
    MAPSS_model_parameters(MAPSSparameterSetName paramSet, float in_wue);
    ~MAPSS_model_parameters() { }

    void default_parameters();
    void ORNL_model_parameters();
    void US10km_model_parameters(MAPSSparameterSetName paramSet);
    void Wies_Yose_model_parameters();
    void adjust_parameters();
    
    float ice_boundary; // growing degree days base 0 C
    float frost;	// threshold (C) for beginning, end of growing season
    float  n_decid_bound; // northern temp threshold, deciduous
    float  s_decid_bound;  // southern temp threshold, deciduous
    float  evergreen_productivity_tree_lad; 
    float  evergreen_productivity_shrub_lad;
    float  evergreen_gdd; //	growing degree days for evergreen (frost based)
    float  evergreen_events;
    float  evergreen_events_grow_std; //	(1.5 to 2.0)
    float  evergreen_events_grow_min;
    float  evergreen_events_grow_cov;
    float  evergreen_events_grow_rng;	
    float  evergreen_events_grow;
    float  evergreen_events_avg;
    float  evergreen_events_grow_avg; //	(6.0 to 7.5)
    float  evergreen_gdd_ratio;
    float  evergreen_gdd_north;
    float  evergreen_gdd_south;
    
    float forest_threshold; // lai below which forest -> savanna (classif. only)
    float min_tree_lai[NZONES][2];
    float LaiUpperBoundsLifeform;
    float LaiUpperBoundsEnergy[NUMBER_OF_LAI_ZONES][NCL];
    float fire_threshold_1;	/* Threshold for Grass Sum LAI */
    float fire_threshold_2;	/* Threshold for Shrub LAI */
    float fire_threshold_5;	/* Threshold for pet */
    float fire_threshold_6;	/* Threshold for ppt in High Month */
    float max_grass_shrub_threshold;
    float taiga_tundra_boundary[NZONES];
    float tundra_boundary[NZONES];
    float short_grass_threshold;
    float broad_ppt_min; //	minimum monthly ppt for broadleaf, growing season
    float spring_grass; // maximum grass lai in 1st month of growing season
    float grass_snowpack; // max depth of snow expressed as mm of water under which grass will grow 
    float pj_xeric_threshold;
    float pj_max_lai_continental;
    float pj_max_lai_maritime;
    float north_hard_threshold;
    float dry_trop_threshold;
    float maritime_boundary[NZONES];
    float semi_desert_threshold;
    float tall_grass_threshold;
    float desert_grass_sum_threshold;
    float desert_shrub_threshold;
    float desert_grass_threshold;
    float cool_grass_threshold; 
    float max_grass_threshold;
    float xeric_savanna_threshold;
    float mediterranean_savanna_threshold;
    float snow0; //		temp (C) above which snow fraction equals 0 
    float snow1; //		temp (C) below which snow fraction equals 1
    float no_melt; //		temp (C) below which no snow melt occurs
    float melt_slope; //	temp coefficient for snow melt rate (mm)    
    float cond_min_nominal[NZONES][NCL];
    float cond_max_nominal[NZONES][NCL];
    float cond_min[NZONES][NCL][2];
    float cond_max[NZONES][NCL][2];
    float normal_cond_max[NZONES];
    int cond_min_months[NCL];
    float b[NZONES][NCL][2];
    float c[NZONES][NCL][2];    
    float a_slope[NCL];
    const static int k_factor_winter_boundary_upper = 8;
    const static int k_factor_winter_boundary_lower = 4;
    float event_ppt; //	nominal # of rainfall events per unit of rain (i.e. 1/(rain/event))
    float event_pet; //	pet threshold (mm/mo) for max_events determination
    float max_events[2]; //	maximum number of events at pet <= event_pet and maximum number of events at pet > event_pet
    float interc_lai; //	maximum ppt interception per event (mm)
    float full_attenuation_lai; // woody lai threshold for light attenuation
    float no_attenuation_lai; 
    float cond_surface_max[NZONES][NCL];
    float k_transp[NZONES][NCL][2][2];
    float k_trans_addend[NCL][2];
    float k_factor_slope[NZONES][2];
    float k_factor_constraint;
    float k_factor_pet_boundary[NZONES][2];
    float tsg_threshold;
    float wp[NCL];
    float cond_lai_max[NZONES][NCL][2];
    float tree_pet_factor;
    float tallgrass_pet_factor;
    float evergreen_productivity_tree; // tree productivity for growing season
    float evergreen_productivity_shrub; // used for evergreen/decid decision
    int WhichSoils;
    float mapss_subalpine_threshold; // When GDD0>taiga_tundra_boundary and <=subalpine_threshold, EN forests are classified as ForestSubalpine
          // and tree savannas are classified as TreeSavannaSubalpine
    float wue; // water use efficiency, unitless scalar, default value is 1.0
}; // end of class MAPSS_model_parameters


class Soil_MAPSS
{
  public:
    friend class MAPSS_BiogeogModel;
    friend class CENTURY_BiogeochemModel;
    
    Soil_MAPSS(void);
    Soil_MAPSS(SoilDataStruct in_soil_data);
    ~Soil_MAPSS (void);
    
  protected:
    double sand[NUMBER_SOIL_LAYERS],
					clay[NUMBER_SOIL_LAYERS],
					organic[NUMBER_SOIL_LAYERS],
					rock_frag_mineral[NUMBER_SOIL_LAYERS],
					rock_frag_organic[NUMBER_SOIL_LAYERS],
					eff_thickness[NUMBER_SOIL_LAYERS],
					theta_m[NUMBER_SOIL_LAYERS];
    float mineral_depth;
    float organic_depth;
    float	swhc[NUMBER_SOIL_LAYERS];
		float	wilt_pt[NUMBER_SOIL_LAYERS];
    float field_cap[NUMBER_SOIL_LAYERS];
    float awc[NUMBER_SOIL_LAYERS];
    float bulk_density; // bulk density isn't used in the MAPSS model, but it is used in
                        // the CENTURY biogeochemistry model
    
  private:
    bool valid;
    double calc_hycon(double sand, double clay, double h2o, int layer);
    double calc_swp(double pct_soil_h2o, int layer);
    void saturated_drainage(State * site, int layer, float * destination /* h2o in destination layer */ );
    void unsaturated_drainage(State * site, int layer, float * destination /* h2o in destination layer */ );
        
    double a;
    double b;
    double c;
    double d;
    double e;
    double f;
    double g;
    double h;
    double j;
    double k;
    double m;
    double n;
    float rock_frag_max;
    float k_surfrun; //	coeff of surface runoff (increase to reduce runoff)
    
    float matrix_pot; // MPa
    float field_pot;
    
    float thickness[NSL]; // mm
    double pp[NSL];
    double qq[NSL];
    double rr[NSL];
    double tt[NSL];
    double uu[NSL];
    double vv[NSL];
    float KK[2][NSL];
    float exp_perc[2][NSL];

    double aa[NSL];
    double bb[NSL];
    double h2o_sat[NSL]; // amount of water per unit effective thickness at saturation
    float sat2mat[NSL];
    
}; // end of class Soil_MAPSS declaration


class MAPSS_BiogeogModel: public ProcessModel, public MAPSS_model_parameters
   {
   public:
      MAPSS_BiogeogModel() {};
      ~MAPSS_BiogeogModel(void) { }
      
#ifdef UNIX_MAPSS
      MAPSS_BiogeogModel(Simulation * pRun);
      MAPSS_BiogeogModel(MAPSSparameterSetName paramSet, SoilDataStruct soil_data, float * pptP, float * tmpP, 
                         float * vprP, float * wndP, float in_elev, float in_wue);   
      bool runModelAndCollectData(const int year_index, const int row_index, const int col_index);
#endif
   
      void MAPSSmodelInit(MAPSSparameterSetName paramSet, SoilDataStruct soil_data, float * pptP, float * tmpP, float * vprP, 
                          float * wndP, float in_elev);
      MAPSSvegClass mapss_model(MAPSSdataOutputs * in_mapss_outputP);
         // Calls satvp (Sci), CalcTrbxfr (0),  
         // GrowingDegreeDays (Sci), FindRatio (Sci), C3ProdPct (Sci), 
         // BroadDecid (3), initialize_lai (0), LaiCycle (3).
      
      int m_zone; // climate zone
      MAPSSdataOutputs * m_mapss_outputP;
      // float m_c3pct;
      
  private:
    // This first batch are the routines which call other methods of the MAPSS_BiogeogModel class.
    // 16 routines
    void AdjustMaxlai(float excess, State * site, float lai[][2],  
        float grass_lai[], bool deciduous, State curr_mo_State[]); // Calls LightAttenuate (1). 2X
    bool AnnualCycleEquilibrium(State * curr_mo_State, int mo, int * months);
        // Calls CheckStateValues (0). 1X
    void BroadDecid(Rulecat * cats, float productivity, bool IsTree, bool SecondTime);
        // Calls ProcessSeason (2). 2X
    bool CheckEvergreen(State site[], State site_begin[], float offset_begin, float lai_values[], 
        State ** curr_yearP, State ** prev_yearP); 
        // Calls BroadDecid (3), FindMeanLai (0), ClassifyStation (1), 
        // LaiCycle (3). 1X
    MAPSSvegClass ClassifyStation(float lai_values[], C3C4Dominance canopy_type);
        // Calls HeatLimitedRules (0), ForestRules (0), TreeSavannaRules (0), ShrubRules (0), GrassRules (0),
        // DesertRules (0). Called 2X.
    void GrassAlone(State site[], float conductance[][2]);  // Calls LightAttenuate (1),
        // InitSoilWater (0), InitSnow (0), SetCurrPrev (0), DistributePpt (0), transpiration (2). 2X
    void GrassWoodyPart1(State curr_year[], State prev_year[], State site_begin[], float offset_begin);
        // Calls LightAttenuate (1), InitSoilWater (0), InitSnow (0), SetCurrPrev (0), DistributePpt (0),
        // transpiration (2), GrassWoodyEquilibrium (3), FindMeanLai (0), pet_adjust. 2X
    bool GrassWoodyEquilibrium(State site[], int mo, int * months, bool * shrubs, float deficit[], 
        bool deciduous, State curr_mo_State[]);
        // Calls AdjustMaxlai (2), minimum_reserve (0), CheckStateValues (0), NoDeficits (0),
        // LightAttenuate (1), pet_adjust (0). 1X
    void LaiCycle(float * offset, State site_begin[], float offset_begin, State ** curr_yearP, State ** prev_yearP);
        // Calls InitConductance (0), InitSoilWater (0), InitSnow (0), SetCurrPrev (0), DistributePpt (0),
        // transpiration (2), AnnualCycleEquilibrium (1), AdjustLaiIfUnbalanced (1), SwapYears (0), NewMonthlyLai (0), 
        // pet_adjust (0). 2X
    void LightAttenuate(float maxlai, float grass_lai[], float lai[][2], bool deciduous, State site[]); 
        // Calls NewMonthlyLai (0). 5X
    void ProcessSeason(float frost_line, float * ppt_ptr); // Calls TimeInterval (1), MonthMin (0). 1X
    bool AdjustLaiIfUnbalanced(State * * curr_year, State * * prev_year, int * inc_mode, int lifeform); 
        // Calls minimum_reserve (0), SwapYears (0). 1X
    void TimeInterval(float curr_tmp, float next_tmp, float frost_line, Interval * period);
        // Calls OneThreshold (0). 1X
    float transpiration(State * curr_mo_State, float conductance[][2], float lai[], float grass_lai[], 
        int mo); // Calls infiltrate (0), Soil_MAPSS methods, TranspireStep (1). 3X
    float TranspireStep(State * curr_mo_State, float conductance[][2], float lai[], float * grass_lai, int mo);
        // Calls StomatalConductance (0), AccumulateTranspiration (0), GrassPotAtSatisfied (0). 1X

    // The next batch of routines are primitives in the sense that they do not call any other 
    // MAPSS_BiogeogModel methods (although they may call ScienceFcns or Soil_MAPSS methods).
    // 18 routines
    float AccumulateTranspiration(State * curr_mo_State, float * pot_transp, float water, int soil); // Called 4X.
    void CalcTrbxfr(int lifeform); // Called 4X.
    bool CheckStateValues(State * curr_mo_State, State lastyear[], float threshold, int mo, float base); // 3X
    MAPSSvegClass DesertRules(int m_zone); // 3X
    void DistributePpt(State * curr_mo_State, float * lai, int mo); // 3X
    float FindMeanLai(float * total_lai, int lifeform); // 6X
    MAPSSvegClass ForestRules(float lai_values[], int m_zone); // 1X
    bool GrassPotAtSatisfied(float pot_at, State * curr_mo_State, float conductance, float lai[], 
        float * grass_lai, float begin_h2o, float pot_transp); // 1X
    MAPSSvegClass GrassRules(float lai_values[], int m_zone, C3C4Dominance canopy_type); // 3X
    MAPSSvegClass HeatLimitedRules(int m_zone); // 1X
    void infiltrate(State * site); // 3X
    void InitConductance(float conductance[][2], float grass, float woody); // 2X
    float initialize_lai(); // 2X
    void InitSnow(State * curr_year); // 3X
    void InitSoilWater(State * curr_year); // 3X
    float minimum_reserve(State * site); // 3X
    void MonthMin(float current, float next, float prior, float * min, Interval * period); // 1X
    void NewMonthlyLai(float maxlai, float grass_lai[], float grasslai, float lai[][2],
        bool deciduous, State site[]); // 2X
    bool NoDeficits(float deficit[]); // 1X
    float OneThreshold(float curr_tmp, float next_tmp, float frost_line); // 2X
    void pet_adjust(); // 6X 
    void PS_CenturyNew(float lai[][2], PsPath * biomass); // 1X
    int SetCurrPrev(int months, State site[], State * curr_mo_State[], State * prev_mo[], float * * currlai); // 3X
    MAPSSvegClass ShrubRules(float lai_values[], int m_zone); // 3X
    float StomatalConductanceMC2(State * site, float pet, int lifeform, const float conductance[][2], 
        int mo, bool shrubs, bool broadleaf); // 1X
    void SwapYears(State * * curr_year, State * * prev_year); // 2X
    MAPSSvegClass TreeSavannaRules(float lai_values[], int m_zone); // 1X
    
    float m_at_woody_ppt;
    float	m_at_woody_ppt_norm;
    C3C4Dominance m_canopy_type;
    float m_capacity; // awc in top 2 soil layers
    Rulecat m_cats; // deciduous, broadleaf, mclass, ...
    float	m_conductance[MONTHS][2];
    float m_elevation;
    float m_equilibrium_factor;
    float m_events_summer_avg; // average # of rainfall events per month during the growing season
    float m_excess_h2o[NSL];
    float m_grass_lai[MONTHS];
    float m_growing_deg_days_frost;
    float m_growing_deg_days_zero;
    float m_inadeq_h2o[NSL];
    float m_lai[MONTHS][2];
    float	m_lai_values[NCL*2];
    float m_k_factor;
    float m_LaiUpperBoundEnergy[NCL];
    int m_length_winter;
    float m_maxlai;
    float m_max_grass;
    float m_max_tmp; // mean temperature of warmest month 
    int m_max_tmp_mo; // warmest month (Jan=0)
    float m_min_tmp; // mean temperature of coolest month
    float m_mix_ratio; // Gets compared to evergreen_gdd_ratio
    float m_mix_slope; // Used to calculate m_mix_ratio.
    float m_mix_y_intercept; // Used to calculate m_mix_ratio.
    Month m_months[MONTHS]; // monthly storm events, snow, melt_rate
    float	m_new_tree_lai;
    float m_offset;
    float	m_old_tree_lai;
    float m_pet_adj; // pot. ET fudge factor for the current veg (e.g. 1.22 for tallgrass, 2.55 for trees)
		float m_pet_array[NCL][MONTHS];
    float * m_petP; // Contains the address of the first of 12 monthly values for a particular "lifeform"
        // in m_pet_array[NCL][MONTHS]
    float m_ppt[MONTHS];
    float m_ppt_winter; // Used in at_woody calculations, etc.
    Soil_MAPSS m_soil;
    int m_thaw; // month of spring thaw (Jan=1)
    float m_tmp[MONTHS];
    float m_total_saturation[NSL]; // [GRASS] and [WOODY]: amount of water available to shallow-rooted
        // plants and deep-rooted plants when soil is saturated
    float m_unavail; // total water remaining in the top 2 soil layers at wilting point 
    float m_vp_sat[MONTHS];
    float m_vpr[MONTHS];
    float m_woody_alone_at;
    float m_wnd[MONTHS];

}; // end of class MAPSS_BiogeogModel declaration


