/*
 *  MCfire.h
 *  mc2
 */

/* 
#ifndef Boolean
#define Boolean int
#endif
*/ 

#define NUM_EM_SPECIES 7 
typedef enum {PMem=0, PM10em, PM2P5em, COem, CO2em, CH4em, NMHCem} EmissionSpeciesEnum; 

#define NUM_EM_SOURCES 2 
typedef enum {FLAMING=0, SMOLDERING} EmissionSourceEnum;

#define NUM_EM_FUEL_TYPES 8 
typedef enum {DOUGLAS_FIR_SLASH=0, HARDWOODS_SLASH, PONDEROSA_LODGEPOLE_PINE_SLASH,
	MIXED_CONIFER_SLASH, JUNIPER_SLASH, SAGEBRUSH, CHAPARRAL, MEAN_VALUES} EmissionFuelTypeEnum;


typedef struct 
{
	float sg1, sg10, sg100, s1000, sgherb, sgwood;
	float hl, hd;
	float mxd;
	float wndfac;
} FuelCharacteristics;

typedef struct
{
	float snw0, snw1, no_melt, melt_b; // parameters for daily snow model
	float	mc_tree_min, mc_tree_max, mc_grass_min, mc_grass_max; // parameters for live fuel moisture content
	float	slp; // parameter for slope of terrain
	float prob_thres; // probability threshold for complete mortality
} FireParam;

typedef struct 
{
	float mc_1hr[365], mc_10hr[365], mc_100hr[365], mc_1000hr[365];
} FuelDat;


class MC_FireModel: public ProcessModel
{
	public:
		// Use this constructor when instantiating a fire model object with no previous state.
		MC_FireModel(Simulation * pRun, FireState * fireStateP, InputDataClass * inputP);

		// Use this constructor when instantiating a fire model object and restoring a previous state.
		MC_FireModel(Simulation * pRun, float cenOutvars[], InputDataClass * inputP);

		~MC_FireModel() {}

		bool runModelAndCollectData(const int year_index, const int row_index, const int col_index);

		void FireData(float * pptP, float * tmpP, float * wndP, float * vprP, float * petP, 
				float * tminP, float * tmaxP);
		void daily_ppt();
		void daily_dat(float monthly_vals[], float daily_vals[]);
		int daily_values(float monthly_vals[], int mo, float t0, float no_days, float daily_vals[]);
		void initializeFireModelConstantsAndOutputs();
		void snow_cond();
		void kbdi();
		void DFuelMC();
		void fuel_mc();
		void fwi();
		void min_mc_month();
		// void flammability();
		float liveFuelMoistureContent(float stress, float min_mc, float max_mc);
		void FuelLoad(int vtype, float clai, int month_length);
		// void erc_G(int month_length);
		void get_daily(float prev_mo_val, float curr_mo_val, int mo_len, float daily_vals[]);
		void FireOccur(int mo, int yrs_since_fire, int vtype, int * fire_domP, int * fire_doyP);
		void ros_dom(int dom, int doy, float * rosP, float * sgbrtP); 
		void erc_dom(int dom, int doy, float * ercP, float * tauP, float *ireP);
		void fuel_characteristics(int vtype);
		float part_burn(int doy, int dom, int burn_count, bool * suppressedFireFlagP, int vtype);
		void fire_return_interval(int vtype, float * min_friP, float * max_friP);
		void FireEffect(int fire_dom, int fire_doy, int vtype);
		float emfac(int vtype, EmissionSpeciesEnum em_species, EmissionSourceEnum em_source);
		void ann_effect(float cen_outvars[600], float cen_state[200], float part_burn);

	public:
		float c_stl, c_std, c_etasd, c_etasl;
		float m_mc_duff, m_mc_grass, m_mc_tree;
		FireParam fireParams; 
		float m_lbranch;   
		float m_lgras; // g DM m-2
		float m_lleaf; // g DM m-2
		float m_lwod1; // g DM m-2
		float m_lwod100; // g DM m-2   
		float m_dstnd; // g DM m-2
		float m_mlittr; // g DM m-2
		float m_slittr; // g DM m-2
		float m_dwod1; // g DM m-2
		float m_dwod100; // g DM m-2
		float m_d1hr, m_d10hr, m_d100hr, m_d1000hr;
		float m_fuel_depth;
		float m_fine_dead_fuel_yr;
		float m_coarse_dead_fuel_yr;
		int m_month_of_min_mc;
		float m_erc_g;
		float m_consumed[7];
		float m_killed[4];
		float m_consume_totbio, m_consume_live, m_consume_dead, m_blkc_totbio, m_death_totbio;
		float m_ffmc_ann_max;
		float m_bui_ann_max;
		float m_ffmc_threshold;
		float m_bui_threshold;
		float m_tree_ht_m; // current month's tree height in meters
		float m_snow_doy[365];

		// emissions from fire event
		float m_em_co2; // g m-2
		float m_em_co; // g m-2
		float m_em_ch4; // g m-2
		float m_em_nmhc; // g m-2
		float m_em_pm; // g m-2
		float m_em_pm2p5; // g m-2
		float m_em_pm10; // g m-2
		float m_blkc[7]; // black carbon
		float m_fire_fli; // fireline intensity in Btu/ft/sec
		float m_fire_ros, m_fire_erc;
		float m_bui_doy[365];
		float m_ffmc_doy[365];
		float m_fire_margin_doy[365];

	private:
		FuelDat fuel_data;
		float m_ppt[12]; // mm H2O
		float m_ppt_doy[365]; // mm H2O
		float m_ppt_doy_in[365]; // inches H2O
		float m_tmp[12]; // deg C
		float m_tmp_doy[365]; // deg C
		float m_tmp_doy_F[365]; // deg F
		float m_tmin[12];
		float m_tmin_doy_F[365]; // deg F
		float m_tmax[12];
		float m_tmax_doy_F[365]; // deg F
		float m_pet[12];
		// float m_pet_doy[365];
		float m_rh[12];
		float m_rh_doy[365];
		float m_rhmin[12];
		float m_rhmin_doy[365];
		float m_rhmax[12];
		float m_rhmax_doy[365];
		float m_wnd[12];
		float m_wnd_doy[365];
		float m_ppt_rat_in_per_hr[12];
		float m_ppt_rat_doy_in_per_hr[365];
		float m_kbdi_doy[365];
		float m_kbdi_max;
		float m_last_yr_snow;
		float m_tot_snowfall;
		float m_ave_ann_ppt;
		float m_dc_doy[365];
		float m_temp_corr_doy[365];
		float m_min_mc_yr;

		float m_fire_tau, m_fire_sgbrt, m_fire_ire;
		float m_crown_len_m; // current month's crown length in meters
		float m_bark_thick; // current month's bark thickness
		float m_w1p_dom[31]; // lbs per sq ft
		float m_w10_dom[31]; // lbs per sq ft
		float m_w100_dom[31]; // lbs per sq ft
		float m_w1000_dom[31]; // lbs per sq ft
		float m_wherbp_dom[31]; // lbs per sq ft
		// float m_wwood_dom[31]; // lbs per sq ft
		float m_depth_dom[31];
		float m_mc_grass_dom[31], m_mc_tree_dom[31];
		float m_fine_fuel_frac;

		// current fuel characteristics
		float m_sg1, m_sg10, m_sg100, m_sg1000, m_sgherb, m_sgwood, m_hl, m_hd;
		float m_mxd_frac, m_wndfac, m_mxl_frac;

}; // end of class MC_FireModel


