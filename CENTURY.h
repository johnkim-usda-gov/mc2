/*
 *  CENTURY.h
 *  mc2
 */

// forward declarations
class MC_FireModel;
class MC_BiogeographyModel;


class CENTURY_BiogeochemModel: public ProcessModel
{
  public:
    CENTURY_BiogeochemModel() {};
    CENTURY_BiogeochemModel(Simulation * pRun);
    ~CENTURY_BiogeochemModel() {};
    
    bool runModelAndCollectData(const int year_index, const int row_index, const int col_index);
    
    bool century_eq_model();
    bool mc2_spinup();
    bool mc2_transient();
    bool simulate1pt1yr(
         MC_FireModel * pF,
         MC_BiogeographyModel * pB,
         VCategory * prev_vclassP,
         int * vtypeP,
         TreeType * tree_typP,
         int yrNdx,
         int climateYrNdx,
         float	cen_state[CEN_STATE_SIZE],
         float tmp_smoothed[12],
         float ppt_smoothed[12],
         float initial_tmp_index,
         float initial_ppt_index);
    bool simulate1pt1yrInner(
         MC_FireModel * pF,
         MC_BiogeographyModel * pB,
         VCategory * prev_vclassP,
         int * vtypeP,
         TreeType * tree_typP,
         int yrNdx,
         float	cen_state[CEN_STATE_SIZE],
         int mapss_climate_zone,
         float mix_index,
         float tmp_index,
         float ppt_index,
         float ppt_smoothed[],
         float tmp_smoothed[],
         float ppt_raw[],
         float tmp_raw[],
         float tmax_raw[],
         float tmin_raw[],
         float vpr_raw[],
         float rh_raw[]);
    int vemap_aggregate_vegtype(unsigned int vtype, BaseCalibrationEnum biogeog);
    
  private:
    float m_min_smoothed_tmp; // lowest tmp_smoothed of the year, deg C
    float m_max_smoothed_tmp; // highest tmp_smoothed of the year, deg C
    float m_cont_index; // continentality index, deg C
    ClimateZone m_zone4biogeog; // climate zone for the biogeography model
    float m_gdd_zero; // growing degree-days, referenced to 0 deg C
    float m_nidx; // needleleaf-broadleaf index
    float m_eidx; // evergreen-deciduous index
    float m_c3pct; // fraction of total production which is from C3 photosynthesis
    float m_smoothed_max_tree_lai; // smoothed annual max woody LAI
    float m_smoothed_max_grass_lai; // smoothed annual max herbaceous LAI
    float m_burn_count; // number of years since last simulated fire
    float m_npp_yr; // net primary production, g C m-2 yr-1
    float m_npp_mo[12]; // net primary production, g C m-2 mo-1
    float m_nep_yr; // net ecosystem production, g C m-2 yr-1
    float m_nbp_yr; // net biome production, g C m-2 yr-1
    float m_rsp_yr; // heterotrophic respiration, g C m-2 yr-1
    float m_bio_consume_century; // C in biomass consumed by fire, from CENTURY, g C m-2 yr-1
    float m_aglivc_mo[12]; // monthly aboveground live grass carbon, g C m-2
    float m_bglivc_mo[12]; // monthly belowground live grass carbon, g C m-2
    float m_frstc_mo[12]; // monthly live forest carbon, g C m-2
    float m_aflivc_mo[12]; // monthly aboveground live forest carbon, g C m-2
    float m_bflivc_mo[12]; // monthly belowground live forest carbon, g C m-2
    float m_ddecid_mo[12]; // monthly drought deciduous index
    float m_vegc_mo[12];
    float m_fprd_ppt_mo[12];
    float m_gprd_ppt_mo[12];
    float m_fprd_tmp_mo[12];
    float m_gprd_tmp_mo[12];
    float m_aet_mo[12]; // monthly actual evapotranspiration, cm H2O
    float m_max_gfrac; // max monthly grass fraction of total live veg C, fraction
    bool m_no_soil_flag;
    bool m_fire_last_month;
    int m_nlayer;
    float m_frost_index_this_year;
    float m_tmmin; // lowest of the unsmoothed mean monthly temperatures for the year, deg C
    
    float m_wnd[12]; // mean monthly wind with no interannual variation
}; // end of class CENTURY_BiogeochemModel


