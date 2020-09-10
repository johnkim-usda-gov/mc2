/*
 *  MCbiogeog.h
 *  mc2
 */

struct BiogeographyInputData
{
	public:
		TreeType tree_typ;
		float npp_yr;
		float mean_treec;
		float max_grassc;
		float max_grass_frac;
		float mean_vegc;
		ClimateZone zone;
		float c3pct;
		float gdd0;
		float cont_index;
		float min_smoothed_tmp;
		float ppt_yr;
		float tmp_yr;
		float elev;
		float tot_summer_ppt; // mmH2O, JJA in N hemi, 92/90 x DJF in S hemi

		// next batch is for W_WA calibration
		float fog; // fog factor
		float topomoist; // topographic moisture factor
		float deltaTsl; // difference between mean air temperature and mean air temperature at sea level, deg F
		float aspect; // aspect
		float sw; // shortwave
		float cad; // cold air drainage
		int subregion; // modeling subregion ("ecoregion" sensu Henderson PNV model)

		// next 8 are for Land Carbon biogeography 
		float * ddecid_mo;
		float tree_lai;
		float grass_lai;
		float mean_grass_frac;
		float aflivc;
		float tmmin;
		int nlayer;
		int lulcType; // added for land-use processing

};


class MC_BiogeographyModel: public ProcessModel
{
	MC_BiogeographyModel() {}
	public:
	MC_BiogeographyModel(Simulation * sP, RunParamsClass * rP, ModelParamsClass * mP); 
	~MC_BiogeographyModel() {}

	bool runModelAndCollectData(const int year_index, const int row_index, const int col_index); 

	bool BiogeogMC2(BiogeographyInputData inVals, 
			BiomeType * biomeP, PhysiognomicClass * physiognomic_classP, MC2VegType * vtypeP);
	MC2VegType BiogeogLC(BiogeographyInputData inVals);
	// void BiogeogUS50km();
	// void BiogeogLYNXandWWETAC();
	int GrouseModel(int vtype, float smrpre, float augmaxt, float anntmp); // Returns 0 for no grouse,
	// otherwise returns vtype.

	BiomeType m_biome;
	PhysiognomicClass m_physiognomic_class;
	ModelParamsClass * modelParamsP;
	BaseCalibrationEnum m_baseCalibration;
	float m_taiga_tundra_threshold;

}; // end of class MC_BiogeographyModel


