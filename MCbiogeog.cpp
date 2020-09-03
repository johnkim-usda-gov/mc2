/*
 *  MCbiogeog.cpp
 *  mc2
 *
 */

#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <math.h> 
#include <vector>

#include "netcdf.h"
#include "category_bgc.h" 
#include "assert.h"

#include "MAPSSvegClasses.h"
#include "ScienceFcns.h"

#include "ProcessModel.h"

#include "MC2.h"
#include "MAPSSbiogeographyModel.h"
#include "MCfire.h"
#include "MCbiogeog.h"
#include "CENTURY.h"


MC_BiogeographyModel::MC_BiogeographyModel(Simulation * sP, RunParamsClass * rP, ModelParamsClass * mP)
{ 
	ProcessModelInit(sP);
	MAPSS_model_parameters MAPSSparams(mP->MAPSSparameterSet, 1.0f);

	modelParamsP = mP; 
	runParamsP = rP;
	m_taiga_tundra_threshold = MAPSSparams.taiga_tundra_boundary[mapssBOREAL];

	m_psl = m_SSZupperBound = m_PSFZlowerBound = m_SAFZlowerBound = -9999.f;

} // end of constructor for MC_BiogeographyModel


bool MC_BiogeographyModel::runModelAndCollectData(const int year_index, const int row_index, const int col_index) 
{
	return(true);
}; // end of MC_BiogeographyModel::runModelAndCollectData()


void MC_BiogeographyModel::BiogeogMC2(BiogeographyInputData inVals, 
		BiomeType * biomeP, PhysiognomicClass * physiognomic_classP, MC2VegType * vtypeP)
{
	ModelParamsClass * mP = modelParamsP;

	MC2VegType vtype;

	C3C4Dominance grass_typ;

	const float subalp_thres = modelParamsP->subalpine_threshold;      /* UPPER GDD LIMIT FOR SUBALPINE */
	const float mari_thres = modelParamsP->maritime_threshold;          /* UPPER LIMIT OF CONTINENTAL INDEX FOR MARITIME FOREST */

	const float forest_thres = modelParamsP->m_forest_thres_C;         /* LOWER LIMIT OF TOTAL WOODY C FOR FOREST */
	const float woodl_thres  = modelParamsP->m_woodl_thres_C;         /* LOWER LIMIT OF TOTAL WOODY C FOR WOODLAND, g C m-2 */
	const float c3_threshold = modelParamsP->c3_threshold; // % of total from C3 photosynthesis
	const float grassfrac_thres = modelParamsP->grassfrac_thres; // frac of live carbon
	const float desert_treec_max = 27.; // Classify as desert when both treec<desert_treec_max
	const float desert_grassc_max = modelParamsP->desert_grass_C_max; // and grassc<desert_grassc_max.
	const float unvegetated_thres = 5; // g C m-2
	const float savanna_thres = 0.10; // In the GRASSLAND_SAVANNA biome, when woody_frac>=savanna_thres, it's a SAVANNApclass.
	const float taiga_tundra_C_min = 400.0f; // gC m-2

	TreeType tree_typ = inVals.tree_typ;
	float treec = inVals.mean_treec; // mean live tree carbon
	float max_grassc = inVals.max_grassc; // max live grass carbon
	float gdd_zero = inVals.gdd0;
	float cont_index = inVals.cont_index;
	ClimateZone zone = inVals.zone;
	float c3 = inVals.c3pct;
	float max_grass_frac = inVals.max_grass_frac; 
	float min_woody_frac = 1. - max_grass_frac; // min fraction of live vegetation carbon which is in woody plants

	// Next batch are per US Forest Service PNW-GTR-841, for WCR base calibration
	float psl = NC_FILL_FLOAT; // precipitation at sea level, mmH2O 
	float fog_effect = NC_FILL_FLOAT; // fog effect, each unit = 20"H2O
	float fog = NC_FILL_FLOAT; // fog effect, mmH2O
	float topomoist = NC_FILL_FLOAT; // topographic moisture factor
	float matsl = NC_FILL_FLOAT; // mean air temperature at sea level, deg C
	float aspect = NC_FILL_FLOAT; // aspect, compass degrees
	float sw = NC_FILL_FLOAT; // shortwave, mean daily clear-sky shortwave radiation, in kilojoules per square meter per day
	float cad_ft = NC_FILL_FLOAT; // cold air drainage factor, ft
	float cad = NC_FILL_FLOAT; // cold air drainage factor, m
	double elev = inVals.elev; // elevation, m
	double pslM = 3.14066e-7; // M
	double pslN = 0.99908387; // N
	double pslP = -139015.01; // P
	double pslQ = 1922.4891; // Q
	double pslR = 561677190.; // R
	double pslS = -13825146.; // S

	if (runParamsP->baseCalibration==mc2WCR)
	{ 
		assert(inVals.fog>-100.f);
		fog_effect = inVals.fog; 
		fog = sciFn.in_to_mm(20.f*fog_effect);

		assert(inVals.topomoist>-100.f);
		topomoist = inVals.topomoist;

		assert(inVals.deltaTsl>-100.f);
		matsl = inVals.tmp_yr + inVals.deltaTsl;

		assert(inVals.aspect>-100.f);
		aspect = inVals.aspect; // aspect, deg

		assert(inVals.sw>-100.f); 
		sw = inVals.sw; 

		assert(inVals.cad>-100.f);
		cad_ft = inVals.cad;  
		cad = sciFn.ft_to_m(cad_ft);  

		double tap = inVals.ppt_yr; // total annual precipitation, mmH2O
		psl = (float)(1./(pslM + pslN/tap + (1./(pslP + pslQ*tap))*elev + (1./(pslR + pslS*tap))*elev*elev));
		assert(psl>=0. && psl<10000.);
		m_psl = psl;

		float pslfog = psl + fog; 
		float sine_term = sin((aspect - 120.f)*(PI/180.f));
		float sw_term = (sw - 14930.f)/3380.4f;

		m_SSZupperBound = // eq. E2 in Table 3, p. 25 of PNW-GTR-841
			// A + B*(PSL + FOG) + C*TM
			mP->SSZ_ub_coeffs[0] + mP->SSZ_ub_coeffs[1]*pslfog + mP->SSZ_ub_coeffs[2]*topomoist
			// + D*MATSL + G*FOG + H*CAD
			+ mP->SSZ_ub_coeffs[3]*matsl + mP->SSZ_ub_coeffs[6]*fog + mP->SSZ_ub_coeffs[7]*cad
			// + {P*E*SIN(ASPECT - 120) + [(1 - P)*E*(SW - 14930)/3380.4]}
			+ (mP->SSZ_ub_coeffs[8]*mP->SSZ_ub_coeffs[4]*sine_term
					+ ((1.f - mP->SSZ_ub_coeffs[8])*mP->SSZ_ub_coeffs[4]*sw_term));
		if (m_SSZupperBound<-10000.f || m_SSZupperBound>10000.f)
			assert(0);

		m_PSFZlowerBound = // eq. E5
			// A + B*(PSL + FOG) + C*TM
			mP->PSFZ_lb_coeffs[0] + mP->PSFZ_lb_coeffs[1]*pslfog + mP->PSFZ_lb_coeffs[2]*topomoist 
			// + D*MATSL + G*FOG + H*CAD
			+ mP->PSFZ_lb_coeffs[3]*matsl + mP->PSFZ_lb_coeffs[6]*fog + mP->PSFZ_lb_coeffs[7]*cad
			// + {P*[E + F*(PSL + FOG)]*SIN(ASPECT - 120) + [(1 - P)*E*(SW - 14930)/3380.4]}
			+ (mP->PSFZ_lb_coeffs[8]*(mP->PSFZ_lb_coeffs[4] + mP->PSFZ_lb_coeffs[5]*pslfog)*sine_term
					+ ((1.f - mP->PSFZ_lb_coeffs[8])*mP->PSFZ_lb_coeffs[4]*sw_term));

		m_SAFZlowerBound = // eq. E1
			// A + B*(PSL + FOG) + C*TM
			mP->SAFZ_lb_coeffs[0] + mP->SAFZ_lb_coeffs[1]*pslfog + mP->SAFZ_lb_coeffs[2]*topomoist
			// + D*MATSL 
			+ mP->SAFZ_lb_coeffs[3]*matsl 
			// + {P*E*SIN(ASPECT - 120) + [(1 - P)*E*(SW - 14930)/3380.4]}
			+ (mP->SAFZ_lb_coeffs[8]*mP->SAFZ_lb_coeffs[4]*sine_term
					+ ((1.f - mP->SAFZ_lb_coeffs[8])*mP->SAFZ_lb_coeffs[4]*sw_term));

		m_PKLZlowerBound = // eq. E1
			// A + B*(PSL + FOG) + C*TM
			mP->PKLZ_lb_coeffs[0] + mP->PKLZ_lb_coeffs[1]*pslfog + mP->PKLZ_lb_coeffs[2]*topomoist
			// + D*MATSL 
			+ mP->PKLZ_lb_coeffs[3]*matsl 
			// + {P*E*SIN(ASPECT - 120) + [(1 - P)*E*(SW - 14930)/3380.4]}
			+ (mP->PKLZ_lb_coeffs[8]*mP->PKLZ_lb_coeffs[4]*sine_term
					+ ((1.f - mP->PKLZ_lb_coeffs[8])*mP->PKLZ_lb_coeffs[4]*sw_term));
	}
	pS->VarDict[MATSL].save(matsl, YEAR_INTERVAL);

	grass_typ = c3>=c3_threshold ? C3Dominance : C3C4Mixed; /* CLASSIFY GRASS LIFEFORM */

	/* Biogeography rule set */

	// Determine biome: one of desert, shrubland, grassland/savanna, woodland, or forest.
	m_biome = UNKNOWNbiome;
	if (inVals.npp_yr <= 0. || (treec < desert_treec_max && max_grassc < desert_grassc_max)) 
		m_biome = DESERTbiome;
	else if (treec < woodl_thres) 
	{
		if (max_grass_frac < grassfrac_thres)
			m_biome = SHRUBLANDbiome;
		else
			m_biome = GRASSLAND_SAVANNAbiome;
	}
	else if (treec < forest_thres) 
		m_biome = WOODLANDbiome;
	else 
		m_biome = FORESTbiome;

	// Break down the biomes into physiognomic classes: one of 
	// unvegetated land, semi-desert grassland, semi-desert shrubland, grassland, savanna, 
	// 3 kinds of woodland (evergreen, deciduous, and mixed deciduous/evergreen)
	// 3 kinds of forest (evergreen, deciduous, and mixed deciduous/evergreen)
	switch (m_biome)
	{
		case DESERTbiome: 
			if (inVals.mean_vegc<unvegetated_thres) m_physiognomic_class = UNVEGETATEDpclass;
			else m_physiognomic_class = max_grass_frac>=grassfrac_thres ? SEMIDESERT_GRASSLANDpclass : SEMIDESERT_SHRUBLANDpclass;
			break;
		case SHRUBLANDbiome: m_physiognomic_class = SHRUBLANDpclass; break;
		case GRASSLAND_SAVANNAbiome: m_physiognomic_class = min_woody_frac>savanna_thres ? SAVANNApclass : GRASSLANDpclass; break;
		case WOODLANDbiome: switch (tree_typ)
					    {
						    case EN_TREES:
						    case EB_TREES:
						    case EN_EB_TREES:
							    m_physiognomic_class = EVERG_WOODLANDpclass;
							    break;
						    case DN_TREES:
						    case DB_TREES:
							    m_physiognomic_class = DECID_WOODLANDpclass;
							    break;
						    case DN_EN_TREES:
						    case DB_EB_TREES:
						    case EN_DB_TREES:
							    m_physiognomic_class = MIXED_WOODLANDpclass;
							    break;
						    default: assert(0); break;
					    }
				    break;
		case FORESTbiome: switch (tree_typ)
					  {
						  case EN_TREES:
						  case EB_TREES:
						  case EN_EB_TREES:
							  m_physiognomic_class = EVERG_FORESTpclass;
							  break;
						  case DN_TREES:
						  case DB_TREES:
							  m_physiognomic_class = DECID_FORESTpclass;
							  break;
						  case DN_EN_TREES:
						  case DB_EB_TREES:
						  case EN_DB_TREES:
							  m_physiognomic_class = MIXED_FORESTpclass;
							  break;
						  default: assert(0); break;
					  }
				  break;
		default: assert(0); break;
	} // end of switch (m_biome)

	// Now break down the physiognomic classes into potential vegetation types.   
	vtype = UNKNOWNveg;
	switch (m_physiognomic_class)
	{    
		case UNVEGETATEDpclass:
			if (inVals.npp_yr<=0.) vtype = NATURAL_BARRENveg; 
			else switch (zone)
			{
				case ARCTICzone: 
				case BOREALzone: vtype = COLD_BARRENveg; break;
				case TEMPERATEzone: vtype = TEMPERATE_DESERTveg; break;
				case SUBTROPICALzone: vtype = SUBTROPICAL_DESERTveg; break;
				case TROPICALzone: vtype = TROPICAL_DESERTveg; break;
				default: assert(0); break;
			}
			break;
		case SEMIDESERT_SHRUBLANDpclass: 
			switch (zone)
			{
				case ARCTICzone: 
				case BOREALzone: 
					vtype = TUNDRAveg; 
					break;
				case TEMPERATEzone: 
				case SUBTROPICALzone: vtype = c3>=c3_threshold ? C3SHRUBveg : C4SHRUBveg; break;
				case TROPICALzone: assert(c3<c3_threshold); vtype = C4SHRUBveg; break;
				default: assert(0); break;
			}
			break;
		case SHRUBLANDpclass:
			switch (zone)
			{
				case ARCTICzone: 
				case BOREALzone: 
					if (gdd_zero>=m_taiga_tundra_threshold) vtype = C3SHRUBveg;  
					else vtype = (inVals.mean_vegc>=taiga_tundra_C_min) ? TAIGA_TUNDRAveg : TUNDRAveg; 
					break;
				case TEMPERATEzone: 
				case SUBTROPICALzone: vtype = c3>=c3_threshold ? C3SHRUBveg : C4SHRUBveg; break;
				case TROPICALzone: assert(c3<c3_threshold); vtype = C4SHRUBveg; break;
				default: assert(0); break;
			}
			break;
		case SEMIDESERT_GRASSLANDpclass:
			switch (zone)
			{
				case ARCTICzone: 
				case BOREALzone: vtype = TUNDRAveg; break;
				case TEMPERATEzone: vtype = TEMPERATE_DESERTveg; break;
				case SUBTROPICALzone: vtype = SUBTROPICAL_DESERTveg; break;
				case TROPICALzone: assert(c3<c3_threshold); vtype = TROPICAL_DESERTveg; break;
				default: assert(0); break;
			}
			break;
		case GRASSLANDpclass:
			switch (zone)
			{
				case ARCTICzone: 
				case BOREALzone: vtype = TUNDRAveg; break;
				case TEMPERATEzone: 
				case SUBTROPICALzone: vtype = c3>=c3_threshold ? C3GRASSveg : C4GRASSveg; break;
				case TROPICALzone: assert(c3<c3_threshold); vtype = C4GRASSveg; break;
				default: assert(0); break;
			}
			break;
		case SAVANNApclass:
			switch (zone)
			{
				case ARCTICzone: 
				case BOREALzone: 
					if (gdd_zero>=m_taiga_tundra_threshold) vtype = C3GRASSveg;  
					else vtype = TAIGA_TUNDRAveg; 
					break;
				case TEMPERATEzone: 
				case SUBTROPICALzone: vtype = c3>=c3_threshold ? C3GRASSveg : C4GRASSveg; break;
				case TROPICALzone: vtype = TROPICAL_SAVANNAveg; break;
				default: assert(0); break;
			}
			break;
		case EVERG_WOODLANDpclass:
			switch (zone)
			{
				case ARCTICzone: 
				case BOREALzone: vtype = gdd_zero>=subalp_thres ? TEMPERATE_EN_WOODLANDveg : BOREAL_WOODLANDveg; break;
				case TEMPERATEzone: 
						 /*
						    if (tree_typ==EN_TREES) vtype = 
						    inVals.ppt_yr>mP->dry_temperate_threshold ? TEMPERATE_EN_WOODLANDveg : XERIC_NEEDLELEAF_WOODLANDveg;
						    else vtype = TEMPERATE_WARM_MIXED_WOODLANDveg;
						  */
						 vtype = tree_typ==EN_TREES ? TEMPERATE_EN_WOODLANDveg : TEMPERATE_WARM_MIXED_WOODLANDveg;
						 break;
				case SUBTROPICALzone: vtype = SUBTROPICAL_EB_WOODLANDveg; break;
				case TROPICALzone: vtype = max_grass_frac >= grassfrac_thres ? TROPICAL_SAVANNAveg : TROPICAL_EB_FORESTveg; break;
				default: assert(0); break;
			}
			break;
		case DECID_WOODLANDpclass:
			switch (zone)
			{
				case ARCTICzone: 
				case BOREALzone: vtype = gdd_zero>=subalp_thres ? TEMPERATE_DB_WOODLANDveg : BOREAL_WOODLANDveg; break;
				case TEMPERATEzone: vtype = TEMPERATE_DB_WOODLANDveg; break;
				case SUBTROPICALzone: vtype = SUBTROPICAL_DB_WOODLANDveg; break;
				case TROPICALzone: vtype = TROPICAL_DECIDUOUS_WOODLANDveg; break;
				default: assert(0); break;
			}
			break;
		case MIXED_WOODLANDpclass:
			switch (zone)
			{
				case ARCTICzone: 
				case BOREALzone: vtype = gdd_zero>=subalp_thres ? TEMPERATE_COOL_MIXED_WOODLANDveg : BOREAL_WOODLANDveg; break;
				case TEMPERATEzone: 
						 vtype = tree_typ==EN_DB_TREES ? TEMPERATE_COOL_MIXED_WOODLANDveg : TEMPERATE_WARM_MIXED_WOODLANDveg; 
						 break;
				case SUBTROPICALzone: vtype = SUBTROPICAL_DB_WOODLANDveg; break;
				case TROPICALzone: vtype = TROPICAL_DECIDUOUS_WOODLANDveg; break;
				default: assert(0); break;
			}
			break;
		case EVERG_FORESTpclass:
			switch (zone)
			{
				case ARCTICzone: 
				case BOREALzone: 
					switch (tree_typ)
					{
						case EN_TREES: 
						case DN_EN_TREES:
						case DN_TREES: 
							if (gdd_zero<subalp_thres) vtype = BOREAL_NEEDLELEAF_FORESTveg;
							else if (inVals.ppt_yr>mP->moist_temperate_threshold) vtype = MOIST_TEMPERATE_NEEDLELEAF_FORESTveg;
							else if (inVals.ppt_yr>mP->dry_temperate_threshold) vtype = MESIC_TEMPERATE_NEEDLELEAF_FORESTveg;
							else vtype = DRY_TEMPERATE_NEEDLELEAF_FORESTveg;
							break;
						case EN_DB_TREES: vtype = COOL_MIXED_FORESTveg; break;
						default: assert(0); break;
					}
					break;
				case TEMPERATEzone: 
					switch (tree_typ)
					{
						case EN_TREES:
							if (gdd_zero<=subalp_thres) 
							{ // SUBALPINE_FORESTveg; /* Subalpine */
								if (runParamsP->baseCalibration==mc2WCR)
								{
									if (elev<m_SAFZlowerBound) vtype = MHZveg;
									else if (elev<m_PKLZlowerBound) vtype = SAFZveg;
									else vtype = PKLZveg;
								}
								else vtype = SUBALPINE_FORESTveg;
							}
							/*
							   else if (cont_index<=mari_thres) 
							   {              

							   if (inVals.min_smoothed_tmp<modelParamsP->tmmin_threshold)
							   {
							   if (inVals.ppt_yr>mP->moist_temperate_threshold) vtype = MOIST_TEMPERATE_NEEDLELEAF_FORESTveg;
							   else if (inVals.ppt_yr>mP->dry_temperate_threshold) vtype = MESIC_TEMPERATE_NEEDLELEAF_FORESTveg;
							   else vtype = DRY_TEMPERATE_NEEDLELEAF_FORESTveg;
							   }
							   else
							   { // MARITIME_EN_FORESTveg;
							   if (runParamsP->baseCalibration==mc2WCR)
							   {
							   if (elev<=m_SSZupperBound) vtype = SSZveg;
							   else if (elev<m_PSFZlowerBound) vtype = WHZveg;
							   else vtype = PSFZveg;
							   }
							   else vtype = MARITIME_EN_FORESTveg;
							   }
							   }
							   else vtype = MESIC_TEMPERATE_NEEDLELEAF_FORESTveg; 
							   break;
							 */
							else if (cont_index<=mari_thres && inVals.min_smoothed_tmp>=modelParamsP->tmmin_threshold) 
							{ // MARITIME_EN_FORESTveg;
								if (runParamsP->baseCalibration==mc2WCR)
								{
									if (elev<=m_SSZupperBound) vtype = SSZveg;
									else if (elev<m_PSFZlowerBound) vtype = WHZveg;
									else vtype = PSFZveg;
								}
								else vtype = MARITIME_EN_FORESTveg;
							}
							else if (inVals.ppt_yr>mP->moist_temperate_threshold) vtype = MOIST_TEMPERATE_NEEDLELEAF_FORESTveg;
							else if (inVals.ppt_yr>mP->dry_temperate_threshold) vtype = MESIC_TEMPERATE_NEEDLELEAF_FORESTveg;
							else vtype = DRY_TEMPERATE_NEEDLELEAF_FORESTveg;
							break;              
						case EB_TREES: vtype = WARM_EB_FORESTveg; break;
						case EN_EB_TREES: vtype = TEMPERATE_WARM_MIXED_FORESTveg; break;
						default: assert(0); break;
					}
					break;
				case SUBTROPICALzone: 
					switch (tree_typ)
					{
						case EN_TREES: vtype = cont_index<=mari_thres ? MARITIME_EN_FORESTveg : SUBTROPICAL_EN_FORESTveg; break;
						case EN_EB_TREES: vtype = SUBTROPICAL_MIXED_FORESTveg; break;
						case EB_TREES: vtype = WARM_EB_FORESTveg; break;
						default: assert(0); break;
					}
					break;
				case TROPICALzone: 
					assert(tree_typ!=DN_TREES && tree_typ!=DN_EN_TREES);
					vtype = TROPICAL_EB_FORESTveg; 
					break;
				default: assert(0); break;
			}
			break;
		case DECID_FORESTpclass:
			switch (zone)
			{
				case ARCTICzone: 
				case BOREALzone: 
					switch (tree_typ)
					{
						case DN_EN_TREES:
						case DN_TREES: vtype = LARCH_FORESTveg; break;
						case EN_DB_TREES: vtype = COOL_MIXED_FORESTveg; break;
						case DB_TREES: vtype = TEMPERATE_DB_FORESTveg; break;
						default: assert(0); break;
					}
					break;
				case TEMPERATEzone: 
					switch (tree_typ)
					{
						case DB_TREES: vtype = TEMPERATE_DB_FORESTveg; break;
						case EN_DB_TREES: vtype = COOL_MIXED_FORESTveg; break;
						case DB_EB_TREES: vtype = TEMPERATE_WARM_MIXED_FORESTveg; break;
						default: assert(0); break;
					}
					break;
				case SUBTROPICALzone: 
					assert(tree_typ==DB_EB_TREES);
					vtype = SUBTROPICAL_DB_FORESTveg; 
					break;
				case TROPICALzone: 
					assert(tree_typ==DB_EB_TREES);
					vtype = TROPICAL_EB_FORESTveg; 
					break;
				default: assert(0); break;
			}
			break;
		case MIXED_FORESTpclass:
			switch (zone)
			{
				case ARCTICzone: 
				case BOREALzone: 
					switch (tree_typ)
					{
						case DN_EN_TREES: vtype = LARCH_FORESTveg; break; // does happen at global cell row 33, col 591
						case EN_DB_TREES: vtype = COOL_MIXED_FORESTveg; break;
						default: assert(0); break;
					}
					break;
				case TEMPERATEzone: 
					switch (tree_typ)
					{
						case EN_DB_TREES: vtype = COOL_MIXED_FORESTveg; break;
						case DB_EB_TREES:
						case EN_EB_TREES: vtype = TEMPERATE_WARM_MIXED_FORESTveg; break;
						default: assert(0); break;
					}
					break;
				case SUBTROPICALzone: 
					switch (tree_typ)
					{
						case EN_DB_TREES: 
						case EN_EB_TREES: vtype = SUBTROPICAL_MIXED_FORESTveg; break;
						case DB_EB_TREES: vtype = SUBTROPICAL_DB_FORESTveg; break;
						default: 
								  printf("*** BiogeogMC2(): tree_typ = %d in SUBTROPICALzone\n", tree_typ);
								  assert(0);              
								  break;
					}
					break;
				case TROPICALzone: 
					assert(tree_typ!=DN_TREES && tree_typ!=DN_EN_TREES);
					vtype = TROPICAL_EB_FORESTveg; 
					break;
				default: assert(0); break;
			}
			break;
		default: assert(0); break;
	} // end of switch (m_physiognomic_class)

	assert(vtype!=UNKNOWNveg);  
	assert(vtype<=MAX_VTYPE);

	*biomeP = m_biome;
	*physiognomic_classP = m_physiognomic_class;
	*vtypeP = vtype;

} // end of BiogeogMC2()


#define LYNX_BIOGEOG 2
#define WWETAC_BIOGEOG 5

MC2VegType MC_BiogeographyModel::BiogeogLC(const BiogeographyInputData inVals)
{
	MC2VegType     vtype;
	// Check to see if land use dictates overriding computing biogeography
	if (inVals.lulcType) {
		vtype = Lynx_Undefined;
		switch(inVals.lulcType) {
			case LULC_Undefined:
			case LULC_Default:
			case LULC_Mechanical_Disturb:
				// Do nothing
				break;
			case LULC_Agriculture:
				vtype = Lynx_AgricultureGrazing;
				break;
			case LULC_Developed:
				vtype = Lynx_Developed;
				break;
			case LULC_Mining:
				vtype = Lynx_Mining;
				break;
			default:
				// do nothing
				break;
		} // switch(inVals.lulcType) {

		if(vtype != Lynx_Undefined)
			return vtype;
	}


	// Use WWETAC biogeography
	//const int biogeog_option = WWETAC_BIOGEOG;
	const int biogeog_option = LYNX_BIOGEOG;

	/* ZONE THRESHOLDS */
	// Don't need these, as zone is supplied in inVals
	//      const float     az_thres = 1000.;     /* UPPER GDD LIMIT FOR ARCTIC aka
	//                                             * ALPINE ZONE */
	//      const float     bz_thres = -13.0;     /* UPPER MIN TEMP LIMIT FOR BOREAL
	//                                             * ZONE */
	//      const float     tz_thres = 7.75;      /* UPPER MIN TEMP LIMIT FOR TEMPERATE
	//                                             * ZONE */
	//      const float     stz_thres = 18.0;     /* UPPER MIN TEMP LIMIT FOR
	//                                             * SUBTROPICAL ZONE */
	//    
	/* VTYPE THRESHOLDS */

	const float     tt_thres = 1330.;     /* UPPER GDD LIMIT FOR TAIGA-TUNDRA */
	const float     subalpine_thres = modelParamsP->subalpine_threshold; /* UPPER GDD LIMIT FOR
									      * SUBALPINE */
	const float     mari_thres = modelParamsP->maritime_threshold;        /* UPPER LIMIT OF
									       * CONTINENTAL INDEX FOR
									       * MARITIME FOREST */
	const float     ddecid_thres = .45;   /* UPPER LIMIT OF DROUGHT-DECID INDEX
					       * FOR TROPICAL WOODLAND */

	const float     forest_thres = modelParamsP->m_forest_thres_C;        /* LOWER LIMIT OF TOTAL
									       * WOODY C FOR FOREST */
	const float savanna_threshold = 1150.;
	const float woodl_thres = 1150.;

	// tjs 2012.11.26 changing threshold for temperate region
	// this should prevent temperate deserts.
	const float     shrubl_thres_temperate = 1.0;   /* LOWER LIMIT OF TOTAL WOODY C FOR
							 * TEMPERATE SHRUBLAND */
	const float     shrubl_thres = 5.0;   /* LOWER LIMIT OF TOTAL WOODY C FOR
					       * SHRUBLAND */
	const float     grass_thres = 200.0;  /* LOWER LIMIT OF TOTAL GRASS C FOR
					       * GRASSLAND */
	const float     tree_grass_thres = 80.0;      /* Upper Limit of Tree C for
						       * Grassland in the temperate
						       * zone */
	const float     c3_threshold = modelParamsP->c3_threshold;
	//%of total from C3 photosynthesis

	const float     grassfrac_thres = modelParamsP->grassfrac_thres;
	//frac of live carbon

	const float tmmin_threshold = modelParamsP->tmmin_threshold;

	/* ASSIGNMENTS FROM DATAPOINT */

	//needle = data_point->vemap2_mo.nidx;
	//everg = data_point->vemap2_mo.eidx;

	const float     gdd_zero = inVals.gdd0;

	const float     treec = inVals.mean_treec;

	const float grassc = inVals.max_grassc;

	const float npp = inVals.npp_yr;

	const float aflivc = inVals.aflivc; //mean monthly aboveground live forest carbon, g C m - 2

	// Tree lifeform is now passed in via inVals, so we don 't need to calculate it here

	TreeType tree_typ = inVals.tree_typ;

	/* CLASSIFY GRASS LIFEFORM */

	C3C4Dominance   grass_typ;
	const float     c3 = inVals.c3pct;

	if (biogeog_option == WWETAC_BIOGEOG)
	{ if (c3 >= c3_threshold)
		grass_typ = C3Dominance;
		else
			grass_typ = C4Dominance;
	}
	else
	{ if (c3 > 66.)
		grass_typ = C3Dominance;
		else if (c3 < 33.)
			grass_typ = C4Dominance;
		else
			grass_typ = C3C4Mixed;
	}

	/* Calculate average tropical drought-decid index */
	float           sum = 0.;
	for (int mo = 0; mo < 12; mo++)
		sum += inVals.ddecid_mo[mo];
	const float     ddecid = sum / 12.;

	/* CALCULATE CONTINENTAL INDEX */

	const float     cont_index = inVals.cont_index;

	/* DETERMINE ZONE */

	//This is already in inVals...

	const ClimateZone zone = inVals.zone;

	vtype = Lynx_Undefined;


	if (biogeog_option == WWETAC_BIOGEOG
			// 2012.12.27 tjs dealing with temperate Lynx_Subtropical_Shrubland in temperate section
			&& zone != TEMPERATEzone
			//
			&& gdd_zero > tt_thres
			&& aflivc < savanna_threshold
			&& grass_typ != C3Dominance
			&& inVals.mean_grass_frac < grassfrac_thres)
		// 2012.12.21 tjs getting rid of coniferous xeromorphic woodland and
		// using subtropical as C4
		//vtype = Lynx_Coniferous_Xeromorphic_Woodland;
		vtype = Lynx_Subtropical_Shrubland;
	// tjs 2013.01.14 disabling subalpine meadow
#if 0
	else if (biogeog_option == WWETAC_BIOGEOG
			&& tt_thres < gdd_zero
			&& gdd_zero <= subalpine_thres
			&& tree_typ == EN_TREES
			&& tree_lai < subalpine_meadow_thres)
		vtype = Lynx_Subalpine_Meadow;
#endif
	else
		switch (zone)
		{ case UNKNOWNzone:
			vtype = Lynx_Undefined;
			break;

			case ARCTICzone:
			if (gdd_zero <= 0.)
				vtype = Lynx_Barren;
			else
				vtype = Lynx_Tundra;
			break;

			case BOREALzone:
			if (gdd_zero <= tt_thres)
				vtype = Lynx_Taiga_Tundra;

			// tjs 2013.01.04 Grassland into boreal zone
			else if (grassc >= grass_thres && treec <= tree_grass_thres)
			{ if (biogeog_option == WWETAC_BIOGEOG && grass_typ != C3Dominance)
				vtype = Lynx_Subtropical_Grassland;
				else
					vtype = Lynx_Temperate_Grassland;
			}
			//

			else if (treec >= forest_thres)
				vtype = Lynx_Boreal_Evergreen_Needleleaf_Forest;
			else
				vtype = Lynx_Boreal_Mixed_Woodland;
			break;

			case TEMPERATEzone:
			if (gdd_zero <= subalpine_thres)
				vtype = Lynx_Subalpine;
			else if (grassc >= grass_thres && treec <= tree_grass_thres)
			{ if (biogeog_option == WWETAC_BIOGEOG && grass_typ != C3Dominance)
				vtype = Lynx_Subtropical_Grassland;
				else
					vtype = Lynx_Temperate_Grassland;
			}
			else if (treec >= forest_thres)
			{ if (tree_typ == EN_TREES)
				{ if (cont_index <= mari_thres)
					{ if (biogeog_option == WWETAC_BIOGEOG && inVals.tmmin < tmmin_threshold)
						vtype = Lynx_Cool_Needleleaf_Forest;
						else
							vtype = Lynx_Maritime_Evergreen_Needleleaf_Forest;
					}
					else
						vtype =  Lynx_Temperate_Evergreen_Needleleaf_Forest;
				}
				else if (tree_typ == DB_TREES)
					vtype =  Lynx_Temperate_Deciduous_Broadleaf_Forest;
				else if (tree_typ == EN_DB_TREES)
					vtype = Lynx_Temperate_Cool_Mixed_Forest;
				else
					vtype = Lynx_Temperate_Warm_Mixed_Forest;
			}
			else if (treec >= woodl_thres)
			{ if (tree_typ == EN_TREES)
				vtype = Lynx_Temperate_Evergreen_Needleleaf_Woodland;
				else if (tree_typ == DB_TREES)
					vtype = Lynx_Temperate_Deciduous_Broadleaf_Woodland;
				else if (tree_typ == EN_DB_TREES)
					vtype = Lynx_Temperate_Cool_Mixed_Woodland;
				else
					vtype = Lynx_Temperate_Warm_Mixed_Woodland;
			}
			else if (treec >= shrubl_thres_temperate)

				// 2012.12.27 tjs moving subtropical shrubland into temperate section
				if (biogeog_option == WWETAC_BIOGEOG && grass_typ != C3Dominance)
					vtype = Lynx_Subtropical_Shrubland;
				else
					vtype = Lynx_Temperate_Shrubland;
			// previous code
			//vtype = Lynx_Temperate_Shrubland;

			else
				vtype = Lynx_Temperate_Desert;
			break;

			case SUBTROPICALzone:
			if (grassc >= grass_thres && treec <= grass_thres)
			{ if (biogeog_option == WWETAC_BIOGEOG && grass_typ == C3Dominance)
				vtype = Lynx_Temperate_Grassland;
				else
					vtype = Lynx_Subtropical_Grassland;
			}
			else if (treec >= forest_thres)
			{ if (tree_typ == EN_TREES)
				{ if (cont_index <= mari_thres)
					vtype = Lynx_Maritime_Evergreen_Needleleaf_Forest;
					else
						vtype = Lynx_Subtropical_Evergreen_Needleleaf_Forest;
				}
				else if (tree_typ == DB_TREES)
					vtype = Lynx_Subtropical_Deciduous_Broadleaf_Forest;
				else if (tree_typ == EB_TREES)
					vtype = Lynx_Subtropical_Evergreen_Broadleaf_Forest;
				else
					vtype = Lynx_Subtropical_Mixed_Forest;
			}
			else if (treec >= woodl_thres)
			{ if (tree_typ == EN_TREES)
				vtype = Lynx_Subtropical_Evergreen_Needleleaf_Woodland;
				else if (tree_typ == DB_TREES)
					vtype = Lynx_Subtropical_Deciduous_Broadleaf_Woodland;
				else if (tree_typ == EB_TREES)
					vtype = Lynx_Subtropical_Evergreen_Broadleaf_Woodland;
				else
					vtype = Lynx_Subtropical_Mixed_Woodland;
			}
			else if (treec >= shrubl_thres)
				vtype = Lynx_Subtropical_Shrubland;
			else
				vtype = Lynx_Subtropical_Desert;
			break;
			case TROPICALzone:
			if (grassc >= grass_thres && treec <= grass_thres)
				vtype = Lynx_Tropical_Grassland;
			else if (treec >= forest_thres)
				vtype = Lynx_Tropical_Evergreen_Broadleaf_Forest;
			else if (treec >= woodl_thres)
			{ if (ddecid <= ddecid_thres)
				vtype = Lynx_Tropical_Deciduous_Woodland;
				else
					vtype = Lynx_Tropical_Savanna;
			}
			else if (treec >= shrubl_thres)
				vtype = Lynx_Tropical_Shrubland;
			else
				vtype = Lynx_Tropical_Desert;
			break;
		}

	if (inVals.nlayer == 0 || npp <= 0.)
		vtype = Lynx_Barren;

	assert(vtype <= MAX_VTYPE);

	return vtype;
	}


	/*
	   int vemap_agg_vtype_from_vtypeVINCERA[] =
	   {0,  6,  1, 1, 1, 3, 3, 2, 4, 4,  5,  5,  5, 5, 5, 5,  5, 6, 6, 6, 6, 7, 5, 2};
	//0  1   2  3  4  5  6  7  8  9  10  11  12 13 14 15  16 17 18 19 20 21 22 23

	void CENTURY_BiogeochemModel::BiogeogLYNXandWWETAC()
	{

	...

	if (m_nlayer == 0. ||
	m_npp_yr <= 0.)
	vtype = 1; // Barren 

	assert(vtype<=MAX_VTYPE);

	// 
	VTYPE; PHYSIOGNOMIC_CLASS Key

	1. ice; 2. Barren
	2. tundra aka alpine; 10. Shrub
	3. taiga-tundra; 4. Conifer Woodland
	4. boreal evergreen needleleaf forest; 3. Conifer Forest
	5. boreal mixed woodland; 4. Conifer Woodland
	6. subalpine; 3. Conifer Forest
	7. maritime evergreen needleleaf forest; 3. Conifer Forest
	8. temperate evergreen needleleaf forest; 3. Conifer Forest
	9. temperate deciduous broadleaf forest; 7. Hardwood Forest
	10. temperate cool mixed forest; 7. Hardwood Forest
	11. temperate warm mixed forest; 7. Hardwood Forest
	12. temperate evergreen needleleaf woodland; 4. Conifer Woodland 
	13. temperate deciduous broadleaf woodland; 8. Hardwood Woodland
	14. temperate cool mixed woodland; 8. Hardwood Woodland
	15. temperate warm mixed woodland; 8. Hardwood Woodland
	16. temperate shrubland; 10. Shrub
	17. temperate grassland; 9. Herbaceous
	18. temperate desert; 5. Desert Shrub
	19. subtropical evergreen needleleaf forest; 3. Conifer Forest
	20. subtropical deciduous broadleaf forest; 7. Hardwood Forest
	21. subtropical evergreen broadleaf forest; 7. Hardwood Forest
	22. subtropical mixed forest; 7. Hardwood Forest
	23. subtropical evergreen needleleaf woodland; 4. Conifer Woodland
	24. subtropical deciduous broadleaf woodland; 8. Hardwood Woodland
	25. subtropical evergreen broadleaf woodland; 8. Hardwood Woodland
	26. subtropical mixed woodland; 8. Hardwood Woodland
	27. subtropical shrubland; 10. Shrub
	28. subtropical grassland; 9. Herbaceous
	29. subtropical desert; 5. Desert Shrub
	30. tropical evergreen broadleaf forest; 7. Hardwood Forest
	31. tropical deciduous woodland; 8. Hardwood Woodland
	32. tropical savanna; 8. Hardwood Woodland
	33. tropical shrubland; 10. Shrub
	34. tropical grassland; 9. Herbaceous
	35. tropical desert; 5. Desert Shrub
	36. cool needleleaf forest; 3. Conifer Forest
	37. coniferous xeromorphic woodland; 4. Conifer Woodland
	38. subalpine meadow; 9. Herbaceous
	39. water; 2. Barren
	40. natural barren; 2. Barren
	41. developed; 2. Barren


	PHYSIOGNOMIC_CLASS Key 

	1. (not used)
	2. Barren 
	3. Conifer Forest
	4. Conifer Woodland
	5. Desert Shrub
	6. Desert Woodland
	7. Hardwood Forest
	8. Hardwood Woodland
	9. Herbaceous
		10. Shrub



		{ int agg_vtype_from_vtypeLYNXandWWETAC[] =
			{0,  2, 10, 4, 3, 4, 3, 3, 3, 7,  7,  7,  4, 8, 8, 8, 10, 9, 5, 3, 7,
				//0  1   2  3  4  5  6  7  8  9  10  11  12 13 14 15  16 17 18 19 20
				7,  7, 4, 8, 8, 8, 10, 9, 5,  7,  8,  8, 10, 9, 5, 3, 4, 9, 2, 2, 2};
			//  21  22 23 24 25 26  27 28 29  30  31  32  33 34 35 36 37 38 39 40 41
			assert(0<=vtype && vtype<sizeof(agg_vtype_from_vtypeLYNXandWWETAC)/sizeof(int));
			outputData.intOutvars[OUTagg_vtype] = agg_vtype_from_vtypeLYNXandWWETAC[vtype];
		}

	outputData.intOutvars[OUTvtype] = vtype;
	outputData.intOutvars[OUTzone] = zone;
	outputData.intOutvars[OUTtree_typ] = tree_typ;
	outputData.intOutvars[OUTgrass_typ] = grass_typ;
	outputData.intOutvars[OUTddecid] = ddecid;

} // end of BiogeogLYNXandWWETAC()


int CENTURY_BiogeochemModel::vemap_aggregate_vegtype(unsigned int vtype, BaseCalibrationEnum biogeog)
	// Convert vegetation type to VEMAP aggregated vegetation type
	// Aggregated vegetation types
	// from Bachelet et al. Simulating past and future dynamics of natural ecosystems in the United States
	// Global Biogeochemical Cycles, 17(2):1045 (2003)
	// 1 coniferous forests
	// 2 winter deciduous forests
	// 3 mixed forests
	// 4 broadleaf and evergreen dought deciduous forests
	// 5 savannas and woodlands
	// 6 grasslands and shrublands
	// 7 deserts

	// LYNX, WWETAC, and US50km crosswalk to VEMAP aggregated vegetation types

	1. ice; 7. desert
	2. tundra aka alpine; 6. grasslands and shrublands
	3. taiga-tundra; 5. savannas and woodlands
	4. boreal evergreen needleleaf forest; 1. coniferous forests
	5. boreal mixed woodland; 5. savannas and woodlands
	6. subalpine; 1. coniferous forests
	7. maritime evergreen needleleaf forest; 1. coniferous forests
	8. temperate evergreen needleleaf forest; 1. coniferous forests
	9. temperate deciduous broadleaf forest; 2. winter deciduous forests
	10. temperate cool mixed forest; 3. mixed forests
	11. temperate warm mixed forest; 3. mixed forests
	12. temperate evergreen needleleaf woodland; 5. savannas and woodlands 
	13. temperate deciduous broadleaf woodland; 5. savannas and woodlands
	14. temperate cool mixed woodland; 5. savannas and woodlands
	15. temperate warm mixed woodland; 5. savannas and woodlands
	16. temperate shrubland; 6. grasslands and shrublands
	17. temperate grassland; 6. grasslands and shrublands
	18. temperate desert; 7. deserts
	19. subtropical evergreen needleleaf forest; 1. coniferous forests
	20. subtropical deciduous broadleaf forest; 4 broadleaf and evergreen dought deciduous forests
	21. subtropical evergreen broadleaf forest;  4 broadleaf and evergreen dought deciduous forests
	22. subtropical mixed forest; 3 mixed forests
	23. subtropical evergreen needleleaf woodland; 5 savannas and woodlands
	24. subtropical deciduous broadleaf woodland; 5 savannas and woodlands
	25. subtropical evergreen broadleaf woodland; 5 savannas and woodlands
	26. subtropical mixed woodland; 5 savannas and woodlands
	27. subtropical shrubland; 6 grasslands and shrublands
	28. subtropical grassland; 6 grasslands and shrublands
	29. subtropical desert; 7 deserts
	30. tropical evergreen broadleaf forest; 4 broadleaf and evergreen dought deciduous forests
	31. tropical deciduous woodland; 5 savannas and woodlands
	32. tropical savanna; 5 savannas and woodlands
	33. tropical shrubland; 6 grasslands and shrublands
	34. tropical grassland; 6 grasslands and shrublands
	35. tropical desert; 7 deserts
	36. cool needleleaf forest; 1. coniferous forests
	37. coniferous xeromorphic woodland; 5 savannas and woodlands
	38. subalpine meadow; 6 grasslands and shrublands
	39. water; 7. deserts
	40. natural barren; 7. deserts
	41. developed; 7. deserts

	// AGG_VTYPE crosswalk to VEMAP aggregated vegetation types
	AGG_VTYPE Key 

1. (not used)
	2. Barren; 7 desert
	3. Conifer Forest; 1 coniferous forests
	4. Conifer Woodland; 5 savannas and woodlands
	5. Desert Shrub; 7 desert
	6. Desert Woodland; 7 desert
	7. Hardwood Forest; 2 winter deciduous forests
	8. Hardwood Woodland; 5 savannas and woodlands
	9. Herbaceous; 6 grasslands and shrublands
	10. Shrub; 6 grasslands and shrublands

{
	int rtnval;
	unsigned int vemap_agg_vtype_from_agg_vtype[] =
	{0, 0, 7, 1, 5, 7, 7, 2, 5, 6, 6};
	//0  1  2  3  4  5  6  7  8  9 10 
	unsigned int agg_vtype;

	switch (biogeog)
	{
		case mc2NA8km:
		case mc2WWETAC:
		case mc2US50km:
			{ unsigned int vemap_agg_vtype_from_vtypeLYNXandWWETAC[] =
				{0,  7,  6, 6, 1, 5, 1, 1, 1, 2,  3,  3,  5, 5, 5, 5,  6, 6, 7, 1, 4,
					//0  1   2  3  4  5  6  7  8  9  10  11  12 13 14 15  16 17 18 19 20
					4,  3, 5, 5, 5, 5,  6, 6, 7,  4,  5,  5,  6, 6, 7, 1, 5, 6, 7, 7, 7};
				//  21  22 23 24 25 26  27 28 29  30  31  32  33 34 35 36 37 38 39 40 41
				assert(0<=vtype && vtype<sizeof(vemap_agg_vtype_from_vtypeLYNXandWWETAC)/sizeof(int));
				rtnval = vemap_agg_vtype_from_vtypeLYNXandWWETAC[vtype];
			}
			break;
		case mc2CA08:
			{ int agg_vtype_from_vtypeCA08[] =
				{0,  2,  3, 3, 3, 7, 7, 4, 8, 8, 10, 10, 10, 9, 5, 2, 2, 2, 9};
				//0  1   2  3  4  5  6  7  8  9  10  11  12 13 14 15 16 17 18
				assert(0<=vtype && vtype<sizeof(agg_vtype_from_vtypeCA08)/sizeof(int));
				agg_vtype = agg_vtype_from_vtypeCA08[vtype];
				assert(0<=agg_vtype && agg_vtype<sizeof(vemap_agg_vtype_from_agg_vtype)/sizeof(int));
				rtnval = vemap_agg_vtype_from_agg_vtype[agg_vtype];
			}
			break;
		case mc2YOSE:
			{ int agg_vtype_from_vtypeYOSE[] =
				{0, 2,  2, 4, 4, 5, 5, 4, 5, 5, 10, 10, 10, 13, 10, 15, 15, 15, 2};
				//0  1   2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18
				assert(0<=vtype && vtype<sizeof(agg_vtype_from_vtypeYOSE)/sizeof(int));
				agg_vtype = agg_vtype_from_vtypeYOSE[vtype];
				assert(0<=agg_vtype && agg_vtype<sizeof(vemap_agg_vtype_from_agg_vtype)/sizeof(int));
				rtnval = vemap_agg_vtype_from_agg_vtype[agg_vtype];
			}
			break;
		case mc2VINCERA:
			assert(0<=vtype && vtype<sizeof(vemap_agg_vtype_from_vtypeVINCERA)/sizeof(int));
			rtnval = vemap_agg_vtype_from_vtypeVINCERA[vtype];
			break;
		default: assert(0);
	}

	return(rtnval);
} // end of vemap_aggregate_vegtype()

*/
