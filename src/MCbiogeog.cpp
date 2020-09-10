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
#include "PNVmodel.h"


MC_BiogeographyModel::MC_BiogeographyModel(Simulation * sP, RunParamsClass * rP, ModelParamsClass * mP)
{ 
	ProcessModelInit(sP);
	MAPSS_model_parameters MAPSSparams(mP->MAPSSparameterSet, 1.0f);

	modelParamsP = mP; 
	runParamsP = rP;
	m_taiga_tundra_threshold = MAPSSparams.taiga_tundra_boundary[mapssBOREAL];

} // end of constructor for MC_BiogeographyModel


bool MC_BiogeographyModel::runModelAndCollectData(const int year_index, const int row_index, const int col_index) 
{
	return(true);
}; // end of MC_BiogeographyModel::runModelAndCollectData()


bool MC_BiogeographyModel::BiogeogMC2(BiogeographyInputData inVals, 
		BiomeType * biomeP, PhysiognomicClass * physiognomic_classP, MC2VegType * vtypeP)
{
	bool rtnFlag = true;
	ModelParamsClass * mP = modelParamsP;

	MC2VegType vtype;

	// C3C4Dominance grass_typ;

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
	float smoothed_mean_annual_ppt = inVals.ppt_yr;

	// grass_typ = c3>=c3_threshold ? C3Dominance : C3C4Mixed; /* CLASSIFY GRASS LIFEFORM */

	/* Biogeography rule set */

	// Determine biome: one of desert, shrubland, grassland/savanna, woodland, or forest.
	m_biome = UNKNOWNbiome;
	if (inVals.npp_yr<=0. || (treec<desert_treec_max && max_grassc<desert_grassc_max)) 
		m_biome = DESERTbiome;
	else if (treec<woodl_thres) m_biome = max_grass_frac<grassfrac_thres ? SHRUBLANDbiome : GRASSLAND_SAVANNAbiome;
	else if (treec<forest_thres) m_biome = WOODLANDbiome;
	else m_biome = FORESTbiome;

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
		case SHRUBLANDbiome: m_physiognomic_class = SHRUB_STEPPEpclass; break;
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
	if (m_physiognomic_class<=SAVANNApclass) switch (m_physiognomic_class)
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
		case SEMIDESERT_SHRUBLANDpclass: // shrubland without a significant grass component
			switch (zone)
		{
			case ARCTICzone: 
			case BOREALzone: 
				vtype = TUNDRAveg; 
				break;
			case TEMPERATEzone: 
			case SUBTROPICALzone:
				if (modelParamsP->code_flags[REGRESSION_TEST_2B101_FLAG])
					vtype = c3>=c3_threshold ? SHRUB_STEPPEveg : DRY_SHRUB_STEPPEveg;
				else vtype = SEMIDESERT_SHRUBLANDveg; 
				break;
			case TROPICALzone: 
				assert(c3<c3_threshold); 
				vtype = modelParamsP->code_flags[REGRESSION_TEST_2B101_FLAG] ? DRY_SHRUB_STEPPEveg : TROPICAL_SHRUBLANDveg; break;
			default: assert(0); break;
		}
			break;
		case SHRUB_STEPPEpclass:
			switch (zone)
		{
			case ARCTICzone: 
			case BOREALzone: 
				if (gdd_zero>=m_taiga_tundra_threshold) 
					vtype = modelParamsP->code_flags[REGRESSION_TEST_2B101_FLAG] ? SHRUB_STEPPEveg : BOREAL_SHRUBLANDveg;  
				else vtype = (inVals.mean_vegc>=taiga_tundra_C_min) ? TAIGA_TUNDRAveg : TUNDRAveg; 
				break;
			case TEMPERATEzone: 
			case SUBTROPICALzone: 
				if (modelParamsP->code_flags[REGRESSION_TEST_2B101_FLAG])
					vtype = c3>=c3_threshold ? SHRUB_STEPPEveg : DRY_SHRUB_STEPPEveg;
				else vtype = smoothed_mean_annual_ppt>=modelParamsP->shrub_steppe_precip_threshold ? 
					SHRUB_STEPPEveg : DRY_SHRUB_STEPPEveg; 
				break;
			case TROPICALzone: 
				vtype = modelParamsP->code_flags[REGRESSION_TEST_2B101_FLAG] ? DRY_SHRUB_STEPPEveg : TROPICAL_SHRUBLANDveg; 
				break;
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
			case SUBTROPICALzone: vtype = (c3>=c3_threshold || inVals.tot_summer_ppt<modelParamsP->c4grass_min_summer_precip) ?
C3GRASSveg : C4GRASSveg; break;
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
			case SUBTROPICALzone: vtype = (c3>=c3_threshold || inVals.tot_summer_ppt<modelParamsP->c4grass_min_summer_precip) ? 
C3GRASSveg : C4GRASSveg; break;
			case TROPICALzone: vtype = TROPICAL_SAVANNAveg; break;
			default: assert(0); break;
		}
			break;
		default: assert(0); break;
	} // end of if (m_physiognomic_class<=SAVANNApclass) switch (m_physiognomic_class)
	else if (runParamsP->baseCalibration==mc2W_WA) 
	{
		if (zone<TROPICALzone)
		{
			MC2VegType vtypeFromPNVcode[] = {
				// vtype, PNVcode, PNV vegzones
				UNKNOWNveg, // 0
				UNKNOWNveg, // 1
				UNKNOWNveg, // 2
				SEMIDESERT_SHRUBLANDveg, // 3,SDZ,Salt Desert Zone
				C3GRASSveg, // 4,GRZ,Grassland Zone
				SHRUB_STEPPEveg, // 5,STZ,Steppe Zone
				TEMPERATE_EN_WOODLANDveg, // 6,WJZ,Western Juniper Zone
				UNKNOWNveg, // 7
				LPPZveg, // 8,LPPZ,Lodgepole Pine Zone
				SSZveg, // 9,SSZ,Sitka Spruce Zone
				PPZveg, // 10,PPZ,Ponderosa Pine Zone
				TEMPERATE_WARM_MIXED_FORESTveg, // 11,OWOZ,Oregon White Oak Zone
				JPZveg, // 12,JPZ,Jeffrey Pine Zone
				POCZveg, // 13,POCZ,Port Orford-Cedar Zone
				TEMPERATE_NEEDLELEAF_FORESTveg, // 14,DFZ,Douglas-fir Zone
				SUBTROPICAL_MIXED_FORESTveg, // 15,TOZ,Tan Oak Zone
				GFZveg, // 16,GFZ,Grand Fir Zone
				WWPZveg, // 17,WWPZ,Western White Pine Zone
				DFZ2veg, // 18,DFZ2,Douglas-fir 2 Zone
				WHZveg, // 19,WHZ,Western Hemlock Zone
				WFZveg, // 20,WFZ,White Fir Zone
				SRFZveg, // 21,SRFZ,Shasta Red Fir Zone
				PSFZveg, // 22,SFZ,Pacific Silver Fir Zone
				MHZveg, // 23,MHZ,Mountain Hemlock Zone
				UNKNOWNveg, // 24
				SAFZveg, // 25,SAFZ,Subalpine Fir Zone
				UNKNOWNveg, // 26
				UNKNOWNveg, // 27
				UNKNOWNveg, // 28
				UNKNOWNveg, // 29
				UNKNOWNveg, // 30
				UNKNOWNveg, // 31
				PKLZveg, // 32,PKLZ,Parkland Zone
				TUNDRAveg}; // 33,ALPZ,Alpine Zone  

			PNVmodel PNVinstance(pS, &inVals);
			int pnv_code_as_int = (int)PNVinstance.pnv(&inVals, mP);
			// printf("*** BiogeogMC2(): pnv_code = %d\n", pnv_code_as_int);
			vtype = (pnv_code_as_int>=0 && pnv_code_as_int<(int)(sizeof(vtypeFromPNVcode)/sizeof(int))) ? 
				vtypeFromPNVcode[pnv_code_as_int] : UNKNOWNveg;
		} 
		else vtype = UNKNOWNveg; // PNVmodel isn't calibrated for TROPICALzone
	} // end of if (runParamsP->baseCalibration==mc2W_WA) when pclass>SAVANNApclass
	else switch (m_physiognomic_class)
	{
		case EVERG_WOODLANDpclass:
			switch (zone)
		{
			case ARCTICzone: 
			case BOREALzone: vtype = gdd_zero>=subalp_thres ? TEMPERATE_EN_WOODLANDveg : BOREAL_WOODLANDveg; break;
			case TEMPERATEzone: vtype = tree_typ==EN_TREES ? TEMPERATE_EN_WOODLANDveg : TEMPERATE_WARM_MIXED_WOODLANDveg; break;
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
						else if (modelParamsP->code_flags[REGRESSION_TEST_2B89_FLAG]) vtype = TEMPERATE_NEEDLELEAF_FORESTveg;
						else if (inVals.ppt_yr>mP->moist_temperate_threshold) vtype = MOIST_TEMPERATE_NEEDLELEAF_FORESTveg;
						else if (inVals.ppt_yr>mP->dry_temperate_threshold) vtype = TEMPERATE_NEEDLELEAF_FORESTveg;
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
						if (gdd_zero<=subalp_thres) vtype = SUBALPINE_FORESTveg;

						else if (modelParamsP->code_flags[REGRESSION_TEST_2B89_FLAG])
						{
							if (cont_index<=mari_thres) 
							{              
								if (inVals.min_smoothed_tmp<modelParamsP->tmmin_threshold)
								{
									if (inVals.ppt_yr>mP->moist_temperate_threshold) vtype = MOIST_TEMPERATE_NEEDLELEAF_FORESTveg;
									else if (inVals.ppt_yr>mP->dry_temperate_threshold) vtype = TEMPERATE_NEEDLELEAF_FORESTveg;
									else vtype = DRY_TEMPERATE_NEEDLELEAF_FORESTveg;
								}
								else vtype = MARITIME_EN_FORESTveg;
							}
							else vtype = TEMPERATE_NEEDLELEAF_FORESTveg; 
						} // end of else if (code_flags[REGRESSION_TEST_2B89_FLAG])

						else if (cont_index<=mari_thres && inVals.min_smoothed_tmp>=modelParamsP->tmmin_threshold) 
							vtype = MARITIME_EN_FORESTveg;
						else if (inVals.ppt_yr>mP->moist_temperate_threshold) vtype = MOIST_TEMPERATE_NEEDLELEAF_FORESTveg;
						else if (inVals.ppt_yr>mP->dry_temperate_threshold) vtype = TEMPERATE_NEEDLELEAF_FORESTveg;
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
	} // end of switch (m_physiognomic_class) for woodland and forest classes when base calibration is not W_WA

	*biomeP = m_biome;
	*physiognomic_classP = m_physiognomic_class;
	*vtypeP = vtype;

	rtnFlag = vtype>UNKNOWNveg and vtype<=MAX_VTYPE;
	return(rtnFlag); 

} // end of BiogeogMC2()


int MC_BiogeographyModel::GrouseModel(int vtype, float smrpre, float augmaxt, float anntmp)
	// Returns 0 for no grouse, otherwise returns vtype.
{
	int grouse_habitat_type;
	bool grouse_absence_flag; 

	switch (vtype)
	{
		case SHRUB_STEPPEveg: 
		case C3GRASSveg: 
		case DRY_SHRUB_STEPPEveg: 
		case C4GRASSveg: 
			grouse_absence_flag = smrpre <= modelParamsP->grouse_smrpre_threshold
				&& augmaxt > modelParamsP->grouse_augmaxt_threshold
				&& anntmp > modelParamsP->grouse_anntmp_threshold;
			grouse_habitat_type = grouse_absence_flag ? 0 : vtype;
			break;

		default: 
			grouse_habitat_type = 0;
			break;
	}

	return(grouse_habitat_type);
} // end of GrouseModel()      


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
	} // end of BiogeogLC()


