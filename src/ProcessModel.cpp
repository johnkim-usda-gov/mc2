/*
 *  ProcessModel.cpp
 *  mc2
 *
 *  Created by Dave Conklin on 12/30/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>	// JK
#include <stdlib.h>	// JK
#include <math.h>
#include "assert.h"

#ifndef ASSERT
#define ASSERT assert
#endif

#include "netcdf.h"
#include "ScienceFcns.h"
#include "category_bgc.h"
#include "MAPSSvegClasses.h"

#include "ProcessModel.h"

#include "MAPSSbiogeographyModel.h"
#include "CENTURY.h"
#include "MCfire.h"
#include "MC2.h"

using namespace std;


void ProcessModel::ProcessModelInit(Simulation * pRun)
{
	pS = pRun;
	runParamsP = &(pRun->runParams);
	modelParamsP = &(pRun->modelParams);

	m_portable_next = 1; // random number seed

	// constants
	c_frost = 13.0; // deg C
	c_s_decid_bound = 2.25; // deg C
	c_n_decid_bound = -16.0; // deg C
	c_az_thres = 1000.; // upper GDD limit for arctic aka alpine zone, growing degree-days ref to 0 deg C
	// c_bz_thres = -13.0; // upper min_smoothed_tmp limit for boreal zone, deg C
	c_tz_thres = 7.75; // upper min_smoothed_tmp limit for temperate zone, deg C
	c_stz_thres = 18.; // upper min_smoothed tmp limit for subtropical zone, deg C
	c_z_all = 10.f; // z for the Penman-Monteith calculations for all zones and lifeforms

	// fireState will be initialized to reasonable values when a fire model object is instantiated.
	stateData.fireState.sum_ann_ppt = NC_FILL_FLOAT;
	stateData.fireState.Pprev = NC_FILL_FLOAT;
	stateData.fireState.Kprev = NC_FILL_FLOAT;
	stateData.fireState.prev_l1hr = NC_FILL_FLOAT;
	stateData.fireState.prev_dstand = NC_FILL_FLOAT;
	stateData.fireState.prev_litter = NC_FILL_FLOAT;
	stateData.fireState.prev_1hr = NC_FILL_FLOAT;
	stateData.fireState.prev_10hr = NC_FILL_FLOAT;
	stateData.fireState.prev_100hr = NC_FILL_FLOAT;
	stateData.fireState.prev_1000hr = NC_FILL_FLOAT;
	stateData.fireState.prev_mxd = NC_FILL_FLOAT;
	stateData.fireState.prev_mc_grass = NC_FILL_FLOAT;
	stateData.fireState.prev_mc_tree = NC_FILL_FLOAT;
	stateData.fireState.prev_depth = NC_FILL_FLOAT;
	stateData.fireState.clai_in_midseason = NC_FILL_FLOAT;
	stateData.fireState.prev_snow = NC_FILL_FLOAT;
	stateData.fireState.F1 = NC_FILL_INT;
	stateData.fireState.F2 = NC_FILL_INT;
	stateData.fireState.rand_seed = 0;
	stateData.fireState.yrs_in_sum_ann_ppt = NC_FILL_INT;
	stateData.fireState.ffmc_prev = NC_FILL_FLOAT;
	stateData.fireState.dmc_prev = NC_FILL_FLOAT;
	stateData.fireState.dc_prev = NC_FILL_FLOAT;

} // end of ProcessModelInit()


InputDataClass::InputDataClass() 
{ 
	pptP = tminP = tmaxP = tmpP = vprP = tdmeanP = wndP = inFromMAPSS = inFromEQ = inFromWS = NULL; 
	elev = fog = topomoist = deltaTsl = aspect = sw = cad = NC_FILL_FLOAT;
	lat = lon = northing = easting = NC_FILL_DOUBLE;
	num_months = NC_FILL_INT;
} // end of default constructor for InputDataClass


InputDataClass::~InputDataClass()
{
	if (wndP!=NULL) free(wndP);  
	if (vprP!=NULL) free(vprP);
	if (tdmeanP!=NULL) free(tdmeanP);
	if (tmpP!=NULL) free(tmpP);
	if (tmaxP!=NULL) free(tmaxP);
	if (tminP!=NULL) free(tminP);
	if (pptP!=NULL) free(pptP);
	if (inFromWS!=NULL) free(inFromWS);
	if (inFromEQ!=NULL) free(inFromEQ);
	if (inFromMAPSS!=NULL) free(inFromMAPSS);
} // end of ~InputDataClass()


// Portable random number generation routines
// Adapted from p.337 of P.J. Plauger's book, "The Standard C Library", 1992.
// Used here instead of the rand() and srand() functions from <stdlib.h> because
// rand() and srand() may generate a different sequence of random numbers on different computers.
// Note that portable_rand() is set up to generate a series which is uniformly distributed on the 
// interval [0,32767], regardless of what the local value of RAND_MAX is.  However, the local 
// value of RAND_MAX is constrained by the C standard to be at least 32767.

int ProcessModel::portable_rand(void)
{
	m_portable_next = m_portable_next*1103515245 + 12345;
	return (unsigned int)(m_portable_next/65536)%32768;
} // end of ProcessModel::portable_rand()

void ProcessModel::portable_srand(unsigned int seed)
{
	assert(seed!=0); // seeding with 0 leads to portable_rand() always returning 0
	m_portable_next = seed;
} // end of ProcessModel::portable_srand()


VCategory ProcessModel::convertMAPSStoVEMAP(MAPSSvegClass mclass)
{
	VCategory		vclass = VUnknown;

	switch (mclass) 
	{		
		case DesertBoreal :                             /*500*/
		case Tundra :					/*601*/
			vclass = VTundra;				/*1 NWT ARC*/
			break;
		case TaigaTundra :				/*600*/
			vclass = VTaiga;
			break;
		case ForestDeciduousNeedle :			/*103*/
			// case TreeSavannaDeciduousNeedle:		/*202*/
			vclass = VBorealLarchForest;
			break;
		case ForestEvergreenNeedleTaiga :               /*107*/
		case ForestSubalpine: /* 114 */
		case TreeSavannaSubalpine: /* 214 */
			vclass = VBorealConiferousForest;		/*2 BNZ CPR3*/
			break;
		case ForestMixedWarm_EN :                       /*108*/
		case ForestEvergreenNeedleMaritime :		/*112*/
		case TreeSavannaEvergreenNeedleMaritime :       /*207*/
			vclass = VMaritimeTemperateConiferousForest;	/*3 AND CPR3*/
			break;
		case ForestEvergreenNeedleContinental :		/*113*/
		case TreeSavannaEvergreenNeedleContinental :    /*208*/
			vclass = VContinentalTemperateConiferousForest;	/*4 WPINE CPR3*/
			break;
		case ForestMixedCool :				/*102*/
			vclass = VCoolTemperateMixedForest;		/*5 CWT CPR3*/
			break;
		case ForestMixedWarm_DEB :			/*101*/
		case TreeSavannaMixedWarm_EN :                  /*206*/
		case ForestSeasonalTropical_ED :                /*109*/
		case ForestSavannaDryTropical_ED :              /*110*/
			vclass = VWarmTemperateSubtropicalMixedForest;	/*6 CWT KNZ*/
			break;
		case ForestDeciduousBroadleaf :			/*100*/
		case ForestHardwoodCool :			/*111*/
			vclass = VTemperateDeciduousForest;		/*7 CWT CPR*/
			break;
		case ForestEvergreenBroadleafTropical :		/*105*/
			vclass = VTropicalEvergreenForest;		/*9 LUQ TKNZ*/
			break;
		case TreeSavannaPJMaritime :                    /*210*/
		case ShrublandSubTropicalXeromorphic :		/*311*/
			vclass = VTemperateMixedXeromorphicWoodland;	/*10 THORN JRN*/
			break;
		case TreeSavannaPJContinental :                 /*209*/
			vclass = VTemperateConiferXeromorphicWoodland;	/*11 WPINE JRN*/
			break;
		case ShrubSavannaTropical_EB :			/*305*/
			vclass = VTropicalThornWoodland;		/*12 THORN TKNZ*/
			break;
		case TreeSavannaDeciduousBroadleaf :		/*200*/
		case TreeSavannaMixedWarm_DEB :			/*201*/
		case TreeSavannaMixedCool_EN :                  /*205*/
		case GrassTallC3C4 :                            /*417*/
		case GrassTallC4 :                              /*420*/
			vclass = VTemperateSubtropicalDeciduousSavanna;	/*13 CWT KNZ*/
			break;
		case ShrubSavannaMixedWarm_DEB :		/*303*/
		case ShrubSavannaSubTropicalMixed :		/*310*/
			vclass = VWarmTemperateSubtropicalMixedSavanna;	/*14 CWT KNZ*/
			break;
		case TreeSavannaPJXericContinental :            /*211*/
			vclass = VTemperateConiferSavanna;		/*15 WPINE G3*/
			break;
		case GrassTallC3 :				/*414*/
		case GrassMidC3 :				/*415*/
		case GrassShortC3 :				/*416*/
		case GrassMidC3C4 :                             /*418*/
		case GrassShortC3C4 :                           /*419*/
		case GrassSemiDesertC3 :			/*423*/
		case GrassSemiDesertC3C4 :                      /*424*/
			vclass = VC3Grasslands;				/*17 CWT CPR3*/
			break;
		case GrassMidC4 :                               /*421*/
		case GrassShortC4 :                             /*422*/
			vclass = VC4Grasslands;				/*18 CWT KNZ*/
			break;
		case ShrublandSubTropicalMediterranean :	/*312*/
			vclass = VMediterraneanShrubland;		/*19 KNZ JR*/
			break;
		case OpenShrublandNoGrass :                     /*301*/ 
		case ShrubSavannaDeciduousBroadleaf :           /*302*/
			// case ShrubSavannaEvergreenNeedle :              /*306*/
		case ShrubSavannaMixedCool_EN :                 /*307*/ 
		case ShrubSavannaEvergreenMicro :               /*308*/
		case DesertTemperate :                          /*501*/
			vclass = VTemperateAridShrubland;		/*20 SAGE CPR3*/
			break;

		case DesertSubtropical :			/*502*/
		case DesertTropical :				/*503*/
		case DesertExtreme :				/*504*/
		case GrassSemiDesertC4 :                        /*425*/
			vclass = VSubtropicalAridShrubland;		/*21 THODR JORN*/
			break;
		case Ice :					/*602*/
			vclass = VIce;				/*90*/
			break;
			/*
			   case Wetlands :					// 700
			   vclass = VWetlands;			// 92
			   break;
			 */ 
		case Unknown :
		default :
			// assert(0); cells with bad soil data may come thru here
			break;
	}
	assert(VUnknown<=vclass && vclass<=VWetlands);
	return(vclass);
} // end of convertMAPSStoVEMAP()


float ProcessModel::calcPET(float tmp, float vpr, float vp_sat, float wnd, float elevation, float z, float z0, int days_in_month)
{
	float t, pet;
	t = (vp_sat>vpr) ? (float)sciFn.trbxfr(z, z0, tmp, tmp, vpr, vp_sat, wnd, elevation) : 0.; 

	// trbxfr() returns a positive number when an error occurs.
	if (t>=0.0 || vp_sat<vpr) 
	{
#ifdef UNIX_MAPSS
		// if ((vpr - vp_sat)>50.0f) printf("*** ProcessModel.cpp/calcPET(): vp_sat, vpr t = %f, %f, %f, %f\n", vp_sat, vpr, vpr - vp_sat, t);
#endif
		t = 0.;
		// ASSERT(0);
	}

	pet = t * SECONDS_PER_DAY * days_in_month * -1.0f;

	ASSERT(pet>=0.0);

	return(pet);
} // end of ProcessModel::calcPET()


/*
   void StateVariables::initializeSpinupData()
   {


   } // end of StateVariables::initializeSpinupData()
 */

void ProcessModel::treeType(float everg, float needle, BaseCalibrationEnum baseCalib, 
		TreeType * tree_typP, bool * deciduousP, bool * broadleafP)
// Use evergreen and needleleaf indices to determine tree type.
{
	assert(0.<=everg && everg<=100. && 0.<=needle && needle<=100.);

	*tree_typP = UNKNOWNtree_typ;
	switch (baseCalib)
	{
		default:
		case mc2WCR:
		case mc2GLOBAL:
		case mc2ConUS_LC:
		case mc2ConUS:
		case mc2California:
			if (everg > 90. && needle > 75.) *tree_typP = EN_TREES;
			else if (everg > 66. && needle <= 33.) *tree_typP = EB_TREES;
			else if (everg <= 33. && needle > 66.) *tree_typP = DN_TREES;
			else if (everg <= 33. && needle <= 33.) *tree_typP = DB_TREES;
			else if (everg > 66. && needle > 33. && needle <= 66.) *tree_typP = EN_EB_TREES;
			else if (everg > 33. && everg <= 66. && needle > 75.) *tree_typP = DN_EN_TREES;
			else if (everg > 33. && everg <= 66. && needle <= 33.) *tree_typP = DB_EB_TREES;
			else if ((everg <=66. && needle > 33. && needle <= 66.) 
					|| (everg > 33. && needle > 66. && needle <= 75.)
					|| (everg > 66. && everg <= 90. && needle > 75.)) *tree_typP = EN_DB_TREES; 
			break;
			/*      
				case mc2VEMAP:
				case mc2NA8km:
				case mc2CA08:
				case mc2YOSE:
				case mc2WIES_YOSE:
				case mc2WWETAC:
				case mc2VINCERA:
				case mc2US10kmAlbers:
				if (everg > 66. && needle > 66.) *tree_typP = EN_TREES;         
				else if (everg > 66. && needle < 33.) *tree_typP = EB_TREES;         
				else if (everg < 33. && needle > 66.) *tree_typP = DN_TREES;        
				else if (everg < 33. && needle < 33.)  *tree_typP = DB_TREES;
				else if (everg > 66. && needle >= 33. && needle <= 66.) *tree_typP = EN_EB_TREES;         
				else if (everg >= 33. && everg <= 66. && needle > 66.) *tree_typP = DN_EN_TREES;         
				else if (everg >= 33. && everg <= 66. && needle < 33.) *tree_typP = DB_EB_TREES;         
				else *tree_typP = EN_DB_TREES;
				break;
			 */
	}

	assert(*tree_typP!=UNKNOWNtree_typ);

	switch (*tree_typP)
	{
		case EN_TREES:
			*deciduousP = false;
			*broadleafP = false;
			break;

		case EN_DB_TREES:
		case DB_TREES:
		case DB_EB_TREES:
			*deciduousP = true;
			*broadleafP = true;
			break;

		case EN_EB_TREES:
		case EB_TREES:
			*deciduousP = false;
			*broadleafP = true;
			break;

		case DN_TREES:
		case DN_EN_TREES:
			*deciduousP = true;
			*broadleafP = false;
			break;

		default:
			assert(0);
			break;
	}

} // end of ProcessModel::treeType()


int ProcessModel::convertVEMAPtoVtype(VCategory vclass, BaseCalibrationEnum baseCalib)
{
	int vtype;
	int vtypeGLOBAL[24] = {
		0, // 0 VUnknown -> unknown
		2, // 1 VTundra -> 2 tundra aka alpine
		4, // 2 VBorealConiferousForest -> 4 boreal needleleaf forest
		7, // 3 VMaritimeTemperateConiferousForest -> 7 maritime needleleaf forest
		8, // 4 VContinentalTemperateConiferousForest -> 8 temperate needleleaf forest
		10, // 5 VCoolTemperateMixedForest -> 10 cool mixed forest
		11, // 6 VWarmTemperateSubtropicalMixedForest -> 11 temperate warm mixed forest
		9, // 7 VTemperateDeciduousForest -> 9 temperate deciduous broadleaf forest
		31, // 8 VTropicalDeciduousForest -> 31 tropical deciduous woodland
		30, // 9 VTropicalEvergreenForest -> 30 tropical evergreen broadleaf forest
		15, // 10 VTemperateMixedXeromorphicWoodland -> 15 temperate warm mixed woodland
		16, // 11 VTemperateConiferXeromorphicWoodland -> 16 C3 shrubland
		31, // 12 VTropicalThornWoodland -> 31 tropical deciduous woodland
		13, // 13 VTemperateSubtropicalDeciduousSavanna -> 13 temperate deciduous broadleaf woodland
		15, // 14 VWarmTemperateSubtropicalMixedSavanna -> 15 temperate warm mixed woodland
		12, // 15 VTemperateConiferSavanna -> 12 temperate needleleaf woodland
		31, // 16 VTropicalDeciduousSavanna -> 31 tropical deciduous woodland
		17, // 17 VC3Grasslands -> 17 C3 grassland
		28, // 18 VC4Grasslands -> 28 C4 grassland 
		27, // 19 VMediterraneanShrubland -> 27 C4 shrubland
		18, // 20 VTemperateAridShrubland -> 18 temperate desert and semidesert
		29, // 21 VSubtropicalAridShrubland -> 29 subtropical desert and semidesert
		5, // 22 VTaiga -> 5 boreal woodland
		42}; // 23 VBorealLarchForest -> 42 larch forest

	if (!(vclass>=0 && (unsigned int)vclass<=(int)sizeof(vtypeGLOBAL)/sizeof(int))
			&& vclass!=90)
	{
		printf("*** convertVEMAPtoVtype(): unrecognized vclass = %d\n", vclass);
		assert(0);
	}

	vtype = vclass==90 ? (int)COLD_BARRENveg : vtypeGLOBAL[(int)vclass]; 

	switch (baseCalib)
	{
		case mc2California:
		case mc2ConUS:
		case mc2GLOBAL:
			break;
		case mc2WCR:
			if (vclass==VMaritimeTemperateConiferousForest) vtype = WHZveg;
			else if (vclass==VTaiga) vtype = PKLZveg;
			break;
		case mc2ConUS_LC:
			if (vclass==VTemperateAridShrubland) vtype = Lynx_Temperate_Shrubland;
			else if (vclass==VBorealLarchForest) vtype = Lynx_Boreal_Evergreen_Needleleaf_Forest;
			break;
		default: assert(false); break;
	}

	return(vtype);
} // end of ProcessModel::convertVEMAPtoVtype();


void StateVariables::copyFStoCENoutvars()
{
	eqState.cenOutvars[FSsum_ann_ppt] = fireState.sum_ann_ppt;
	eqState.cenOutvars[FS_Pprev] = fireState.Pprev;
	eqState.cenOutvars[FS_Kprev] = fireState.Kprev;
	eqState.cenOutvars[FS_Ksum] = fireState.Ksum;
	eqState.cenOutvars[FSprev_l1hr] = fireState.prev_l1hr;
	eqState.cenOutvars[FSprev_dstand] = fireState.prev_dstand;
	eqState.cenOutvars[FSprev_litter] = fireState.prev_litter;
	eqState.cenOutvars[FSprev_1hr] = fireState.prev_1hr;
	eqState.cenOutvars[FSprev_10hr] = fireState.prev_10hr;
	eqState.cenOutvars[FSprev_100hr] = fireState.prev_100hr;
	eqState.cenOutvars[FSprev_1000hr] = fireState.prev_1000hr;
	eqState.cenOutvars[FSprev_day_mc_100hr] = fireState.prev_day_mc_100hr;
	eqState.cenOutvars[FSprev_day_mc_1000hr] = fireState.prev_day_mc_1000hr;
	eqState.cenOutvars[FSprev_mxd] = fireState.prev_mxd;
	eqState.cenOutvars[FSprev_mc_grass] = fireState.prev_mc_grass;
	eqState.cenOutvars[FSprev_mc_tree] = fireState.prev_mc_tree;
	eqState.cenOutvars[FSprev_depth] = fireState.prev_depth;
	eqState.cenOutvars[FSclai_in_midseason] = fireState.clai_in_midseason;
	eqState.cenOutvars[FSprev_snow] = fireState.prev_snow;
	eqState.cenOutvars[FS_F1] = fireState.F1;
	eqState.cenOutvars[FS_F2] = fireState.F2;
	eqState.cenOutvars[FSrand_seed] = fireState.rand_seed;
	eqState.cenOutvars[FSyrs_in_sum_ann_ppt] = fireState.yrs_in_sum_ann_ppt;
	eqState.cenOutvars[FSffmc_prev] = fireState.ffmc_prev;
	eqState.cenOutvars[FSdmc_prev] = fireState.dmc_prev;
	eqState.cenOutvars[FSdc_prev] = fireState.dc_prev;
	eqState.cenOutvars[FSyrs_since_fire] = fireState.yrs_since_fire;
} // end of StateVariables::copyFStoCENoutvars


