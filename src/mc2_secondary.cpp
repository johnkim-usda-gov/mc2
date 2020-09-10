/*
 *  mc2_secondary.cpp
 */

#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <cstring> 

#include "netcdf.h"
#include "ScienceFcns.h"
#include "category_bgc.h"
#include "MAPSSvegClasses.h"
#include "ProcessModel.h"
#include "MCfire.h"
#include "MAPSSbiogeographyModel.h"
#include "MC2.h"
#include "assert.h"

using namespace std;

OutVarDescription::OutVarDescription()
{
	name = "spareOutVar";
	descrip = "unused output variable";
	type = NC_FLOAT;
	scale = 1.0;
	interval = SINGLE_INTERVAL;
	show = false;
	categoricalFlag = false;
	max_categories = 0;
	units = "none";
	var_id = -1;
	moOutNdx = yrOutNdx = multiyrOutNdx = singleOutNdx = -1;
	dictNdx = -1;
} // end of constructor for OutVarDescription


void Simulation::initializeVarDict()
{
	for (int i = 0; i<NUM_OPTIONAL_OUTVARS; i++) 
	{
		VarDict[i].interval = i<=LAST_MONTHLY_VAR ? MONTH_INTERVAL : 
			(i<=LAST_YEARLY_VAR ? YEAR_INTERVAL : 
			 (i<=LAST_MULTIYR_VAR ? MULTIYR_INTERVAL : SINGLE_INTERVAL));
		VarDict[i].dictNdx = i;
	}

	VarDict[PPT].name     = "PPT";
	VarDict[PPT].descrip  = "unsmoothed precipitation";
	VarDict[PPT].units    = "mmH2O";

	VarDict[TMP].name     = "TMP";
	VarDict[TMP].descrip  = "mean unsmoothed temperature";
	VarDict[TMP].units    = "deg C";

	VarDict[TMAX].name     = "TMAX";
	VarDict[TMAX].descrip  = "unsmoothed mean value of diurnal max temperature";
	VarDict[TMAX].units    = "deg C";

	VarDict[TMIN].name     = "TMIN";
	VarDict[TMIN].descrip  = "unsmoothed mean value of diurnal min temperature";
	VarDict[TMIN].units    = "deg C";

	VarDict[AET].name     = "AET";
	VarDict[AET].descrip  = "actual evapotranspiration";
	VarDict[AET].units    = "cm H2O";

	VarDict[VPR].name     = "VPR";
	VarDict[VPR].descrip  = "vapor pressure";
	VarDict[VPR].units    = "Pa";

	VarDict[NPP].name     = "NPP";
	VarDict[NPP].descrip  = "net primary production";
	VarDict[NPP].units    = "g C m-2";

	VarDict[VTYPE].name     = "VTYPE";
	VarDict[VTYPE].descrip  = "potential vegetation type" ;
	VarDict[VTYPE].units    = "veg type" ;
	VarDict[VTYPE].categoricalFlag = true;
	VarDict[VTYPE].max_categories = MAX_VTYPE;

	VarDict[PART_BURN].name     = "PART_BURN";
	VarDict[PART_BURN].descrip  = "fraction of cell affected by wildfire" ;
	VarDict[PART_BURN].units    = "fraction" ;

	VarDict[BIOME].name     = "BIOME";
	VarDict[BIOME].descrip  = "biome" ;
	VarDict[BIOME].units    = "biome category" ;
	VarDict[BIOME].categoricalFlag = true;
	VarDict[BIOME].max_categories = MAX_BIOME;

	VarDict[PHYSIOGNOMIC_CLASS].name     = "PHYSIOGNOMIC_CLASS";
	VarDict[PHYSIOGNOMIC_CLASS].descrip  = "physiognomic class" ;
	VarDict[PHYSIOGNOMIC_CLASS].units    = "class" ;
	VarDict[PHYSIOGNOMIC_CLASS].categoricalFlag = true;
	VarDict[PHYSIOGNOMIC_CLASS].max_categories = MAX_PCLASS;

	VarDict[CONSUMED].name     = "CONSUMED";
	VarDict[CONSUMED].descrip  = "C in biomass consumed by fire" ;
	VarDict[CONSUMED].units    = "g C m-2" ;

	VarDict[C_VEG].name     = "C_VEG";
	VarDict[C_VEG].descrip  = "live biomass carbon" ;
	VarDict[C_VEG].units    = "g C m-2" ;

	VarDict[C_FOREST].name     = "C_FOREST";
	VarDict[C_FOREST].descrip  = "live tree carbon" ;
	VarDict[C_FOREST].units    = "g C m-2" ;

	VarDict[C_MAX_LIVE_GRASS_ABOVEGR].name     = "C_MAX_LIVE_GRASS_ABOVEGR";
	VarDict[C_MAX_LIVE_GRASS_ABOVEGR].descrip  = "max aboveground live grass carbon" ;
	VarDict[C_MAX_LIVE_GRASS_ABOVEGR].units    = "g C m-2" ;

	VarDict[FIRE_KILLED].name     = "FIRE_KILLED";
	VarDict[FIRE_KILLED].descrip  = "carbon in biomass killed by fire" ;
	VarDict[FIRE_KILLED].units    = "g C m-2" ;

	VarDict[GFRAC].name     = "GFRAC";
	VarDict[GFRAC].descrip  = "max grass fraction of live carbon" ;
	VarDict[GFRAC].units    = "fraction" ;

	VarDict[NEP].name     = "NEP";
	VarDict[NEP].descrip  = "net ecosystem production" ;
	VarDict[NEP].units    = "g C m-2" ;

	VarDict[NBP].name     = "NBP";
	VarDict[NBP].descrip  = "net biome production" ;
	VarDict[NBP].units    = "g C m-2" ;

	VarDict[RSP].name     = "RSP";
	VarDict[RSP].descrip  = "heterotrophic respiration" ;
	VarDict[RSP].units    = "g C m-2" ;

	VarDict[BIO_CONSUME_CENTURY].name     = "BIO_CONSUME_CENTURY";
	VarDict[BIO_CONSUME_CENTURY].descrip  = "carbon in biomass consumed by fire, per CENTURY model" ;
	VarDict[BIO_CONSUME_CENTURY].units    = "g C m-2" ;

	VarDict[CONSUMED_LIVE].name     = "CONSUMED_LIVE";
	VarDict[CONSUMED_LIVE].descrip  = "carbon in live biomass consumed by fire, per MCfire model" ;
	VarDict[CONSUMED_LIVE].units    = "g C m-2" ;

	VarDict[CONSUMED_DEAD].name     = "CONSUMED_DEAD";
	VarDict[CONSUMED_DEAD].descrip  = "carbon in dead biomass consumed by fire, per MCfire model" ;
	VarDict[CONSUMED_DEAD].units    = "g C m-2" ;

	VarDict[FFMC_ANN_MAX].name     = "FFMC_ANN_MAX";
	VarDict[FFMC_ANN_MAX].descrip  = "max fine fuel moisture content index (Canadian)" ;
	VarDict[FFMC_ANN_MAX].units    = "none" ;

	VarDict[BUI_ANN_MAX].name     = "BUI_ANN_MAX";
	VarDict[BUI_ANN_MAX].descrip  = "max fuel build up index (Canadian)" ;
	VarDict[BUI_ANN_MAX].units    = "none" ;

	VarDict[NPP_TREE].name     = "NPP_TREE";
	VarDict[NPP_TREE].descrip  = "tree net primary production" ;
	VarDict[NPP_TREE].units    = "g C m-2" ;

	VarDict[NPP_GRASS].name     = "NPP_GRASS";
	VarDict[NPP_GRASS].descrip  = "grass net primary production" ;
	VarDict[NPP_GRASS].units    = "g C m-2" ;

	VarDict[FPRD_PPT].name     = "FPRD_PPT";
	VarDict[FPRD_PPT].descrip  = "forest production moisture limitation" ;
	VarDict[FPRD_PPT].units    = "unit scalar" ;

	VarDict[FPRD_TMP].name     = "FPRD_TMP";
	VarDict[FPRD_TMP].descrip  = "forest production temperature limitation" ;
	VarDict[FPRD_TMP].units    = "unit scalar" ;

	VarDict[GPRD_PPT].name     = "GPRD_PPT";
	VarDict[GPRD_PPT].descrip  = "grass production moisture limitation" ;
	VarDict[GPRD_PPT].units    = "unit scalar" ;

	VarDict[GPRD_TMP].name     = "GPRD_TMP";
	VarDict[GPRD_TMP].descrip  = "grass production temperature limitation" ;
	VarDict[GPRD_TMP].units    = "unit scalar" ;

	VarDict[SWHC_TOP].name     = "SWHC_TOP";
	VarDict[SWHC_TOP].descrip  = "soil water holding capacity of top soil layer" ;
	VarDict[SWHC_TOP].units    = "" ;
	// VarDict[SWHC_TOP].interval = SINGLE_INTERVAL;

	VarDict[SWHC_MID].name     = "SWHC_MID";
	VarDict[SWHC_MID].descrip  = "soil water holding capacity of middle soil layer" ;
	VarDict[SWHC_MID].units    = "" ;
	// VarDict[SWHC_MID].interval = SINGLE_INTERVAL;

	VarDict[SWHC_DEEP].name     = "SWHC_DEEP";
	VarDict[SWHC_DEEP].descrip  = "soil water holding capacity of deep soil layer" ;
	VarDict[SWHC_DEEP].units    = "" ;
	// VarDict[SWHC_DEEP].interval = SINGLE_INTERVAL;

	makeVarDictEntry(SOIL_DEPTH, "SOIL_DEPTH", "soil depth", "mm", NC_FLOAT, SINGLE_INTERVAL); 
	makeVarDictEntry(ELEV, "ELEV", "elevation", "m", NC_FLOAT, SINGLE_INTERVAL); 
	makeVarDictEntry(RH, "RH", "relative humidity from unsmoothed temperature and vapor pressure", 
			"%", NC_FLOAT, MONTH_INTERVAL);
	makeVarDictEntry(MONTH_OF_FIRE, "MONTH_OF_FIRE", "month in which a fire was simulated", 
			"month, 1-12; 0 for no fire", NC_INT, YEAR_INTERVAL);
	makeVarDictEntry(MINERL_5, "MINERL_5", "mineral N in 5th soil layer", 
			"g N m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(N_VOLATIL, "N_VOLATIL", "N volatilization loss", "g N m-2", NC_FLOAT, MONTH_INTERVAL);
	makeVarDictEntry(PET, "PET", "potential evapotranspiration", "cmH2O", NC_FLOAT, MONTH_INTERVAL);
	makeVarDictEntry(H2O_STREAM_FLOW, "H2O_STREAM_FLOW", "water stream flow", "cmH2O", NC_FLOAT, MONTH_INTERVAL);
	makeVarDictEntry(SFC_RUNOFF, "SFC_RUNOFF", "surface runoff", "cmH2O", NC_FLOAT, MONTH_INTERVAL);
	makeVarDictEntry(MAX_SFC_RUNOFF,"MAX_SFC_RUNOFF", "maximum surface runoff", "cmH2O", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_ECOSYS, "C_ECOSYS", "total ecosystem C storage", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_ECOSYS_DEC, "C_ECOSYS_DEC", "total ecosystem C storage in December", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_LIVE_ABOVEGR, "C_LIVE_ABOVEGR", "C in live aboveground biomass", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_LIVE_BELOWGR, "C_LIVE_BELOWGR", "C in live belowground biomass", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_DEAD_ABOVEGR, "C_DEAD_ABOVEGR", "C in dead aboveground biomass", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_DEAD_BELOWGR, "C_DEAD_BELOWGR", "C in dead belowground biomass", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_MAX_FOREST_LEAF, "C_MAX_FOREST_LEAF", "highest monthly C in forest leaves", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_FINE_BRANCH, "C_FINE_BRANCH", "C in fine branches", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_BOLE, "C_BOLE", "C in large wood", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_MAX_COARSE_ROOT, "C_MAX_COARSE_ROOT", "highest monthly C in coarse roots", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_MAX_FINE_ROOT, "C_MAX_FINE_ROOT", "highest monthly C in fine roots", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_MAX_LIVE_GRASS_BELOWGR, "C_MAX_LIVE_GRASS_BELOWGR", "highest monthly C in live grass belowground", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_LITTER, "C_LITTER", "aboveground litter C for forest", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_LITTER_METAB, "C_LITTER_METAB", "metabolic C in surface litter", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_LITTER_STRUC, "C_LITTER_STRUC", "structural C in surface litter", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_DEAD_WOOD, "C_DEAD_WOOD", "C in large dead wood", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_MAX_STANDING_DEAD, "C_MAX_STANDING_DEAD", "highest monthly value of C in standing dead grass", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_SOIL_AND_LITTER, "C_SOIL_AND_LITTER", "soil and litter organic C", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_SOM_X_STRUC_METAB, "C_SOM_X_STRUC_METAB", "soil C exluding litter and structural C", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_SOM, "C_SOM", "total soil C including belowground structural and metabolic C", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(D1HR, "D1HR", "dead 1-hr fuel", "g DM m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(D10HR, "D10HR", "dead 10-hr fuel", "g DM m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(D100HR, "D100HR", "dead 100-hr fuel", "g DM m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(D1000HR, "D1000HR", "dead 1000-hr fuel", "g DM m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(EM_CO, "EM_CO", "CO emissions from fire", "g CO m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(EM_CO2, "EM_CO2", "CO2 emissions from fire", "g CO2 m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(EM_CH4, "EM_CH4", "CH4 emissions from fire", "g CH4 m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(EM_NMHC, "EM_NMHC", "non-methane hydrocarbon emissions from fire", "g NMHC m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(EM_PM, "EM_PM", "particulate matter emissions from fire", "g DM m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_GRAIN, "C_GRAIN", "C related to grain production", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C_HARVEST, "C_HARVEST", "C removed thru straw during harvest", "g C m-2", NC_FLOAT, MONTH_INTERVAL);
	makeVarDictEntry(TREE_HT, "TREE_HT", "tree height", "m", NC_FLOAT, YEAR_INTERVAL);         
	makeVarDictEntry(C_MAX_LIVE_GRASS, "C_MAX_LIVE_GRASS", "highest monthly C in live grass", "g C m-2", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(C3_PCT_PROD, "C3_PCT_PROD", "% of production from C3 photosynthesis", "%", NC_FLOAT, YEAR_INTERVAL);
	makeCategoricalVarDictEntry(TREE_TYPE, "TREE_TYPE", "1-8 = EN,EN-DB,DB,DB-EB,EN-EB,EB,DN,DN-EN", "tree type", 8);
	makeVarDictEntry(MIN_SMOOTHED_TMP, "MIN_SMOOTHED_TMP", "lowest smoothed monthly mean temperature", "deg C", NC_FLOAT, YEAR_INTERVAL);
	makeCategoricalVarDictEntry(MC_CLIMATE_ZONE, "MC_CLIMATE_ZONE", "1-5: arctic, boreal, temperate, subtropical, tropical", "zone", 5);  
	makeVarDictEntry(TMP_INDEX, "TMP_INDEX", "index of temperature limitation to growth", "0-100: none to total", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(PPT_INDEX, "PPT_INDEX", "index of precipitation limitation to growth", "0-100: none to total", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(NEEDLE_INDEX, "NEEDLE_INDEX", "index of leaf form", "0-100: broadleaf to needleleaf", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(EVERGREEN_INDEX, "EVERGREEN_INDEX", "index of evergreen-ness", "0-100: deciduous to evergreen", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(PSL, "PSL", "precipitation at sea level per PNW-GTR-841", "inH2O", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(SSZ_UB, "SSZ_UB", "Sitka spruce zone upper bound per PNW-GTR-841", "m", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(PSFZ_LB, "PSFZ_LB", "Pacific silver fir zone lower bound per PNW-GTR-841", "m", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(SAFZ_LB, "SAFZ_LB", "subalpine fir zone lower bound per PNW-GTR-841", "m", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(PKLZ_LB, "PKLZ_LB", "subalpine parkland lower bound per PNW-GTR-841", "m", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(MATSL, "MATSL", "mean air temperature at sea level per PNW-GTR-841", "deg F", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(MAX_GRASS_LAI, "MAX_GRASS_LAI", "unsmoothed maximum of monthly grass leaf area index", "m2 leaf m-2 ground", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(MAX_TREE_LAI, "MAX_TREE_LAI", "unsmoothed maximum of monthly tree leaf area index", "m2 leaf m-2 ground", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(SNOWPACK_DAY91, "SNOWPACK_DAY91", "snowpack on Apr 1 (N hemisphere) or Sep 29 (S hemisphere)", "mm H2O", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(CONTINENTALITY, "CONTINENTALITY", "continentality index", "deg C", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(GDD, "GDD", "growing degree-days, referenced to 0 C", "deg C day", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(SOIL_TMP_MAX, "SOIL_TMP_MAX", "estimated max soil temperature", "deg C", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(SWHC_ALL, "SWHC_ALL", "soil water holding capacity, total in all soil layers", "", NC_FLOAT, SINGLE_INTERVAL);
	makeVarDictEntry(FIRE, "FIRE", "fire count", "fires per year", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(FIRE_UNSUPPRESSED, "FIRE_UNSUPPRESSED", "unsuppressed fire count", "unsuppressed fires per year", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(FIRE_FLI, "FIRE_FLI", "fire line intensity of unsuppressed fire", "Btu ft-1 sec-1", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(FIRE_ROS, "FIRE_ROS", "rate of spread of unsuppressed fire", "ft min-1", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(FIRE_ERC, "FIRE_ERC", "energy release component of unsuppressed fire", "Btu ft-2", NC_FLOAT, YEAR_INTERVAL);
	makeCategoricalVarDictEntry(GROUSE_HABITAT, "GROUSE_HABITAT", "grouse habitat: 0=none, non-zero is veg type", "veg type", MAX_VTYPE);
	makeVarDictEntry(GROUSE_SMRPRE, "GROUSE_SMRPRE", "smoothed summer precip (J+J+A or D+J+F)", "mm H2O", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(GROUSE_AUGMAXT, "GROUSE_AUGMAXT", "smoothed Aug (Feb in S hemi) tmax", "deg C", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(GROUSE_ANNTMP, "GROUSE_ANNTMP", "smoothed annual mean tmp", "deg C", NC_FLOAT, YEAR_INTERVAL);
	makeVarDictEntry(TDMEAN, "TDMEAN", "mean dewpoint temperature", 
			"deg C", NC_FLOAT, MONTH_INTERVAL);
	makeVarDictEntry(TMAX_SMOOTHED, "TMAX_SMOOTHED", "smoothed mean value of diurnal max temperature", 
			"deg C", NC_FLOAT, MONTH_INTERVAL);
	makeVarDictEntry(PPT_SMOOTHED, "PPT_SMOOTHED", "smoothed precipitation", 
			"mm H2O", NC_FLOAT, MONTH_INTERVAL);
	makeVarDictEntry(FIRE_MARGIN, "FIRE_MARGIN", "how close to fire thresholds?", 
			"0:1 scale, 1 = far", NC_FLOAT, MONTH_INTERVAL);
	makeVarDictEntry(FIRE_MARGIN_ANN_MIN, "FIRE_MARGIN_ANN_MIN", "how close to fire thresholds?", 
			"0:1 scale, 1 = far", NC_FLOAT, YEAR_INTERVAL);

	writeOutputVarList("optional_output_vars.txt", VarDict, NUM_OPTIONAL_OUTVARS);   
} // end of Simulation::initializeVarDict()


void Simulation::initializeWSoutputList(OutVarDescription f[], int filter_length)
{
	initializeEQoutputList(f, NUM_EQ_OUTVARS);

	// Add fire model state variables and biogeography model state variables here.
	describeFSvars(f);
	f[BSvtype] = describeVar("vtype", "potential vegetation type", "", NC_FLOAT);
	f[BSsmoothed_max_tree_lai] = describeVar("smoothed_max_tree_lai", "smoothed annual maximum tree LAI", 
			"m2 leaf area m-2 ground area", NC_FLOAT);
	f[BSsmoothed_max_grass_lai] = describeVar("smoothed_max_grass_lai", "smoothed annual maximum grass LAI", 
			"m2 leaf area m-2 ground area", NC_FLOAT);

	// Finally, add the warmstart variables at the tail end of the list.
	f[CENtmax_prev_0] = describeVar("tmax_smoothed_0", "smoothed tmax for 1st month", "deg C", NC_FLOAT);
	f[CENtmax_prev_0+1] = describeVar("tmax_smoothed_1", "smoothed tmax for 2nd month", "deg C", NC_FLOAT);
	f[CENtmax_prev_0+2] = describeVar("tmax_smoothed_2", "smoothed tmax for 3rd month", "deg C", NC_FLOAT);
	f[CENtmax_prev_0+3] = describeVar("tmax_smoothed_3", "smoothed tmax for 4th month", "deg C", NC_FLOAT);
	f[CENtmax_prev_0+4] = describeVar("tmax_smoothed_4", "smoothed tmax for 5th month", "deg C", NC_FLOAT);
	f[CENtmax_prev_0+5] = describeVar("tmax_smoothed_5", "smoothed tmax for 6th month", "deg C", NC_FLOAT);
	f[CENtmax_prev_0+6] = describeVar("tmax_smoothed_6", "smoothed tmax for 7th month", "deg C", NC_FLOAT);
	f[CENtmax_prev_0+7] = describeVar("tmax_smoothed_7", "smoothed tmax for 8th month", "deg C", NC_FLOAT);
	f[CENtmax_prev_0+8] = describeVar("tmax_smoothed_8", "smoothed tmax for 9th month", "deg C", NC_FLOAT);
	f[CENtmax_prev_0+9] = describeVar("tmax_smoothed_9", "smoothed tmax for 10th month", "deg C", NC_FLOAT);
	f[CENtmax_prev_0+10] = describeVar("tmax_smoothed_10", "smoothed tmax for 11th month", "deg C", NC_FLOAT);
	f[CENtmax_prev_0+11] = describeVar("tmax_smoothed_11", "smoothed tmax for 12th month", "deg C", NC_FLOAT);

	f[CENtmin_prev_0] = describeVar("tmin_smoothed_0", "smoothed tmin for 1st month", "deg C", NC_FLOAT);
	f[CENtmin_prev_0+1] = describeVar("tmin_smoothed_1", "smoothed tmin for 2nd month", "deg C", NC_FLOAT);
	f[CENtmin_prev_0+2] = describeVar("tmin_smoothed_2", "smoothed tmin for 3rd month", "deg C", NC_FLOAT);
	f[CENtmin_prev_0+3] = describeVar("tmin_smoothed_3", "smoothed tmin for 4th month", "deg C", NC_FLOAT);
	f[CENtmin_prev_0+4] = describeVar("tmin_smoothed_4", "smoothed tmin for 5th month", "deg C", NC_FLOAT);
	f[CENtmin_prev_0+5] = describeVar("tmin_smoothed_5", "smoothed tmin for 6th month", "deg C", NC_FLOAT);
	f[CENtmin_prev_0+6] = describeVar("tmin_smoothed_6", "smoothed tmin for 7th month", "deg C", NC_FLOAT);
	f[CENtmin_prev_0+7] = describeVar("tmin_smoothed_7", "smoothed tmin for 8th month", "deg C", NC_FLOAT);
	f[CENtmin_prev_0+8] = describeVar("tmin_smoothed_8", "smoothed tmin for 9th month", "deg C", NC_FLOAT);
	f[CENtmin_prev_0+9] = describeVar("tmin_smoothed_9", "smoothed tmin for 10th month", "deg C", NC_FLOAT);
	f[CENtmin_prev_0+10] = describeVar("tmin_smoothed_10", "smoothed tmin for 11th month", "deg C", NC_FLOAT);
	f[CENtmin_prev_0+11] = describeVar("tmin_smoothed_11", "smoothed tmin for 12th month", "deg C", NC_FLOAT);

	f[CENtdmean_prev_0] = describeVar("tdmean_smoothed_0", "smoothed tdmean for 1st month", "deg C", NC_FLOAT);
	f[CENtdmean_prev_0+1] = describeVar("tdmean_smoothed_1", "smoothed tdmean for 2nd month", "deg C", NC_FLOAT);
	f[CENtdmean_prev_0+2] = describeVar("tdmean_smoothed_2", "smoothed tdmean for 3rd month", "deg C", NC_FLOAT);
	f[CENtdmean_prev_0+3] = describeVar("tdmean_smoothed_3", "smoothed tdmean for 4th month", "deg C", NC_FLOAT);
	f[CENtdmean_prev_0+4] = describeVar("tdmean_smoothed_4", "smoothed tdmean for 5th month", "deg C", NC_FLOAT);
	f[CENtdmean_prev_0+5] = describeVar("tdmean_smoothed_5", "smoothed tdmean for 6th month", "deg C", NC_FLOAT);
	f[CENtdmean_prev_0+6] = describeVar("tdmean_smoothed_6", "smoothed tdmean for 7th month", "deg C", NC_FLOAT);
	f[CENtdmean_prev_0+7] = describeVar("tdmean_smoothed_7", "smoothed tdmean for 8th month", "deg C", NC_FLOAT);
	f[CENtdmean_prev_0+8] = describeVar("tdmean_smoothed_8", "smoothed tdmean for 9th month", "deg C", NC_FLOAT);
	f[CENtdmean_prev_0+9] = describeVar("tdmean_smoothed_9", "smoothed tdmean for 10th month", "deg C", NC_FLOAT);
	f[CENtdmean_prev_0+10] = describeVar("tdmean_smoothed_10", "smoothed tdmean for 11th month", "deg C", NC_FLOAT);
	f[CENtdmean_prev_0+11] = describeVar("tdmean_smoothed_11", "smoothed tdmean for 12th month", "deg C", NC_FLOAT);

} // end of Simulation::initializeWSoutputList()


void Simulation::initializeEQoutputList(OutVarDescription f[], int filter_length)
{
	for (int i = 0; i<filter_length; i++) f[i].name = "";

	f[CENmapss_canopy] = describeVar("mapss_canopy", "1=C3Dominance, 2=C3C4Mixed, 3=C4Dominance", "", NC_FLOAT);
	f[CENmapss_mix_index] = describeVar("mapss_mix_index", "MAPSS mix_index", "", NC_FLOAT);
	f[CENmapss_tmp_index] = describeVar("mapss_tmp_index", "MAPSS tmp_index", "", NC_FLOAT);
	f[CENmapss_ppt_index] = describeVar("mapss_ppt_index", "MAPSS ppt_index", "", NC_FLOAT);
	f[CENmapss_zone] = describeVar("mapss_zone", "MAPSS climate zone", "", NC_FLOAT);
	f[CENfrost_index_this_year] = describeVar("mapss_frost_index_this_year", "MAPSS frost index", "", NC_FLOAT);
	f[CENfire_last_month] = describeVar("fire_last_month", "a fire was simulated in the previous month", 
			"boolean", NC_FLOAT);
	f[CENequil_time] = describeVar("equil_time", "equil_time", "years", NC_FLOAT);
	f[CENvclass_mapss] = describeVar("vclass_mapss", "vclass from MAPSS", "", NC_FLOAT);
	f[CENclass_mapss] = describeVar("class_mapss", "class from MAPSS", "", NC_FLOAT);

	f[0].name     = "agcacc" ;
	f[0].descrip  = "agcacc" ;
	f[0].units    = "" ;
	f[0].type     = NC_FLOAT ;
	f[0].scale    = 1.000000 ;
	f[0].interval = MONTH_INTERVAL ;

	f[1].name     = "aglcis_1" ;
	f[1].descrip  = "aglcis_1" ;
	f[1].units    = "" ;
	f[1].type     = NC_FLOAT ;
	f[1].scale    = 1.000000 ;
	f[1].interval = MONTH_INTERVAL ;

	f[2].name     = "aglcis_2" ;
	f[2].descrip  = "aglcis_2" ;
	f[2].units    = "" ;
	f[2].type     = NC_FLOAT ;
	f[2].scale    = 1.000000 ;
	f[2].interval = MONTH_INTERVAL ;

	f[3].name     = "aglivc" ;
	f[3].descrip  = "aglivc" ;
	f[3].units    = "" ;
	f[3].type     = NC_FLOAT ;
	f[3].scale    = 1.000000 ;
	f[3].interval = MONTH_INTERVAL ;

	f[4].name     = "aglive_1" ;
	f[4].descrip  = "aglive_1" ;
	f[4].units    = "" ;
	f[4].type     = NC_FLOAT ;
	f[4].scale    = 1.000000 ;
	f[4].interval = MONTH_INTERVAL ;

	f[5].name     = "aglive_2" ;
	f[5].descrip  = "aglive_2" ;
	f[5].units    = "" ;
	f[5].type     = NC_FLOAT ;
	f[5].scale    = 1.000000 ;
	f[5].interval = MONTH_INTERVAL ;

	f[6].name     = "aglive_3" ;
	f[6].descrip  = "aglive_3" ;
	f[6].units    = "" ;
	f[6].type     = NC_FLOAT ;
	f[6].scale    = 1.000000 ;
	f[6].interval = MONTH_INTERVAL ;

	f[7].name     = "aminrl_1" ;
	f[7].descrip  = "aminrl_1" ;
	f[7].units    = "" ;
	f[7].type     = NC_FLOAT ;
	f[7].scale    = 1.000000 ;
	f[7].interval = MONTH_INTERVAL ;

	f[8].name     = "aminrl_2" ;
	f[8].descrip  = "aminrl_2" ;
	f[8].units    = "" ;
	f[8].type     = NC_FLOAT ;
	f[8].scale    = 1.000000 ;
	f[8].interval = MONTH_INTERVAL ;

	f[9].name     = "aminrl_3" ;
	f[9].descrip  = "aminrl_3" ;
	f[9].units    = "" ;
	f[9].type     = NC_FLOAT ;
	f[9].scale    = 1.000000 ;
	f[9].interval = MONTH_INTERVAL ;

	f[10].name     = "amt1c2" ;
	f[10].descrip  = "amt1c2" ;
	f[10].units    = "" ;
	f[10].type     = NC_FLOAT ;
	f[10].scale    = 1.000000 ;
	f[10].interval = MONTH_INTERVAL ;

	f[11].name     = "amt2c2" ;
	f[11].descrip  = "amt2c2" ;
	f[11].units    = "" ;
	f[11].type     = NC_FLOAT ;
	f[11].scale    = 1.000000 ;
	f[11].interval = MONTH_INTERVAL ;

	f[12].name     = "anerb" ;
	f[12].descrip  = "anerb" ;
	f[12].units    = "" ;
	f[12].type     = NC_FLOAT ;
	f[12].scale    = 1.000000 ;
	f[12].interval = MONTH_INTERVAL ;

	f[13].name     = "as11c2" ;
	f[13].descrip  = "as11c2" ;
	f[13].units    = "" ;
	f[13].type     = NC_FLOAT ;
	f[13].scale    = 1.000000 ;
	f[13].interval = MONTH_INTERVAL ;

	f[14].name     = "as21c2" ;
	f[14].descrip  = "as21c2" ;
	f[14].units    = "" ;
	f[14].type     = NC_FLOAT ;
	f[14].scale    = 1.000000 ;
	f[14].interval = MONTH_INTERVAL ;

	f[15].name     = "as2c2" ;
	f[15].descrip  = "as2c2" ;
	f[15].units    = "" ;
	f[15].type     = NC_FLOAT ;
	f[15].scale    = 1.000000 ;
	f[15].interval = MONTH_INTERVAL ;

	f[16].name     = "as3c2" ;
	f[16].descrip  = "as3c2" ;
	f[16].units    = "" ;
	f[16].type     = NC_FLOAT ;
	f[16].scale    = 1.000000 ;
	f[16].interval = MONTH_INTERVAL ;

	f[17].name     = "asmos_1" ;
	f[17].descrip  = "asmos_1" ;
	f[17].units    = "" ;
	f[17].type     = NC_FLOAT ;
	f[17].scale    = 1.000000 ;
	f[17].interval = MONTH_INTERVAL ;

	f[18].name     = "asmos_2" ;
	f[18].descrip  = "asmos_2" ;
	f[18].units    = "" ;
	f[18].type     = NC_FLOAT ;
	f[18].scale    = 1.000000 ;
	f[18].interval = MONTH_INTERVAL ;

	f[19].name     = "asmos_3" ;
	f[19].descrip  = "asmos_3" ;
	f[19].units    = "" ;
	f[19].type     = NC_FLOAT ;
	f[19].scale    = 1.000000 ;
	f[19].interval = MONTH_INTERVAL ;

	f[20].name     = "asmos_4" ;
	f[20].descrip  = "asmos_4" ;
	f[20].units    = "" ;
	f[20].type     = NC_FLOAT ;
	f[20].scale    = 1.000000 ;
	f[20].interval = MONTH_INTERVAL ;

	f[21].name     = "asmos_5" ;
	f[21].descrip  = "asmos_5" ;
	f[21].units    = "" ;
	f[21].type     = NC_FLOAT ;
	f[21].scale    = 1.000000 ;
	f[21].interval = MONTH_INTERVAL ;

	f[22].name     = "asmos_6" ;
	f[22].descrip  = "asmos_6" ;
	f[22].units    = "" ;
	f[22].type     = NC_FLOAT ;
	f[22].scale    = 1.000000 ;
	f[22].interval = MONTH_INTERVAL ;

	f[23].name     = "asmos_7" ;
	f[23].descrip  = "asmos_7" ;
	f[23].units    = "" ;
	f[23].type     = NC_FLOAT ;
	f[23].scale    = 1.000000 ;
	f[23].interval = MONTH_INTERVAL ;

	f[24].name     = "asmos_8" ;
	f[24].descrip  = "asmos_8" ;
	f[24].units    = "" ;
	f[24].type     = NC_FLOAT ;
	f[24].scale    = 1.000000 ;
	f[24].interval = MONTH_INTERVAL ;

	f[25].name     = "asmos_9" ;
	f[25].descrip  = "asmos_9" ;
	f[25].units    = "" ;
	f[25].type     = NC_FLOAT ;
	f[25].scale    = 1.000000 ;
	f[25].interval = MONTH_INTERVAL ;

	f[26].name     = "asmos_10" ;
	f[26].descrip  = "asmos_10" ;
	f[26].units    = "" ;
	f[26].type     = NC_FLOAT ;
	f[26].scale    = 1.000000 ;
	f[26].interval = MONTH_INTERVAL ;

	f[27].name     = "ast1c2" ;
	f[27].descrip  = "ast1c2" ;
	f[27].units    = "" ;
	f[27].type     = NC_FLOAT ;
	f[27].scale    = 1.000000 ;
	f[27].interval = MONTH_INTERVAL ;

	f[28].name     = "ast2c2" ;
	f[28].descrip  = "ast2c2" ;
	f[28].units    = "" ;
	f[28].type     = NC_FLOAT ;
	f[28].scale    = 1.000000 ;
	f[28].interval = MONTH_INTERVAL ;

	f[29].name     = "avh2o_1" ;
	f[29].descrip  = "soil available water for plant growth" ;
	f[29].units    = "cm H2O" ;
	f[29].type     = NC_FLOAT ;
	f[29].scale    = 1.000000 ;
	f[29].interval = MONTH_INTERVAL ;

	f[30].name     = "avh2o_2" ;
	f[30].descrip  = "soil available water for plant survival" ;
	f[30].units    = "cm H2O" ;
	f[30].type     = NC_FLOAT ;
	f[30].scale    = 1.000000 ;
	f[30].interval = MONTH_INTERVAL ;

	f[31].name     = "avh2o_3" ;
	f[31].descrip  = "avh2o_3" ;
	f[31].units    = "" ;
	f[31].type     = NC_FLOAT ;
	f[31].scale    = 1.000000 ;
	f[31].interval = MONTH_INTERVAL ;

	f[32].name     = "bgcacc" ;
	f[32].descrip  = "bgcacc" ;
	f[32].units    = "" ;
	f[32].type     = NC_FLOAT ;
	f[32].scale    = 1.000000 ;
	f[32].interval = MONTH_INTERVAL ;

	f[33].name     = "bglcis_1" ;
	f[33].descrip  = "bglcis_1" ;
	f[33].units    = "" ;
	f[33].type     = NC_FLOAT ;
	f[33].scale    = 1.000000 ;
	f[33].interval = MONTH_INTERVAL ;

	f[34].name     = "bglcis_2" ;
	f[34].descrip  = "bglcis_2" ;
	f[34].units    = "" ;
	f[34].type     = NC_FLOAT ;
	f[34].scale    = 1.000000 ;
	f[34].interval = MONTH_INTERVAL ;

	f[35].name     = "bglivc" ;
	f[35].descrip  = "bglivc" ;
	f[35].units    = "" ;
	f[35].type     = NC_FLOAT ;
	f[35].scale    = 1.000000 ;
	f[35].interval = MONTH_INTERVAL ;

	f[36].name     = "bglive_1" ;
	f[36].descrip  = "bglive_1" ;
	f[36].units    = "" ;
	f[36].type     = NC_FLOAT ;
	f[36].scale    = 1.000000 ;
	f[36].interval = MONTH_INTERVAL ;

	f[37].name     = "bglive_2" ;
	f[37].descrip  = "bglive_2" ;
	f[37].units    = "" ;
	f[37].type     = NC_FLOAT ;
	f[37].scale    = 1.000000 ;
	f[37].interval = MONTH_INTERVAL ;

	f[38].name     = "bglive_3" ;
	f[38].descrip  = "bglive_3" ;
	f[38].units    = "" ;
	f[38].type     = NC_FLOAT ;
	f[38].scale    = 1.000000 ;
	f[38].interval = MONTH_INTERVAL ;

	f[39].name     = "cgrain" ;
	f[39].descrip  = "cgrain" ;
	f[39].units    = "" ;
	f[39].type     = NC_FLOAT ;
	f[39].scale    = 1.000000 ;
	f[39].interval = MONTH_INTERVAL ;

	f[40].name     = "cinput" ;
	f[40].descrip  = "cinput" ;
	f[40].units    = "" ;
	f[40].type     = NC_FLOAT ;
	f[40].scale    = 1.000000 ;
	f[40].interval = MONTH_INTERVAL ;

	f[41].name     = "clittr_1_1" ;
	f[41].descrip  = "clittr_1_1" ;
	f[41].units    = "" ;
	f[41].type     = NC_FLOAT ;
	f[41].scale    = 1.000000 ;
	f[41].interval = MONTH_INTERVAL ;

	f[42].name     = "clittr_2_1" ;
	f[42].descrip  = "clittr_2_1" ;
	f[42].units    = "" ;
	f[42].type     = NC_FLOAT ;
	f[42].scale    = 1.000000 ;
	f[42].interval = MONTH_INTERVAL ;

	f[43].name     = "clittr_1_2" ;
	f[43].descrip  = "clittr_1_2" ;
	f[43].units    = "" ;
	f[43].type     = NC_FLOAT ;
	f[43].scale    = 1.000000 ;
	f[43].interval = MONTH_INTERVAL ;

	f[44].name     = "clittr_2_2" ;
	f[44].descrip  = "clittr_2_2" ;
	f[44].units    = "" ;
	f[44].type     = NC_FLOAT ;
	f[44].scale    = 1.000000 ;
	f[44].interval = MONTH_INTERVAL ;

	f[45].name     = "co2cce_1_1_1" ;
	f[45].descrip  = "co2cce_1_1_1" ;
	f[45].units    = "" ;
	f[45].type     = NC_FLOAT ;
	f[45].scale    = 1.000000 ;
	f[45].interval = MONTH_INTERVAL ;

	f[46].name     = "co2cce_2_1_1" ;
	f[46].descrip  = "co2cce_2_1_1" ;
	f[46].units    = "" ;
	f[46].type     = NC_FLOAT ;
	f[46].scale    = 1.000000 ;
	f[46].interval = MONTH_INTERVAL ;

	f[47].name     = "co2cce_1_2_1" ;
	f[47].descrip  = "co2cce_1_2_1" ;
	f[47].units    = "" ;
	f[47].type     = NC_FLOAT ;
	f[47].scale    = 1.000000 ;
	f[47].interval = MONTH_INTERVAL ;

	f[48].name     = "co2cce_2_2_1" ;
	f[48].descrip  = "co2cce_2_2_1" ;
	f[48].units    = "" ;
	f[48].type     = NC_FLOAT ;
	f[48].scale    = 1.000000 ;
	f[48].interval = MONTH_INTERVAL ;

	f[49].name     = "co2cce_1_1_2" ;
	f[49].descrip  = "co2cce_1_1_2" ;
	f[49].units    = "" ;
	f[49].type     = NC_FLOAT ;
	f[49].scale    = 1.000000 ;
	f[49].interval = MONTH_INTERVAL ;

	f[50].name     = "co2cce_2_1_2" ;
	f[50].descrip  = "co2cce_2_1_2" ;
	f[50].units    = "" ;
	f[50].type     = NC_FLOAT ;
	f[50].scale    = 1.000000 ;
	f[50].interval = MONTH_INTERVAL ;

	f[51].name     = "co2cce_1_2_2" ;
	f[51].descrip  = "co2cce_1_2_2" ;
	f[51].units    = "" ;
	f[51].type     = NC_FLOAT ;
	f[51].scale    = 1.000000 ;
	f[51].interval = MONTH_INTERVAL ;

	f[52].name     = "co2cce_2_2_2" ;
	f[52].descrip  = "co2cce_2_2_2" ;
	f[52].units    = "" ;
	f[52].type     = NC_FLOAT ;
	f[52].scale    = 1.000000 ;
	f[52].interval = MONTH_INTERVAL ;

	f[53].name     = "co2cce_1_1_3" ;
	f[53].descrip  = "co2cce_1_1_3" ;
	f[53].units    = "" ;
	f[53].type     = NC_FLOAT ;
	f[53].scale    = 1.000000 ;
	f[53].interval = MONTH_INTERVAL ;

	f[54].name     = "co2cce_2_1_3" ;
	f[54].descrip  = "co2cce_2_1_3" ;
	f[54].units    = "" ;
	f[54].type     = NC_FLOAT ;
	f[54].scale    = 1.000000 ;
	f[54].interval = MONTH_INTERVAL ;

	f[55].name     = "co2cce_1_2_3" ;
	f[55].descrip  = "co2cce_1_2_3" ;
	f[55].units    = "" ;
	f[55].type     = NC_FLOAT ;
	f[55].scale    = 1.000000 ;
	f[55].interval = MONTH_INTERVAL ;

	f[56].name     = "co2cce_2_2_3" ;
	f[56].descrip  = "co2cce_2_2_3" ;
	f[56].units    = "" ;
	f[56].type     = NC_FLOAT ;
	f[56].scale    = 1.000000 ;
	f[56].interval = MONTH_INTERVAL ;

	f[57].name     = "co2crs_1" ;
	f[57].descrip  = "co2crs_1" ;
	f[57].units    = "" ;
	f[57].type     = NC_FLOAT ;
	f[57].scale    = 1.000000 ;
	f[57].interval = MONTH_INTERVAL ;

	f[58].name     = "co2crs_2" ;
	f[58].descrip  = "co2crs_2" ;
	f[58].units    = "" ;
	f[58].type     = NC_FLOAT ;
	f[58].scale    = 1.000000 ;
	f[58].interval = MONTH_INTERVAL ;

	f[59].name     = "co2cpr_1" ;
	f[59].descrip  = "co2cpr_1" ;
	f[59].units    = "" ;
	f[59].type     = NC_FLOAT ;
	f[59].scale    = 1.000000 ;
	f[59].interval = MONTH_INTERVAL ;

	f[60].name     = "co2cpr_2" ;
	f[60].descrip  = "co2cpr_2" ;
	f[60].units    = "" ;
	f[60].type     = NC_FLOAT ;
	f[60].scale    = 1.000000 ;
	f[60].interval = MONTH_INTERVAL ;

	f[61].name     = "co2ctr_1" ;
	f[61].descrip  = "co2ctr_1" ;
	f[61].units    = "" ;
	f[61].type     = NC_FLOAT ;
	f[61].scale    = 1.000000 ;
	f[61].interval = MONTH_INTERVAL ;

	f[62].name     = "co2ctr_2" ;
	f[62].descrip  = "co2ctr_2" ;
	f[62].units    = "" ;
	f[62].type     = NC_FLOAT ;
	f[62].scale    = 1.000000 ;
	f[62].interval = MONTH_INTERVAL ;

	f[63].name     = "cproda" ;
	f[63].descrip  = "cproda" ;
	f[63].units    = "" ;
	f[63].type     = NC_FLOAT ;
	f[63].scale    = 1.000000 ;
	f[63].interval = MONTH_INTERVAL ;

	f[64].name     = "cprodc" ;
	f[64].descrip  = "cprodc" ;
	f[64].units    = "" ;
	f[64].type     = NC_FLOAT ;
	f[64].scale    = 1.000000 ;
	f[64].interval = MONTH_INTERVAL ;

	f[65].name     = "cprodf" ;
	f[65].descrip  = "cprodf" ;
	f[65].units    = "" ;
	f[65].type     = NC_FLOAT ;
	f[65].scale    = 1.000000 ;
	f[65].interval = MONTH_INTERVAL ;

	f[66].name     = "creta" ;
	f[66].descrip  = "creta" ;
	f[66].units    = "" ;
	f[66].type     = NC_FLOAT ;
	f[66].scale    = 1.000000 ;
	f[66].interval = MONTH_INTERVAL ;

	f[67].name     = "crmvst" ;
	f[67].descrip  = "crmvst" ;
	f[67].units    = "" ;
	f[67].type     = NC_FLOAT ;
	f[67].scale    = 1.000000 ;
	f[67].interval = MONTH_INTERVAL ;

	f[68].name     = "crpstg_1" ;
	f[68].descrip  = "crpstg_1" ;
	f[68].units    = "" ;
	f[68].type     = NC_FLOAT ;
	f[68].scale    = 1.000000 ;
	f[68].interval = MONTH_INTERVAL ;

	f[69].name     = "crpstg_2" ;
	f[69].descrip  = "crpstg_2" ;
	f[69].units    = "" ;
	f[69].type     = NC_FLOAT ;
	f[69].scale    = 1.000000 ;
	f[69].interval = MONTH_INTERVAL ;

	f[70].name     = "crpstg_3" ;
	f[70].descrip  = "crpstg_3" ;
	f[70].units    = "" ;
	f[70].type     = NC_FLOAT ;
	f[70].scale    = 1.000000 ;
	f[70].interval = MONTH_INTERVAL ;

	f[71].name     = "crpval" ;
	f[71].descrip  = "crpval" ;
	f[71].units    = "" ;
	f[71].type     = NC_FLOAT ;
	f[71].scale    = 1.000000 ;
	f[71].interval = MONTH_INTERVAL ;

	f[72].name     = "defac" ;
	f[72].descrip  = "defac" ;
	f[72].units    = "" ;
	f[72].type     = NC_FLOAT ;
	f[72].scale    = 1.000000 ;
	f[72].interval = MONTH_INTERVAL ;

	f[73].name     = "dsomsc" ;
	f[73].descrip  = "dsomsc" ;
	f[73].units    = "" ;
	f[73].type     = NC_FLOAT ;
	f[73].scale    = 1.000000 ;
	f[73].interval = MONTH_INTERVAL ;

	f[74].name     = "egrain_1" ;
	f[74].descrip  = "egrain_1" ;
	f[74].units    = "" ;
	f[74].type     = NC_FLOAT ;
	f[74].scale    = 1.000000 ;
	f[74].interval = MONTH_INTERVAL ;

	f[75].name     = "egrain_2" ;
	f[75].descrip  = "egrain_2" ;
	f[75].units    = "" ;
	f[75].type     = NC_FLOAT ;
	f[75].scale    = 1.000000 ;
	f[75].interval = MONTH_INTERVAL ;

	f[76].name     = "egrain_3" ;
	f[76].descrip  = "egrain_3" ;
	f[76].units    = "" ;
	f[76].type     = NC_FLOAT ;
	f[76].scale    = 1.000000 ;
	f[76].interval = MONTH_INTERVAL ;

	f[77].name     = "elimit" ;
	f[77].descrip  = "elimit" ;
	f[77].units    = "" ;
	f[77].type     = NC_FLOAT ;
	f[77].scale    = 1.000000 ;
	f[77].interval = MONTH_INTERVAL ;

	f[78].name     = "eprodc_1" ;
	f[78].descrip  = "eprodc_1" ;
	f[78].units    = "" ;
	f[78].type     = NC_FLOAT ;
	f[78].scale    = 1.000000 ;
	f[78].interval = MONTH_INTERVAL ;

	f[79].name     = "eprodc_2" ;
	f[79].descrip  = "eprodc_2" ;
	f[79].units    = "" ;
	f[79].type     = NC_FLOAT ;
	f[79].scale    = 1.000000 ;
	f[79].interval = MONTH_INTERVAL ;

	f[80].name     = "eprodc_3" ;
	f[80].descrip  = "eprodc_3" ;
	f[80].units    = "" ;
	f[80].type     = NC_FLOAT ;
	f[80].scale    = 1.000000 ;
	f[80].interval = MONTH_INTERVAL ;

	f[81].name     = "eprodf_1" ;
	f[81].descrip  = "eprodf_1" ;
	f[81].units    = "" ;
	f[81].type     = NC_FLOAT ;
	f[81].scale    = 1.000000 ;
	f[81].interval = MONTH_INTERVAL ;

	f[82].name     = "eprodf_2" ;
	f[82].descrip  = "eprodf_2" ;
	f[82].units    = "" ;
	f[82].type     = NC_FLOAT ;
	f[82].scale    = 1.000000 ;
	f[82].interval = MONTH_INTERVAL ;

	f[83].name     = "eprodf_3" ;
	f[83].descrip  = "eprodf_3" ;
	f[83].units    = "" ;
	f[83].type     = NC_FLOAT ;
	f[83].scale    = 1.000000 ;
	f[83].interval = MONTH_INTERVAL ;

	f[84].name     = "ermvst_1" ;
	f[84].descrip  = "ermvst_1" ;
	f[84].units    = "" ;
	f[84].type     = NC_FLOAT ;
	f[84].scale    = 1.000000 ;
	f[84].interval = MONTH_INTERVAL ;

	f[85].name     = "ermvst_2" ;
	f[85].descrip  = "ermvst_2" ;
	f[85].units    = "" ;
	f[85].type     = NC_FLOAT ;
	f[85].scale    = 1.000000 ;
	f[85].interval = MONTH_INTERVAL ;

	f[86].name     = "ermvst_3" ;
	f[86].descrip  = "ermvst_3" ;
	f[86].units    = "" ;
	f[86].type     = NC_FLOAT ;
	f[86].scale    = 1.000000 ;
	f[86].interval = MONTH_INTERVAL ;

	f[87].name     = "eupacc_1" ;
	f[87].descrip  = "eupacc_1" ;
	f[87].units    = "" ;
	f[87].type     = NC_FLOAT ;
	f[87].scale    = 1.000000 ;
	f[87].interval = MONTH_INTERVAL ;

	f[88].name     = "eupacc_2" ;
	f[88].descrip  = "eupacc_2" ;
	f[88].units    = "" ;
	f[88].type     = NC_FLOAT ;
	f[88].scale    = 1.000000 ;
	f[88].interval = MONTH_INTERVAL ;

	f[89].name     = "eupacc_3" ;
	f[89].descrip  = "eupacc_3" ;
	f[89].units    = "" ;
	f[89].type     = NC_FLOAT ;
	f[89].scale    = 1.000000 ;
	f[89].interval = MONTH_INTERVAL ;

	f[90].name     = "eupaga_1" ;
	f[90].descrip  = "eupaga_1" ;
	f[90].units    = "" ;
	f[90].type     = NC_FLOAT ;
	f[90].scale    = 1.000000 ;
	f[90].interval = MONTH_INTERVAL ;

	f[91].name     = "eupaga_2" ;
	f[91].descrip  = "eupaga_2" ;
	f[91].units    = "" ;
	f[91].type     = NC_FLOAT ;
	f[91].scale    = 1.000000 ;
	f[91].interval = MONTH_INTERVAL ;

	f[92].name     = "eupaga_3" ;
	f[92].descrip  = "eupaga_3" ;
	f[92].units    = "" ;
	f[92].type     = NC_FLOAT ;
	f[92].scale    = 1.000000 ;
	f[92].interval = MONTH_INTERVAL ;

	f[93].name     = "eupbga_1" ;
	f[93].descrip  = "eupbga_1" ;
	f[93].units    = "" ;
	f[93].type     = NC_FLOAT ;
	f[93].scale    = 1.000000 ;
	f[93].interval = MONTH_INTERVAL ;

	f[94].name     = "eupbga_2" ;
	f[94].descrip  = "eupbga_2" ;
	f[94].units    = "" ;
	f[94].type     = NC_FLOAT ;
	f[94].scale    = 1.000000 ;
	f[94].interval = MONTH_INTERVAL ;

	f[95].name     = "eupbga_3" ;
	f[95].descrip  = "eupbga_3" ;
	f[95].units    = "" ;
	f[95].type     = NC_FLOAT ;
	f[95].scale    = 1.000000 ;
	f[95].interval = MONTH_INTERVAL ;

	f[96].name     = "evap" ;
	f[96].descrip  = "monthly evaporation from CENTURY" ;
	f[96].units    = "cmH2O" ;
	f[96].type     = NC_FLOAT ;
	f[96].scale    = 1.000000 ;
	f[96].interval = MONTH_INTERVAL ;

	f[97].name     = "fertot_1" ;
	f[97].descrip  = "fertot_1" ;
	f[97].units    = "" ;
	f[97].type     = NC_FLOAT ;
	f[97].scale    = 1.000000 ;
	f[97].interval = MONTH_INTERVAL ;

	f[98].name     = "fertot_2" ;
	f[98].descrip  = "fertot_2" ;
	f[98].units    = "" ;
	f[98].type     = NC_FLOAT ;
	f[98].scale    = 1.000000 ;
	f[98].interval = MONTH_INTERVAL ;

	f[99].name     = "fertot_3" ;
	f[99].descrip  = "fertot_3" ;
	f[99].units    = "" ;
	f[99].type     = NC_FLOAT ;
	f[99].scale    = 1.000000 ;
	f[99].interval = MONTH_INTERVAL ;

	f[100].name     = "harmth" ;
	f[100].descrip  = "harmth" ;
	f[100].units    = "" ;
	f[100].type     = NC_FLOAT ;
	f[100].scale    = 1.000000 ;
	f[100].interval = MONTH_INTERVAL ;

	f[101].name     = "hi" ;
	f[101].descrip  = "hi" ;
	f[101].units    = "" ;
	f[101].type     = NC_FLOAT ;
	f[101].scale    = 1.000000 ;
	f[101].interval = MONTH_INTERVAL ;

	f[102].name     = "irract" ;
	f[102].descrip  = "irract" ;
	f[102].units    = "" ;
	f[102].type     = NC_FLOAT ;
	f[102].scale    = 1.000000 ;
	f[102].interval = MONTH_INTERVAL ;

	f[103].name     = "irrtot" ;
	f[103].descrip  = "irrtot" ;
	f[103].units    = "" ;
	f[103].type     = NC_FLOAT ;
	f[103].scale    = 1.000000 ;
	f[103].interval = MONTH_INTERVAL ;

	f[104].name     = "metabc_1" ;
	f[104].descrip  = "metabc_1" ;
	f[104].units    = "" ;
	f[104].type     = NC_FLOAT ;
	f[104].scale    = 1.000000 ;
	f[104].interval = MONTH_INTERVAL ;

	f[105].name     = "metabc_2" ;
	f[105].descrip  = "metabc_2" ;
	f[105].units    = "" ;
	f[105].type     = NC_FLOAT ;
	f[105].scale    = 1.000000 ;
	f[105].interval = MONTH_INTERVAL ;

	f[106].name     = "metabe_1_1" ;
	f[106].descrip  = "metabe_1_1" ;
	f[106].units    = "" ;
	f[106].type     = NC_FLOAT ;
	f[106].scale    = 1.000000 ;
	f[106].interval = MONTH_INTERVAL ;

	f[107].name     = "metabe_2_1" ;
	f[107].descrip  = "metabe_2_1" ;
	f[107].units    = "" ;
	f[107].type     = NC_FLOAT ;
	f[107].scale    = 1.000000 ;
	f[107].interval = MONTH_INTERVAL ;

	f[108].name     = "metabe_1_2" ;
	f[108].descrip  = "metabe_1_2" ;
	f[108].units    = "" ;
	f[108].type     = NC_FLOAT ;
	f[108].scale    = 1.000000 ;
	f[108].interval = MONTH_INTERVAL ;

	f[109].name     = "metabe_2_2" ;
	f[109].descrip  = "metabe_2_2" ;
	f[109].units    = "" ;
	f[109].type     = NC_FLOAT ;
	f[109].scale    = 1.000000 ;
	f[109].interval = MONTH_INTERVAL ;

	f[110].name     = "metabe_1_3" ;
	f[110].descrip  = "metabe_1_3" ;
	f[110].units    = "" ;
	f[110].type     = NC_FLOAT ;
	f[110].scale    = 1.000000 ;
	f[110].interval = MONTH_INTERVAL ;

	f[111].name     = "metabe_2_3" ;
	f[111].descrip  = "metabe_2_3" ;
	f[111].units    = "" ;
	f[111].type     = NC_FLOAT ;
	f[111].scale    = 1.000000 ;
	f[111].interval = MONTH_INTERVAL ;

	f[112].name     = "metcis_1_1" ;
	f[112].descrip  = "metcis_1_1" ;
	f[112].units    = "" ;
	f[112].type     = NC_FLOAT ;
	f[112].scale    = 1.000000 ;
	f[112].interval = MONTH_INTERVAL ;

	f[113].name     = "metcis_2_1" ;
	f[113].descrip  = "metcis_2_1" ;
	f[113].units    = "" ;
	f[113].type     = NC_FLOAT ;
	f[113].scale    = 1.000000 ;
	f[113].interval = MONTH_INTERVAL ;

	f[114].name     = "metcis_1_2" ;
	f[114].descrip  = "metcis_1_2" ;
	f[114].units    = "" ;
	f[114].type     = NC_FLOAT ;
	f[114].scale    = 1.000000 ;
	f[114].interval = MONTH_INTERVAL ;

	f[115].name     = "metcis_2_2" ;
	f[115].descrip  = "metcis_2_2" ;
	f[115].units    = "" ;
	f[115].type     = NC_FLOAT ;
	f[115].scale    = 1.000000 ;
	f[115].interval = MONTH_INTERVAL ;

	f[116].name     = "minerl_1_1" ;
	f[116].descrip  = "minerl_1_1" ;
	f[116].units    = "" ;
	f[116].type     = NC_FLOAT ;
	f[116].scale    = 1.000000 ;
	f[116].interval = MONTH_INTERVAL ;

	f[117].name     = "minerl_2_1" ;
	f[117].descrip  = "minerl_2_1" ;
	f[117].units    = "" ;
	f[117].type     = NC_FLOAT ;
	f[117].scale    = 1.000000 ;
	f[117].interval = MONTH_INTERVAL ;

	f[118].name     = "minerl_3_1" ;
	f[118].descrip  = "minerl_3_1" ;
	f[118].units    = "" ;
	f[118].type     = NC_FLOAT ;
	f[118].scale    = 1.000000 ;
	f[118].interval = MONTH_INTERVAL ;

	f[119].name     = "minerl_4_1" ;
	f[119].descrip  = "minerl_4_1" ;
	f[119].units    = "" ;
	f[119].type     = NC_FLOAT ;
	f[119].scale    = 1.000000 ;
	f[119].interval = MONTH_INTERVAL ;

	f[120].name     = "minerl_5_1" ;
	f[120].descrip  = "minerl_5_1" ;
	f[120].units    = "" ;
	f[120].type     = NC_FLOAT ;
	f[120].scale    = 1.000000 ;
	f[120].interval = MONTH_INTERVAL ;

	f[121].name     = "minerl_6_1" ;
	f[121].descrip  = "minerl_6_1" ;
	f[121].units    = "" ;
	f[121].type     = NC_FLOAT ;
	f[121].scale    = 1.000000 ;
	f[121].interval = MONTH_INTERVAL ;

	f[122].name     = "minerl_7_1" ;
	f[122].descrip  = "minerl_7_1" ;
	f[122].units    = "" ;
	f[122].type     = NC_FLOAT ;
	f[122].scale    = 1.000000 ;
	f[122].interval = MONTH_INTERVAL ;

	f[123].name     = "minerl_8_1" ;
	f[123].descrip  = "minerl_8_1" ;
	f[123].units    = "" ;
	f[123].type     = NC_FLOAT ;
	f[123].scale    = 1.000000 ;
	f[123].interval = MONTH_INTERVAL ;

	f[124].name     = "minerl_9_1" ;
	f[124].descrip  = "minerl_9_1" ;
	f[124].units    = "" ;
	f[124].type     = NC_FLOAT ;
	f[124].scale    = 1.000000 ;
	f[124].interval = MONTH_INTERVAL ;

	f[125].name     = "minerl_10_1" ;
	f[125].descrip  = "minerl_10_1" ;
	f[125].units    = "" ;
	f[125].type     = NC_FLOAT ;
	f[125].scale    = 1.000000 ;
	f[125].interval = MONTH_INTERVAL ;

	f[126].name     = "minerl_11_1" ;
	f[126].descrip  = "minerl_11_1" ;
	f[126].units    = "" ;
	f[126].type     = NC_FLOAT ;
	f[126].scale    = 1.000000 ;
	f[126].interval = MONTH_INTERVAL ;

	f[127].name     = "minerl_1_2" ;
	f[127].descrip  = "minerl_1_2" ;
	f[127].units    = "" ;
	f[127].type     = NC_FLOAT ;
	f[127].scale    = 1.000000 ;
	f[127].interval = MONTH_INTERVAL ;

	f[128].name     = "minerl_2_2" ;
	f[128].descrip  = "minerl_2_2" ;
	f[128].units    = "" ;
	f[128].type     = NC_FLOAT ;
	f[128].scale    = 1.000000 ;
	f[128].interval = MONTH_INTERVAL ;

	f[129].name     = "minerl_3_2" ;
	f[129].descrip  = "minerl_3_2" ;
	f[129].units    = "" ;
	f[129].type     = NC_FLOAT ;
	f[129].scale    = 1.000000 ;
	f[129].interval = MONTH_INTERVAL ;

	f[130].name     = "minerl_4_2" ;
	f[130].descrip  = "minerl_4_2" ;
	f[130].units    = "" ;
	f[130].type     = NC_FLOAT ;
	f[130].scale    = 1.000000 ;
	f[130].interval = MONTH_INTERVAL ;

	f[131].name     = "minerl_5_2" ;
	f[131].descrip  = "minerl_5_2" ;
	f[131].units    = "" ;
	f[131].type     = NC_FLOAT ;
	f[131].scale    = 1.000000 ;
	f[131].interval = MONTH_INTERVAL ;

	f[132].name     = "minerl_6_2" ;
	f[132].descrip  = "minerl_6_2" ;
	f[132].units    = "" ;
	f[132].type     = NC_FLOAT ;
	f[132].scale    = 1.000000 ;
	f[132].interval = MONTH_INTERVAL ;

	f[133].name     = "minerl_7_2" ;
	f[133].descrip  = "minerl_7_2" ;
	f[133].units    = "" ;
	f[133].type     = NC_FLOAT ;
	f[133].scale    = 1.000000 ;
	f[133].interval = MONTH_INTERVAL ;

	f[134].name     = "minerl_8_2" ;
	f[134].descrip  = "minerl_8_2" ;
	f[134].units    = "" ;
	f[134].type     = NC_FLOAT ;
	f[134].scale    = 1.000000 ;
	f[134].interval = MONTH_INTERVAL ;

	f[135].name     = "minerl_9_2" ;
	f[135].descrip  = "minerl_9_2" ;
	f[135].units    = "" ;
	f[135].type     = NC_FLOAT ;
	f[135].scale    = 1.000000 ;
	f[135].interval = MONTH_INTERVAL ;

	f[136].name     = "minerl_10_2" ;
	f[136].descrip  = "minerl_10_2" ;
	f[136].units    = "" ;
	f[136].type     = NC_FLOAT ;
	f[136].scale    = 1.000000 ;
	f[136].interval = MONTH_INTERVAL ;

	f[137].name     = "minerl_11_2" ;
	f[137].descrip  = "minerl_11_2" ;
	f[137].units    = "" ;
	f[137].type     = NC_FLOAT ;
	f[137].scale    = 1.000000 ;
	f[137].interval = MONTH_INTERVAL ;

	// f[138] - f[145] spare

	f[146].name     = "mt1c2_1" ;
	f[146].descrip  = "mt1c2_1" ;
	f[146].units    = "" ;
	f[146].type     = NC_FLOAT ;
	f[146].scale    = 1.000000 ;
	f[146].interval = MONTH_INTERVAL ;

	f[147].name     = "mt1c2_2" ;
	f[147].descrip  = "mt1c2_2" ;
	f[147].units    = "" ;
	f[147].type     = NC_FLOAT ;
	f[147].scale    = 1.000000 ;
	f[147].interval = MONTH_INTERVAL ;

	f[148].name     = "mt2c2_1" ;
	f[148].descrip  = "mt2c2_1" ;
	f[148].units    = "" ;
	f[148].type     = NC_FLOAT ;
	f[148].scale    = 1.000000 ;
	f[148].interval = MONTH_INTERVAL ;

	f[149].name     = "mt2c2_2" ;
	f[149].descrip  = "mt2c2_2" ;
	f[149].units    = "" ;
	f[149].type     = NC_FLOAT ;
	f[149].scale    = 1.000000 ;
	f[149].interval = MONTH_INTERVAL ;

	f[150].name     = "nfix" ;
	f[150].descrip  = "nfix" ;
	f[150].units    = "" ;
	f[150].type     = NC_FLOAT ;
	f[150].scale    = 1.000000 ;
	f[150].interval = MONTH_INTERVAL ;

	f[151].name     = "nfixac" ;
	f[151].descrip  = "nfixac" ;
	f[151].units    = "" ;
	f[151].type     = NC_FLOAT ;
	f[151].scale    = 1.000000 ;
	f[151].interval = MONTH_INTERVAL ;

	f[152].name     = "occlud" ;
	f[152].descrip  = "occlud" ;
	f[152].units    = "" ;
	f[152].type     = NC_FLOAT ;
	f[152].scale    = 1.000000 ;
	f[152].interval = MONTH_INTERVAL ;

	f[153].name     = "parent_1" ;
	f[153].descrip  = "parent_1" ;
	f[153].units    = "" ;
	f[153].type     = NC_FLOAT ;
	f[153].scale    = 1.000000 ;
	f[153].interval = MONTH_INTERVAL ;

	f[154].name     = "parent_2" ;
	f[154].descrip  = "parent_2" ;
	f[154].units    = "" ;
	f[154].type     = NC_FLOAT ;
	f[154].scale    = 1.000000 ;
	f[154].interval = MONTH_INTERVAL ;

	f[155].name     = "parent_3" ;
	f[155].descrip  = "parent_3" ;
	f[155].units    = "" ;
	f[155].type     = NC_FLOAT ;
	f[155].scale    = 1.000000 ;
	f[155].interval = MONTH_INTERVAL ;

	f[156].name     = "pet_century" ;
	f[156].descrip  = "monthly potential evapotranspiration from CENTURY" ;
	f[156].units    = "cm H2O" ;
	f[156].type     = NC_FLOAT ;
	f[156].scale    = 1.000000 ;
	f[156].interval = MONTH_INTERVAL ;

	f[157].name     = "petann" ;
	f[157].descrip  = "year-to-date potential evapotranspiration from CENTURY" ;
	f[157].units    = "cm H2O" ;
	f[157].type     = NC_FLOAT ;
	f[157].scale    = 1.000000 ;
	f[157].interval = MONTH_INTERVAL ;

	f[158].name     = "plabil" ;
	f[158].descrip  = "plabil" ;
	f[158].units    = "" ;
	f[158].type     = NC_FLOAT ;
	f[158].scale    = 1.000000 ;
	f[158].interval = MONTH_INTERVAL ;

	f[159].name     = "prcann" ;
	f[159].descrip  = "prcann" ;
	f[159].units    = "" ;
	f[159].type     = NC_FLOAT ;
	f[159].scale    = 1.000000 ;
	f[159].interval = MONTH_INTERVAL ;

	f[160].name     = "ptagc" ;
	f[160].descrip  = "ptagc" ;
	f[160].units    = "" ;
	f[160].type     = NC_FLOAT ;
	f[160].scale    = 1.000000 ;
	f[160].interval = MONTH_INTERVAL ;

	f[161].name     = "ptbgc" ;
	f[161].descrip  = "ptbgc" ;
	f[161].units    = "" ;
	f[161].type     = NC_FLOAT ;
	f[161].scale    = 1.000000 ;
	f[161].interval = MONTH_INTERVAL ;

	f[162].name     = "pttr" ;
	f[162].descrip  = "pttr" ;
	f[162].units    = "" ;
	f[162].type     = NC_FLOAT ;
	f[162].scale    = 1.000000 ;
	f[162].interval = MONTH_INTERVAL ;

	f[163].name     = "rain" ;
	f[163].descrip  = "rain" ;
	f[163].units    = "" ;
	f[163].type     = NC_FLOAT ;
	f[163].scale    = 1.000000 ;
	f[163].interval = MONTH_INTERVAL ;

	f[164].name     = "resp_1" ;
	f[164].descrip  = "resp_1" ;
	f[164].units    = "" ;
	f[164].type     = NC_FLOAT ;
	f[164].scale    = 1.000000 ;
	f[164].interval = MONTH_INTERVAL ;

	f[165].name     = "resp_2" ;
	f[165].descrip  = "resp_2" ;
	f[165].units    = "" ;
	f[165].type     = NC_FLOAT ;
	f[165].scale    = 1.000000 ;
	f[165].interval = MONTH_INTERVAL ;

	f[166].name     = "relyld" ;
	f[166].descrip  = "relyld" ;
	f[166].units    = "" ;
	f[166].type     = NC_FLOAT ;
	f[166].scale    = 1.000000 ;
	f[166].interval = MONTH_INTERVAL ;

	f[167].name     = "rwcf_1" ;
	f[167].descrip  = "rwcf_1" ;
	f[167].units    = "" ;
	f[167].type     = NC_FLOAT ;
	f[167].scale    = 1.000000 ;
	f[167].interval = MONTH_INTERVAL ;

	f[168].name     = "rwcf_2" ;
	f[168].descrip  = "rwcf_2" ;
	f[168].units    = "" ;
	f[168].type     = NC_FLOAT ;
	f[168].scale    = 1.000000 ;
	f[168].interval = MONTH_INTERVAL ;

	f[169].name     = "rwcf_3" ;
	f[169].descrip  = "rwcf_3" ;
	f[169].units    = "" ;
	f[169].type     = NC_FLOAT ;
	f[169].scale    = 1.000000 ;
	f[169].interval = MONTH_INTERVAL ;

	f[170].name     = "rwcf_4" ;
	f[170].descrip  = "rwcf_4" ;
	f[170].units    = "" ;
	f[170].type     = NC_FLOAT ;
	f[170].scale    = 1.000000 ;
	f[170].interval = MONTH_INTERVAL ;

	f[171].name     = "rwcf_5" ;
	f[171].descrip  = "rwcf_5" ;
	f[171].units    = "" ;
	f[171].type     = NC_FLOAT ;
	f[171].scale    = 1.000000 ;
	f[171].interval = MONTH_INTERVAL ;

	f[172].name     = "rwcf_6" ;
	f[172].descrip  = "rwcf_6" ;
	f[172].units    = "" ;
	f[172].type     = NC_FLOAT ;
	f[172].scale    = 1.000000 ;
	f[172].interval = MONTH_INTERVAL ;

	f[173].name     = "rwcf_7" ;
	f[173].descrip  = "rwcf_7" ;
	f[173].units    = "" ;
	f[173].type     = NC_FLOAT ;
	f[173].scale    = 1.000000 ;
	f[173].interval = MONTH_INTERVAL ;

	f[174].name     = "rwcf_8" ;
	f[174].descrip  = "rwcf_8" ;
	f[174].units    = "" ;
	f[174].type     = NC_FLOAT ;
	f[174].scale    = 1.000000 ;
	f[174].interval = MONTH_INTERVAL ;

	f[175].name     = "rwcf_9" ;
	f[175].descrip  = "rwcf_9" ;
	f[175].units    = "" ;
	f[175].type     = NC_FLOAT ;
	f[175].scale    = 1.000000 ;
	f[175].interval = MONTH_INTERVAL ;

	f[176].name     = "rwcf_10" ;
	f[176].descrip  = "rwcf_10" ;
	f[176].units    = "" ;
	f[176].type     = NC_FLOAT ;
	f[176].scale    = 1.000000 ;
	f[176].interval = MONTH_INTERVAL ;

	f[177].name     = "s11c2_1" ;
	f[177].descrip  = "s11c2_1" ;
	f[177].units    = "" ;
	f[177].type     = NC_FLOAT ;
	f[177].scale    = 1.000000 ;
	f[177].interval = MONTH_INTERVAL ;

	f[178].name     = "s11c2_2" ;
	f[178].descrip  = "s11c2_2" ;
	f[178].units    = "" ;
	f[178].type     = NC_FLOAT ;
	f[178].scale    = 1.000000 ;
	f[178].interval = MONTH_INTERVAL ;

	f[179].name     = "s21c2_1" ;
	f[179].descrip  = "s21c2_1" ;
	f[179].units    = "" ;
	f[179].type     = NC_FLOAT ;
	f[179].scale    = 1.000000 ;
	f[179].interval = MONTH_INTERVAL ;

	f[180].name     = "s21c2_2" ;
	f[180].descrip  = "s21c2_2" ;
	f[180].units    = "" ;
	f[180].type     = NC_FLOAT ;
	f[180].scale    = 1.000000 ;
	f[180].interval = MONTH_INTERVAL ;

	f[181].name     = "s2c2_1" ;
	f[181].descrip  = "s2c2_1" ;
	f[181].units    = "" ;
	f[181].type     = NC_FLOAT ;
	f[181].scale    = 1.000000 ;
	f[181].interval = MONTH_INTERVAL ;

	f[182].name     = "s2c2_2" ;
	f[182].descrip  = "s2c2_2" ;
	f[182].units    = "" ;
	f[182].type     = NC_FLOAT ;
	f[182].scale    = 1.000000 ;
	f[182].interval = MONTH_INTERVAL ;

	f[183].name     = "s3c2_1" ;
	f[183].descrip  = "s3c2_1" ;
	f[183].units    = "" ;
	f[183].type     = NC_FLOAT ;
	f[183].scale    = 1.000000 ;
	f[183].interval = MONTH_INTERVAL ;

	f[184].name     = "s3c2_2" ;
	f[184].descrip  = "s3c2_2" ;
	f[184].units    = "" ;
	f[184].type     = NC_FLOAT ;
	f[184].scale    = 1.000000 ;
	f[184].interval = MONTH_INTERVAL ;

	f[185].name     = "satmac" ;
	f[185].descrip  = "satmac" ;
	f[185].units    = "" ;
	f[185].type     = NC_FLOAT ;
	f[185].scale    = 1.000000 ;
	f[185].interval = MONTH_INTERVAL ;

	f[186].name     = "sclosa" ;
	f[186].descrip  = "sclosa" ;
	f[186].units    = "" ;
	f[186].type     = NC_FLOAT ;
	f[186].scale    = 1.000000 ;
	f[186].interval = MONTH_INTERVAL ;

	f[187].name     = "scloss" ;
	f[187].descrip  = "scloss" ;
	f[187].units    = "" ;
	f[187].type     = NC_FLOAT ;
	f[187].scale    = 1.000000 ;
	f[187].interval = MONTH_INTERVAL ;

	f[188].name     = "sdrema" ;
	f[188].descrip  = "sdrema" ;
	f[188].units    = "" ;
	f[188].type     = NC_FLOAT ;
	f[188].scale    = 1.000000 ;
	f[188].interval = MONTH_INTERVAL ;

	f[189].name     = "secndy_1" ;
	f[189].descrip  = "secndy_1" ;
	f[189].units    = "" ;
	f[189].type     = NC_FLOAT ;
	f[189].scale    = 1.000000 ;
	f[189].interval = MONTH_INTERVAL ;

	f[190].name     = "secndy_2" ;
	f[190].descrip  = "secndy_2" ;
	f[190].units    = "" ;
	f[190].type     = NC_FLOAT ;
	f[190].scale    = 1.000000 ;
	f[190].interval = MONTH_INTERVAL ;

	f[191].name     = "secndy_3" ;
	f[191].descrip  = "secndy_3" ;
	f[191].units    = "" ;
	f[191].type     = NC_FLOAT ;
	f[191].scale    = 1.000000 ;
	f[191].interval = MONTH_INTERVAL ;

	f[192].name     = "shrema" ;
	f[192].descrip  = "shrema" ;
	f[192].units    = "" ;
	f[192].type     = NC_FLOAT ;
	f[192].scale    = 1.000000 ;
	f[192].interval = MONTH_INTERVAL ;

	f[193].name     = "sirrac" ;
	f[193].descrip  = "sirrac" ;
	f[193].units    = "" ;
	f[193].type     = NC_FLOAT ;
	f[193].scale    = 1.000000 ;
	f[193].interval = MONTH_INTERVAL ;

	f[194].name     = "snfxac_1" ;
	f[194].descrip  = "snfxac_1" ;
	f[194].units    = "" ;
	f[194].type     = NC_FLOAT ;
	f[194].scale    = 1.000000 ;
	f[194].interval = MONTH_INTERVAL ;

	f[195].name     = "snfxac_2" ;
	f[195].descrip  = "snfxac_2" ;
	f[195].units    = "" ;
	f[195].type     = NC_FLOAT ;
	f[195].scale    = 1.000000 ;
	f[195].interval = MONTH_INTERVAL ;

	f[196].name     = "snlq" ;
	f[196].descrip  = "snlq" ;
	f[196].units    = "" ;
	f[196].type     = NC_FLOAT ;
	f[196].scale    = 1.000000 ;
	f[196].interval = MONTH_INTERVAL ;

	f[197].name     = "snow" ;
	f[197].descrip  = "snow" ;
	f[197].units    = "" ;
	f[197].type     = NC_FLOAT ;
	f[197].scale    = 1.000000 ;
	f[197].interval = MONTH_INTERVAL ;

	f[198].name     = "soilnm_1" ;
	f[198].descrip  = "soilnm_1" ;
	f[198].units    = "" ;
	f[198].type     = NC_FLOAT ;
	f[198].scale    = 1.000000 ;
	f[198].interval = MONTH_INTERVAL ;

	f[199].name     = "soilnm_2" ;
	f[199].descrip  = "soilnm_2" ;
	f[199].units    = "" ;
	f[199].type     = NC_FLOAT ;
	f[199].scale    = 1.000000 ;
	f[199].interval = MONTH_INTERVAL ;

	f[200].name     = "soilnm_3" ;
	f[200].descrip  = "soilnm_3" ;
	f[200].units    = "" ;
	f[200].type     = NC_FLOAT ;
	f[200].scale    = 1.000000 ;
	f[200].interval = MONTH_INTERVAL ;

	f[201].name     = "somsc" ;
	f[201].descrip  = "somsc" ;
	f[201].units    = "" ;
	f[201].type     = NC_FLOAT ;
	f[201].scale    = 1.000000 ;
	f[201].interval = MONTH_INTERVAL ;

	f[202].name     = "somse_1" ;
	f[202].descrip  = "somse_1" ;
	f[202].units    = "" ;
	f[202].type     = NC_FLOAT ;
	f[202].scale    = 1.000000 ;
	f[202].interval = MONTH_INTERVAL ;

	f[203].name     = "somse_2" ;
	f[203].descrip  = "somse_2" ;
	f[203].units    = "" ;
	f[203].type     = NC_FLOAT ;
	f[203].scale    = 1.000000 ;
	f[203].interval = MONTH_INTERVAL ;

	f[204].name     = "somse_3" ;
	f[204].descrip  = "somse_3" ;
	f[204].units    = "" ;
	f[204].type     = NC_FLOAT ;
	f[204].scale    = 1.000000 ;
	f[204].interval = MONTH_INTERVAL ;

	f[205].name     = "somtc" ;
	f[205].descrip  = "somtc" ;
	f[205].units    = "" ;
	f[205].type     = NC_FLOAT ;
	f[205].scale    = 1.000000 ;
	f[205].interval = MONTH_INTERVAL ;

	f[206].name     = "som1c_1" ;
	f[206].descrip  = "som1c_1" ;
	f[206].units    = "" ;
	f[206].type     = NC_FLOAT ;
	f[206].scale    = 1.000000 ;
	f[206].interval = MONTH_INTERVAL ;

	f[207].name     = "som1c_2" ;
	f[207].descrip  = "som1c_2" ;
	f[207].units    = "" ;
	f[207].type     = NC_FLOAT ;
	f[207].scale    = 1.000000 ;
	f[207].interval = MONTH_INTERVAL ;

	f[208].name     = "som1ci_1_1" ;
	f[208].descrip  = "som1ci_1_1" ;
	f[208].units    = "" ;
	f[208].type     = NC_FLOAT ;
	f[208].scale    = 1.000000 ;
	f[208].interval = MONTH_INTERVAL ;

	f[209].name     = "som1ci_2_1" ;
	f[209].descrip  = "som1ci_2_1" ;
	f[209].units    = "" ;
	f[209].type     = NC_FLOAT ;
	f[209].scale    = 1.000000 ;
	f[209].interval = MONTH_INTERVAL ;

	f[210].name     = "som1ci_1_2" ;
	f[210].descrip  = "som1ci_1_2" ;
	f[210].units    = "" ;
	f[210].type     = NC_FLOAT ;
	f[210].scale    = 1.000000 ;
	f[210].interval = MONTH_INTERVAL ;

	f[211].name     = "som1ci_2_2" ;
	f[211].descrip  = "som1ci_2_2" ;
	f[211].units    = "" ;
	f[211].type     = NC_FLOAT ;
	f[211].scale    = 1.000000 ;
	f[211].interval = MONTH_INTERVAL ;

	f[212].name     = "som1e_1_1" ;
	f[212].descrip  = "som1e_1_1" ;
	f[212].units    = "" ;
	f[212].type     = NC_FLOAT ;
	f[212].scale    = 1.000000 ;
	f[212].interval = MONTH_INTERVAL ;

	f[213].name     = "som1e_2_1" ;
	f[213].descrip  = "som1e_2_1" ;
	f[213].units    = "" ;
	f[213].type     = NC_FLOAT ;
	f[213].scale    = 1.000000 ;
	f[213].interval = MONTH_INTERVAL ;

	f[214].name     = "som1e_1_2" ;
	f[214].descrip  = "som1e_1_2" ;
	f[214].units    = "" ;
	f[214].type     = NC_FLOAT ;
	f[214].scale    = 1.000000 ;
	f[214].interval = MONTH_INTERVAL ;

	f[215].name     = "som1e_2_2" ;
	f[215].descrip  = "som1e_2_2" ;
	f[215].units    = "" ;
	f[215].type     = NC_FLOAT ;
	f[215].scale    = 1.000000 ;
	f[215].interval = MONTH_INTERVAL ;

	f[216].name     = "som1e_1_3" ;
	f[216].descrip  = "som1e_1_3" ;
	f[216].units    = "" ;
	f[216].type     = NC_FLOAT ;
	f[216].scale    = 1.000000 ;
	f[216].interval = MONTH_INTERVAL ;

	f[217].name     = "som1e_2_3" ;
	f[217].descrip  = "som1e_2_3" ;
	f[217].units    = "" ;
	f[217].type     = NC_FLOAT ;
	f[217].scale    = 1.000000 ;
	f[217].interval = MONTH_INTERVAL ;

	f[218].name     = "som2c" ;
	f[218].descrip  = "som2c" ;
	f[218].units    = "" ;
	f[218].type     = NC_FLOAT ;
	f[218].scale    = 1.000000 ;
	f[218].interval = MONTH_INTERVAL ;

	f[219].name     = "som2ci_1" ;
	f[219].descrip  = "som2ci_1" ;
	f[219].units    = "" ;
	f[219].type     = NC_FLOAT ;
	f[219].scale    = 1.000000 ;
	f[219].interval = MONTH_INTERVAL ;

	f[220].name     = "som2ci_2" ;
	f[220].descrip  = "som2ci_2" ;
	f[220].units    = "" ;
	f[220].type     = NC_FLOAT ;
	f[220].scale    = 1.000000 ;
	f[220].interval = MONTH_INTERVAL ;

	f[221].name     = "som2e_1" ;
	f[221].descrip  = "som2e_1" ;
	f[221].units    = "" ;
	f[221].type     = NC_FLOAT ;
	f[221].scale    = 1.000000 ;
	f[221].interval = MONTH_INTERVAL ;

	f[222].name     = "som2e_2" ;
	f[222].descrip  = "som2e_2" ;
	f[222].units    = "" ;
	f[222].type     = NC_FLOAT ;
	f[222].scale    = 1.000000 ;
	f[222].interval = MONTH_INTERVAL ;

	f[223].name     = "som2e_3" ;
	f[223].descrip  = "som2e_3" ;
	f[223].units    = "" ;
	f[223].type     = NC_FLOAT ;
	f[223].scale    = 1.000000 ;
	f[223].interval = MONTH_INTERVAL ;

	f[224].name     = "som3c" ;
	f[224].descrip  = "som3c" ;
	f[224].units    = "" ;
	f[224].type     = NC_FLOAT ;
	f[224].scale    = 1.000000 ;
	f[224].interval = MONTH_INTERVAL ;

	f[225].name     = "som3ci_1" ;
	f[225].descrip  = "som3ci_1" ;
	f[225].units    = "" ;
	f[225].type     = NC_FLOAT ;
	f[225].scale    = 1.000000 ;
	f[225].interval = MONTH_INTERVAL ;

	f[226].name     = "som3ci_2" ;
	f[226].descrip  = "som3ci_2" ;
	f[226].units    = "" ;
	f[226].type     = NC_FLOAT ;
	f[226].scale    = 1.000000 ;
	f[226].interval = MONTH_INTERVAL ;

	f[227].name     = "som3e_1" ;
	f[227].descrip  = "som3e_1" ;
	f[227].units    = "" ;
	f[227].type     = NC_FLOAT ;
	f[227].scale    = 1.000000 ;
	f[227].interval = MONTH_INTERVAL ;

	f[228].name     = "som3e_2" ;
	f[228].descrip  = "som3e_2" ;
	f[228].units    = "" ;
	f[228].type     = NC_FLOAT ;
	f[228].scale    = 1.000000 ;
	f[228].interval = MONTH_INTERVAL ;

	f[229].name     = "som3e_3" ;
	f[229].descrip  = "som3e_3" ;
	f[229].units    = "" ;
	f[229].type     = NC_FLOAT ;
	f[229].scale    = 1.000000 ;
	f[229].interval = MONTH_INTERVAL ;

	f[230].name     = "stdcis_1" ;
	f[230].descrip  = "stdcis_1" ;
	f[230].units    = "" ;
	f[230].type     = NC_FLOAT ;
	f[230].scale    = 1.000000 ;
	f[230].interval = MONTH_INTERVAL ;

	f[231].name     = "stdcis_2" ;
	f[231].descrip  = "stdcis_2" ;
	f[231].units    = "" ;
	f[231].type     = NC_FLOAT ;
	f[231].scale    = 1.000000 ;
	f[231].interval = MONTH_INTERVAL ;

	f[232].name     = "st1c2_1" ;
	f[232].descrip  = "st1c2_1" ;
	f[232].units    = "" ;
	f[232].type     = NC_FLOAT ;
	f[232].scale    = 1.000000 ;
	f[232].interval = MONTH_INTERVAL ;

	f[233].name     = "st1c2_2" ;
	f[233].descrip  = "st1c2_2" ;
	f[233].units    = "" ;
	f[233].type     = NC_FLOAT ;
	f[233].scale    = 1.000000 ;
	f[233].interval = MONTH_INTERVAL ;

	f[234].name     = "st2c2_1" ;
	f[234].descrip  = "st2c2_1" ;
	f[234].units    = "" ;
	f[234].type     = NC_FLOAT ;
	f[234].scale    = 1.000000 ;
	f[234].interval = MONTH_INTERVAL ;

	f[235].name     = "st2c2_2" ;
	f[235].descrip  = "st2c2_2" ;
	f[235].units    = "" ;
	f[235].type     = NC_FLOAT ;
	f[235].scale    = 1.000000 ;
	f[235].interval = MONTH_INTERVAL ;

	f[236].name     = "stemp" ;
	f[236].descrip  = "stemp" ;
	f[236].units    = "" ;
	f[236].type     = NC_FLOAT ;
	f[236].scale    = 1.000000 ;
	f[236].interval = MONTH_INTERVAL ;

	f[237].name     = "strcis_1_1" ;
	f[237].descrip  = "strcis_1_1" ;
	f[237].units    = "" ;
	f[237].type     = NC_FLOAT ;
	f[237].scale    = 1.000000 ;
	f[237].interval = MONTH_INTERVAL ;

	f[238].name     = "strcis_2_1" ;
	f[238].descrip  = "strcis_2_1" ;
	f[238].units    = "" ;
	f[238].type     = NC_FLOAT ;
	f[238].scale    = 1.000000 ;
	f[238].interval = MONTH_INTERVAL ;

	f[239].name     = "strcis_1_2" ;
	f[239].descrip  = "strcis_1_2" ;
	f[239].units    = "" ;
	f[239].type     = NC_FLOAT ;
	f[239].scale    = 1.000000 ;
	f[239].interval = MONTH_INTERVAL ;

	f[240].name     = "strcis_2_2" ;
	f[240].descrip  = "strcis_2_2" ;
	f[240].units    = "" ;
	f[240].type     = NC_FLOAT ;
	f[240].scale    = 1.000000 ;
	f[240].interval = MONTH_INTERVAL ;

	f[241].name     = "stream_1" ;
	f[241].descrip  = "stream_1" ;
	f[241].units    = "" ;
	f[241].type     = NC_FLOAT ;
	f[241].scale    = 1.000000 ;
	f[241].interval = MONTH_INTERVAL ;

	f[242].name     = "stream_2" ;
	f[242].descrip  = "stream_2" ;
	f[242].units    = "" ;
	f[242].type     = NC_FLOAT ;
	f[242].scale    = 1.000000 ;
	f[242].interval = MONTH_INTERVAL ;

	f[243].name     = "stream_3" ;
	f[243].descrip  = "stream_3" ;
	f[243].units    = "" ;
	f[243].type     = NC_FLOAT ;
	f[243].scale    = 1.000000 ;
	f[243].interval = MONTH_INTERVAL ;

	f[244].name     = "stream_4" ;
	f[244].descrip  = "stream_4" ;
	f[244].units    = "" ;
	f[244].type     = NC_FLOAT ;
	f[244].scale    = 1.000000 ;
	f[244].interval = MONTH_INTERVAL ;

	f[245].name     = "stream_5" ;
	f[245].descrip  = "stream_5" ;
	f[245].units    = "" ;
	f[245].type     = NC_FLOAT ;
	f[245].scale    = 1.000000 ;
	f[245].interval = MONTH_INTERVAL ;

	f[246].name     = "stream_6" ;
	f[246].descrip  = "stream_6" ;
	f[246].units    = "" ;
	f[246].type     = NC_FLOAT ;
	f[246].scale    = 1.000000 ;
	f[246].interval = MONTH_INTERVAL ;

	f[247].name     = "stream_7" ;
	f[247].descrip  = "stream_7" ;
	f[247].units    = "" ;
	f[247].type     = NC_FLOAT ;
	f[247].scale    = 1.000000 ;
	f[247].interval = MONTH_INTERVAL ;

	f[248].name     = "stream_8" ;
	f[248].descrip  = "stream_8" ;
	f[248].units    = "" ;
	f[248].type     = NC_FLOAT ;
	f[248].scale    = 1.000000 ;
	f[248].interval = MONTH_INTERVAL ;

	f[249].name     = "strlig_1" ;
	f[249].descrip  = "strlig_1" ;
	f[249].units    = "" ;
	f[249].type     = NC_FLOAT ;
	f[249].scale    = 1.000000 ;
	f[249].interval = MONTH_INTERVAL ;

	f[250].name     = "strlig_2" ;
	f[250].descrip  = "strlig_2" ;
	f[250].units    = "" ;
	f[250].type     = NC_FLOAT ;
	f[250].scale    = 1.000000 ;
	f[250].interval = MONTH_INTERVAL ;

	f[251].name     = "strucc_1" ;
	f[251].descrip  = "strucc_1" ;
	f[251].units    = "" ;
	f[251].type     = NC_FLOAT ;
	f[251].scale    = 1.000000 ;
	f[251].interval = MONTH_INTERVAL ;

	f[252].name     = "strucc_2" ;
	f[252].descrip  = "strucc_2" ;
	f[252].units    = "" ;
	f[252].type     = NC_FLOAT ;
	f[252].scale    = 1.000000 ;
	f[252].interval = MONTH_INTERVAL ;

	f[253].name     = "struce_1_1" ;
	f[253].descrip  = "struce_1_1" ;
	f[253].units    = "" ;
	f[253].type     = NC_FLOAT ;
	f[253].scale    = 1.000000 ;
	f[253].interval = MONTH_INTERVAL ;

	f[254].name     = "struce_2_1" ;
	f[254].descrip  = "struce_2_1" ;
	f[254].units    = "" ;
	f[254].type     = NC_FLOAT ;
	f[254].scale    = 1.000000 ;
	f[254].interval = MONTH_INTERVAL ;

	f[255].name     = "struce_1_2" ;
	f[255].descrip  = "struce_1_2" ;
	f[255].units    = "" ;
	f[255].type     = NC_FLOAT ;
	f[255].scale    = 1.000000 ;
	f[255].interval = MONTH_INTERVAL ;

	f[256].name     = "struce_2_2" ;
	f[256].descrip  = "struce_2_2" ;
	f[256].units    = "" ;
	f[256].type     = NC_FLOAT ;
	f[256].scale    = 1.000000 ;
	f[256].interval = MONTH_INTERVAL ;

	f[257].name     = "struce_1_3" ;
	f[257].descrip  = "struce_1_3" ;
	f[257].units    = "" ;
	f[257].type     = NC_FLOAT ;
	f[257].scale    = 1.000000 ;
	f[257].interval = MONTH_INTERVAL ;

	f[258].name     = "struce_2_3" ;
	f[258].descrip  = "struce_2_3" ;
	f[258].units    = "" ;
	f[258].type     = NC_FLOAT ;
	f[258].scale    = 1.000000 ;
	f[258].interval = MONTH_INTERVAL ;

	f[259].name     = "sumnrs_1" ;
	f[259].descrip  = "sumnrs_1" ;
	f[259].units    = "" ;
	f[259].type     = NC_FLOAT ;
	f[259].scale    = 1.000000 ;
	f[259].interval = MONTH_INTERVAL ;

	f[260].name     = "sumnrs_2" ;
	f[260].descrip  = "sumnrs_2" ;
	f[260].units    = "" ;
	f[260].type     = NC_FLOAT ;
	f[260].scale    = 1.000000 ;
	f[260].interval = MONTH_INTERVAL ;

	f[261].name     = "sumnrs_3" ;
	f[261].descrip  = "sumnrs_3" ;
	f[261].units    = "" ;
	f[261].type     = NC_FLOAT ;
	f[261].scale    = 1.000000 ;
	f[261].interval = MONTH_INTERVAL ;

	f[262].name     = "stdedc" ;
	f[262].descrip  = "stdedc" ;
	f[262].units    = "" ;
	f[262].type     = NC_FLOAT ;
	f[262].scale    = 1.000000 ;
	f[262].interval = MONTH_INTERVAL ;

	f[263].name     = "stdede_1" ;
	f[263].descrip  = "stdede_1" ;
	f[263].units    = "" ;
	f[263].type     = NC_FLOAT ;
	f[263].scale    = 1.000000 ;
	f[263].interval = MONTH_INTERVAL ;

	f[264].name     = "stdede_2" ;
	f[264].descrip  = "stdede_2" ;
	f[264].units    = "" ;
	f[264].type     = NC_FLOAT ;
	f[264].scale    = 1.000000 ;
	f[264].interval = MONTH_INTERVAL ;

	f[265].name     = "stdede_3" ;
	f[265].descrip  = "stdede_3" ;
	f[265].units    = "" ;
	f[265].type     = NC_FLOAT ;
	f[265].scale    = 1.000000 ;
	f[265].interval = MONTH_INTERVAL ;

	f[266].name     = "tave" ;
	f[266].descrip  = "tave" ;
	f[266].units    = "" ;
	f[266].type     = NC_FLOAT ;
	f[266].scale    = 1.000000 ;
	f[266].interval = MONTH_INTERVAL ;

	f[267].name     = "tminrl_1" ;
	f[267].descrip  = "tminrl_1" ;
	f[267].units    = "" ;
	f[267].type     = NC_FLOAT ;
	f[267].scale    = 1.000000 ;
	f[267].interval = MONTH_INTERVAL ;

	f[268].name     = "tminrl_2" ;
	f[268].descrip  = "tminrl_2" ;
	f[268].units    = "" ;
	f[268].type     = NC_FLOAT ;
	f[268].scale    = 1.000000 ;
	f[268].interval = MONTH_INTERVAL ;

	f[269].name     = "tminrl_3" ;
	f[269].descrip  = "tminrl_3" ;
	f[269].units    = "" ;
	f[269].type     = NC_FLOAT ;
	f[269].scale    = 1.000000 ;
	f[269].interval = MONTH_INTERVAL ;

	f[270].name     = "tnetmn_1" ;
	f[270].descrip  = "tnetmn_1" ;
	f[270].units    = "" ;
	f[270].type     = NC_FLOAT ;
	f[270].scale    = 1.000000 ;
	f[270].interval = MONTH_INTERVAL ;

	f[271].name     = "tnetmn_2" ;
	f[271].descrip  = "tnetmn_2" ;
	f[271].units    = "" ;
	f[271].type     = NC_FLOAT ;
	f[271].scale    = 1.000000 ;
	f[271].interval = MONTH_INTERVAL ;

	f[272].name     = "tnetmn_3" ;
	f[272].descrip  = "tnetmn_3" ;
	f[272].units    = "" ;
	f[272].type     = NC_FLOAT ;
	f[272].scale    = 1.000000 ;
	f[272].interval = MONTH_INTERVAL ;

	f[273].name     = "totc" ;
	f[273].descrip  = "Century's totc (NOT total ecosystem C)" ;
	f[273].units    = "g C m-2" ;
	f[273].type     = NC_FLOAT ;
	f[273].scale    = 1.000000 ;
	f[273].interval = MONTH_INTERVAL ;

	f[274].name     = "tran" ;
	f[274].descrip  = "monthly transpiration" ;
	f[274].units    = "cmH2O" ;
	f[274].type     = NC_FLOAT ;
	f[274].scale    = 1.000000 ;
	f[274].interval = MONTH_INTERVAL ;

	f[275].name     = "volgma" ;
	f[275].descrip  = "volgma" ;
	f[275].units    = "" ;
	f[275].type     = NC_FLOAT ;
	f[275].scale    = 1.000000 ;
	f[275].interval = MONTH_INTERVAL ;

	f[276].name     = "volexa" ;
	f[276].descrip  = "volexa" ;
	f[276].units    = "" ;
	f[276].type     = NC_FLOAT ;
	f[276].scale    = 1.000000 ;
	f[276].interval = MONTH_INTERVAL ;

	f[277].name     = "volpla" ;
	f[277].descrip  = "volpla" ;
	f[277].units    = "" ;
	f[277].type     = NC_FLOAT ;
	f[277].scale    = 1.000000 ;
	f[277].interval = MONTH_INTERVAL ;

	f[278].name     = "wdfxaa" ;
	f[278].descrip  = "wdfxaa" ;
	f[278].units    = "" ;
	f[278].type     = NC_FLOAT ;
	f[278].scale    = 1.000000 ;
	f[278].interval = MONTH_INTERVAL ;

	f[279].name     = "wdfxas" ;
	f[279].descrip  = "wdfxas" ;
	f[279].units    = "" ;
	f[279].type     = NC_FLOAT ;
	f[279].scale    = 1.000000 ;
	f[279].interval = MONTH_INTERVAL ;

	f[280].name     = "accrst" ;
	f[280].descrip  = "accrst" ;
	f[280].units    = "" ;
	f[280].type     = NC_FLOAT ;
	f[280].scale    = 1.000000 ;
	f[280].interval = MONTH_INTERVAL ;

	f[281].name     = "adefac" ;
	f[281].descrip  = "adefac" ;
	f[281].units    = "" ;
	f[281].type     = NC_FLOAT ;
	f[281].scale    = 1.000000 ;
	f[281].interval = MONTH_INTERVAL ;

	f[282].name     = "agcisa_1" ;
	f[282].descrip  = "agcisa_1" ;
	f[282].units    = "" ;
	f[282].type     = NC_FLOAT ;
	f[282].scale    = 1.000000 ;
	f[282].interval = MONTH_INTERVAL ;

	f[283].name     = "agcisa_2" ;
	f[283].descrip  = "agcisa_2" ;
	f[283].units    = "" ;
	f[283].type     = NC_FLOAT ;
	f[283].scale    = 1.000000 ;
	f[283].interval = MONTH_INTERVAL ;

	f[284].name     = "aglcn" ;
	f[284].descrip  = "aglcn" ;
	f[284].units    = "" ;
	f[284].type     = NC_FLOAT ;
	f[284].scale    = 1.000000 ;
	f[284].interval = MONTH_INTERVAL ;

	f[285].name     = "bgcisa_1" ;
	f[285].descrip  = "bgcisa_1" ;
	f[285].units    = "" ;
	f[285].type     = NC_FLOAT ;
	f[285].scale    = 1.000000 ;
	f[285].interval = MONTH_INTERVAL ;

	f[286].name     = "bgcisa_2" ;
	f[286].descrip  = "bgcisa_2" ;
	f[286].units    = "" ;
	f[286].type     = NC_FLOAT ;
	f[286].scale    = 1.000000 ;
	f[286].interval = MONTH_INTERVAL ;

	f[287].name     = "bglcn" ;
	f[287].descrip  = "bglcn" ;
	f[287].units    = "" ;
	f[287].type     = NC_FLOAT ;
	f[287].scale    = 1.000000 ;
	f[287].interval = MONTH_INTERVAL ;

	f[288].name     = "cgracc" ;
	f[288].descrip  = "cgracc" ;
	f[288].units    = "" ;
	f[288].type     = NC_FLOAT ;
	f[288].scale    = 1.000000 ;
	f[288].interval = MONTH_INTERVAL ;

	f[289].name     = "cisgra_1" ;
	f[289].descrip  = "cisgra_1" ;
	f[289].units    = "" ;
	f[289].type     = NC_FLOAT ;
	f[289].scale    = 1.000000 ;
	f[289].interval = MONTH_INTERVAL ;

	f[290].name     = "cisgra_2" ;
	f[290].descrip  = "cisgra_2" ;
	f[290].units    = "" ;
	f[290].type     = NC_FLOAT ;
	f[290].scale    = 1.000000 ;
	f[290].interval = MONTH_INTERVAL ;

	f[291].name     = "cltfac_1" ;
	f[291].descrip  = "cltfac_1" ;
	f[291].units    = "" ;
	f[291].type     = NC_FLOAT ;
	f[291].scale    = 1.000000 ;
	f[291].interval = MONTH_INTERVAL ;

	f[292].name     = "cltfac_2" ;
	f[292].descrip  = "cltfac_2" ;
	f[292].units    = "" ;
	f[292].type     = NC_FLOAT ;
	f[292].scale    = 1.000000 ;
	f[292].interval = MONTH_INTERVAL ;

	f[293].name     = "cltfac_3" ;
	f[293].descrip  = "cltfac_3" ;
	f[293].units    = "" ;
	f[293].type     = NC_FLOAT ;
	f[293].scale    = 1.000000 ;
	f[293].interval = MONTH_INTERVAL ;

	f[294].name     = "cltfac_4" ;
	f[294].descrip  = "cltfac_4" ;
	f[294].units    = "" ;
	f[294].type     = NC_FLOAT ;
	f[294].scale    = 1.000000 ;
	f[294].interval = MONTH_INTERVAL ;

	f[295].name     = "csrsnk_1" ;
	f[295].descrip  = "csrsnk_1" ;
	f[295].units    = "" ;
	f[295].type     = NC_FLOAT ;
	f[295].scale    = 1.000000 ;
	f[295].interval = MONTH_INTERVAL ;

	f[296].name     = "csrsnk_2" ;
	f[296].descrip  = "csrsnk_2" ;
	f[296].units    = "" ;
	f[296].type     = NC_FLOAT ;
	f[296].scale    = 1.000000 ;
	f[296].interval = MONTH_INTERVAL ;

	f[297].name     = "dblit" ;
	f[297].descrip  = "dblit" ;
	f[297].units    = "" ;
	f[297].type     = NC_FLOAT ;
	f[297].scale    = 1.000000 ;
	f[297].interval = MONTH_INTERVAL ;

	f[298].name     = "dmetc_1" ;
	f[298].descrip  = "dmetc_1" ;
	f[298].units    = "" ;
	f[298].type     = NC_FLOAT ;
	f[298].scale    = 1.000000 ;
	f[298].interval = MONTH_INTERVAL ;

	f[299].name     = "dmetc_2" ;
	f[299].descrip  = "dmetc_2" ;
	f[299].units    = "" ;
	f[299].type     = NC_FLOAT ;
	f[299].scale    = 1.000000 ;
	f[299].interval = MONTH_INTERVAL ;

	f[300].name     = "dslit" ;
	f[300].descrip  = "dslit" ;
	f[300].units    = "" ;
	f[300].type     = NC_FLOAT ;
	f[300].scale    = 1.000000 ;
	f[300].interval = MONTH_INTERVAL ;

	f[301].name     = "dsom1c_1" ;
	f[301].descrip  = "dsom1c_1" ;
	f[301].units    = "" ;
	f[301].type     = NC_FLOAT ;
	f[301].scale    = 1.000000 ;
	f[301].interval = MONTH_INTERVAL ;

	f[302].name     = "dsom1c_2" ;
	f[302].descrip  = "dsom1c_2" ;
	f[302].units    = "" ;
	f[302].type     = NC_FLOAT ;
	f[302].scale    = 1.000000 ;
	f[302].interval = MONTH_INTERVAL ;

	f[303].name     = "dsom2c" ;
	f[303].descrip  = "dsom2c" ;
	f[303].units    = "" ;
	f[303].type     = NC_FLOAT ;
	f[303].scale    = 1.000000 ;
	f[303].interval = MONTH_INTERVAL ;

	f[304].name     = "dsom3c" ;
	f[304].descrip  = "dsom3c" ;
	f[304].units    = "" ;
	f[304].type     = NC_FLOAT ;
	f[304].scale    = 1.000000 ;
	f[304].interval = MONTH_INTERVAL ;

	f[305].name     = "dsomtc" ;
	f[305].descrip  = "dsomtc" ;
	f[305].units    = "" ;
	f[305].type     = NC_FLOAT ;
	f[305].scale    = 1.000000 ;
	f[305].interval = MONTH_INTERVAL ;

	f[306].name     = "dstruc_1" ;
	f[306].descrip  = "dstruc_1" ;
	f[306].units    = "" ;
	f[306].type     = NC_FLOAT ;
	f[306].scale    = 1.000000 ;
	f[306].interval = MONTH_INTERVAL ;

	f[307].name     = "dstruc_2" ;
	f[307].descrip  = "dstruc_2" ;
	f[307].units    = "" ;
	f[307].type     = NC_FLOAT ;
	f[307].scale    = 1.000000 ;
	f[307].interval = MONTH_INTERVAL ;

	f[308].name     = "egracc_1" ;
	f[308].descrip  = "egracc_1" ;
	f[308].units    = "" ;
	f[308].type     = NC_FLOAT ;
	f[308].scale    = 1.000000 ;
	f[308].interval = MONTH_INTERVAL ;

	f[309].name     = "egracc_2" ;
	f[309].descrip  = "egracc_2" ;
	f[309].units    = "" ;
	f[309].type     = NC_FLOAT ;
	f[309].scale    = 1.000000 ;
	f[309].interval = MONTH_INTERVAL ;

	f[310].name     = "egracc_3" ;
	f[310].descrip  = "egracc_3" ;
	f[310].units    = "" ;
	f[310].type     = NC_FLOAT ;
	f[310].scale    = 1.000000 ;
	f[310].interval = MONTH_INTERVAL ;

	f[311].name     = "ereta_1" ;
	f[311].descrip  = "ereta_1" ;
	f[311].units    = "" ;
	f[311].type     = NC_FLOAT ;
	f[311].scale    = 1.000000 ;
	f[311].interval = MONTH_INTERVAL ;

	f[312].name     = "ereta_2" ;
	f[312].descrip  = "ereta_2" ;
	f[312].units    = "" ;
	f[312].type     = NC_FLOAT ;
	f[312].scale    = 1.000000 ;
	f[312].interval = MONTH_INTERVAL ;

	f[313].name     = "ereta_3" ;
	f[313].descrip  = "ereta_3" ;
	f[313].units    = "" ;
	f[313].type     = NC_FLOAT ;
	f[313].scale    = 1.000000 ;
	f[313].interval = MONTH_INTERVAL ;

	f[314].name     = "esrsnk_1" ;
	f[314].descrip  = "esrsnk_1" ;
	f[314].units    = "" ;
	f[314].type     = NC_FLOAT ;
	f[314].scale    = 1.000000 ;
	f[314].interval = MONTH_INTERVAL ;

	f[315].name     = "esrsnk_2" ;
	f[315].descrip  = "esrsnk_2" ;
	f[315].units    = "" ;
	f[315].type     = NC_FLOAT ;
	f[315].scale    = 1.000000 ;
	f[315].interval = MONTH_INTERVAL ;

	f[316].name     = "esrsnk_3" ;
	f[316].descrip  = "esrsnk_3" ;
	f[316].units    = "" ;
	f[316].type     = NC_FLOAT ;
	f[316].scale    = 1.000000 ;
	f[316].interval = MONTH_INTERVAL ;

	f[317].name     = "gromin_1" ;
	f[317].descrip  = "gromin_1" ;
	f[317].units    = "" ;
	f[317].type     = NC_FLOAT ;
	f[317].scale    = 1.000000 ;
	f[317].interval = MONTH_INTERVAL ;

	f[318].name     = "gromin_2" ;
	f[318].descrip  = "gromin_2" ;
	f[318].units    = "" ;
	f[318].type     = NC_FLOAT ;
	f[318].scale    = 1.000000 ;
	f[318].interval = MONTH_INTERVAL ;

	f[319].name     = "gromin_3" ;
	f[319].descrip  = "gromin_3" ;
	f[319].units    = "" ;
	f[319].type     = NC_FLOAT ;
	f[319].scale    = 1.000000 ;
	f[319].interval = MONTH_INTERVAL ;

	f[320].name     = "lhzcac" ;
	f[320].descrip  = "lhzcac" ;
	f[320].units    = "" ;
	f[320].type     = NC_FLOAT ;
	f[320].scale    = 1.000000 ;
	f[320].interval = MONTH_INTERVAL ;

	f[321].name     = "lhzeac_1" ;
	f[321].descrip  = "lhzeac_1" ;
	f[321].units    = "" ;
	f[321].type     = NC_FLOAT ;
	f[321].scale    = 1.000000 ;
	f[321].interval = MONTH_INTERVAL ;

	f[322].name     = "lhzeac_2" ;
	f[322].descrip  = "lhzeac_2" ;
	f[322].units    = "" ;
	f[322].type     = NC_FLOAT ;
	f[322].scale    = 1.000000 ;
	f[322].interval = MONTH_INTERVAL ;

	f[323].name     = "lhzeac_3" ;
	f[323].descrip  = "lhzeac_3" ;
	f[323].units    = "" ;
	f[323].type     = NC_FLOAT ;
	f[323].scale    = 1.000000 ;
	f[323].interval = MONTH_INTERVAL ;

	f[324].name     = "metmnr_1_1" ;
	f[324].descrip  = "metmnr_1_1" ;
	f[324].units    = "" ;
	f[324].type     = NC_FLOAT ;
	f[324].scale    = 1.000000 ;
	f[324].interval = MONTH_INTERVAL ;

	f[325].name     = "metmnr_2_1" ;
	f[325].descrip  = "metmnr_2_1" ;
	f[325].units    = "" ;
	f[325].type     = NC_FLOAT ;
	f[325].scale    = 1.000000 ;
	f[325].interval = MONTH_INTERVAL ;

	f[326].name     = "metmnr_1_2" ;
	f[326].descrip  = "metmnr_1_2" ;
	f[326].units    = "" ;
	f[326].type     = NC_FLOAT ;
	f[326].scale    = 1.000000 ;
	f[326].interval = MONTH_INTERVAL ;

	f[327].name     = "metmnr_2_2" ;
	f[327].descrip  = "metmnr_2_2" ;
	f[327].units    = "" ;
	f[327].type     = NC_FLOAT ;
	f[327].scale    = 1.000000 ;
	f[327].interval = MONTH_INTERVAL ;

	f[328].name     = "metmnr_1_3" ;
	f[328].descrip  = "metmnr_1_3" ;
	f[328].units    = "" ;
	f[328].type     = NC_FLOAT ;
	f[328].scale    = 1.000000 ;
	f[328].interval = MONTH_INTERVAL ;

	f[329].name     = "metmnr_2_3" ;
	f[329].descrip  = "metmnr_2_3" ;
	f[329].units    = "" ;
	f[329].type     = NC_FLOAT ;
	f[329].scale    = 1.000000 ;
	f[329].interval = MONTH_INTERVAL ;

	f[330].name     = "prcfal" ;
	f[330].descrip  = "prcfal" ;
	f[330].units    = "" ;
	f[330].type     = NC_FLOAT ;
	f[330].scale    = 1.000000 ;
	f[330].interval = MONTH_INTERVAL ;

	f[331].name     = "rmvsti_1" ;
	f[331].descrip  = "rmvsti_1" ;
	f[331].units    = "" ;
	f[331].type     = NC_FLOAT ;
	f[331].scale    = 1.000000 ;
	f[331].interval = MONTH_INTERVAL ;

	f[332].name     = "rmvsti_2" ;
	f[332].descrip  = "rmvsti_2" ;
	f[332].units    = "" ;
	f[332].type     = NC_FLOAT ;
	f[332].scale    = 1.000000 ;
	f[332].interval = MONTH_INTERVAL ;

	f[333].name     = "rnpml1" ;
	f[333].descrip  = "rnpml1" ;
	f[333].units    = "" ;
	f[333].type     = NC_FLOAT ;
	f[333].scale    = 1.000000 ;
	f[333].interval = MONTH_INTERVAL ;

	f[334].name     = "sdrmae_1" ;
	f[334].descrip  = "sdrmae_1" ;
	f[334].units    = "" ;
	f[334].type     = NC_FLOAT ;
	f[334].scale    = 1.000000 ;
	f[334].interval = MONTH_INTERVAL ;

	f[335].name     = "sdrmae_2" ;
	f[335].descrip  = "sdrmae_2" ;
	f[335].units    = "" ;
	f[335].type     = NC_FLOAT ;
	f[335].scale    = 1.000000 ;
	f[335].interval = MONTH_INTERVAL ;

	f[336].name     = "sdrmae_3" ;
	f[336].descrip  = "sdrmae_3" ;
	f[336].units    = "" ;
	f[336].type     = NC_FLOAT ;
	f[336].scale    = 1.000000 ;
	f[336].interval = MONTH_INTERVAL ;

	f[337].name     = "sdrmai_1" ;
	f[337].descrip  = "sdrmai_1" ;
	f[337].units    = "" ;
	f[337].type     = NC_FLOAT ;
	f[337].scale    = 1.000000 ;
	f[337].interval = MONTH_INTERVAL ;

	f[338].name     = "sdrmai_2" ;
	f[338].descrip  = "sdrmai_2" ;
	f[338].units    = "" ;
	f[338].type     = NC_FLOAT ;
	f[338].scale    = 1.000000 ;
	f[338].interval = MONTH_INTERVAL ;

	f[339].name     = "shrmai_1" ;
	f[339].descrip  = "shrmai_1" ;
	f[339].units    = "" ;
	f[339].type     = NC_FLOAT ;
	f[339].scale    = 1.000000 ;
	f[339].interval = MONTH_INTERVAL ;

	f[340].name     = "shrmai_2" ;
	f[340].descrip  = "shrmai_2" ;
	f[340].units    = "" ;
	f[340].type     = NC_FLOAT ;
	f[340].scale    = 1.000000 ;
	f[340].interval = MONTH_INTERVAL ;

	f[341].name     = "shrmae_1" ;
	f[341].descrip  = "shrmae_1" ;
	f[341].units    = "" ;
	f[341].type     = NC_FLOAT ;
	f[341].scale    = 1.000000 ;
	f[341].interval = MONTH_INTERVAL ;

	f[342].name     = "shrmae_2" ;
	f[342].descrip  = "shrmae_2" ;
	f[342].units    = "" ;
	f[342].type     = NC_FLOAT ;
	f[342].scale    = 1.000000 ;
	f[342].interval = MONTH_INTERVAL ;

	f[343].name     = "shrmae_3" ;
	f[343].descrip  = "shrmae_3" ;
	f[343].units    = "" ;
	f[343].type     = NC_FLOAT ;
	f[343].scale    = 1.000000 ;
	f[343].interval = MONTH_INTERVAL ;

	f[344].name     = "somsci_1" ;
	f[344].descrip  = "somsci_1" ;
	f[344].units    = "" ;
	f[344].type     = NC_FLOAT ;
	f[344].scale    = 1.000000 ;
	f[344].interval = MONTH_INTERVAL ;

	f[345].name     = "somsci_2" ;
	f[345].descrip  = "somsci_2" ;
	f[345].units    = "" ;
	f[345].type     = NC_FLOAT ;
	f[345].scale    = 1.000000 ;
	f[345].interval = MONTH_INTERVAL ;

	f[346].name     = "somtci_1" ;
	f[346].descrip  = "somtci_1" ;
	f[346].units    = "" ;
	f[346].type     = NC_FLOAT ;
	f[346].scale    = 1.000000 ;
	f[346].interval = MONTH_INTERVAL ;

	f[347].name     = "somtci_2" ;
	f[347].descrip  = "somtci_2" ;
	f[347].units    = "" ;
	f[347].type     = NC_FLOAT ;
	f[347].scale    = 1.000000 ;
	f[347].interval = MONTH_INTERVAL ;

	f[348].name     = "somte_1" ;
	f[348].descrip  = "somte_1" ;
	f[348].units    = "" ;
	f[348].type     = NC_FLOAT ;
	f[348].scale    = 1.000000 ;
	f[348].interval = MONTH_INTERVAL ;

	f[349].name     = "somte_2" ;
	f[349].descrip  = "somte_2" ;
	f[349].units    = "" ;
	f[349].type     = NC_FLOAT ;
	f[349].scale    = 1.000000 ;
	f[349].interval = MONTH_INTERVAL ;

	f[350].name     = "somte_3" ;
	f[350].descrip  = "somte_3" ;
	f[350].units    = "" ;
	f[350].type     = NC_FLOAT ;
	f[350].scale    = 1.000000 ;
	f[350].interval = MONTH_INTERVAL ;

	f[351].name     = "strmnr_1_1" ;
	f[351].descrip  = "strmnr_1_1" ;
	f[351].units    = "" ;
	f[351].type     = NC_FLOAT ;
	f[351].scale    = 1.000000 ;
	f[351].interval = MONTH_INTERVAL ;

	f[352].name     = "strmnr_2_1" ;
	f[352].descrip  = "strmnr_2_1" ;
	f[352].units    = "" ;
	f[352].type     = NC_FLOAT ;
	f[352].scale    = 1.000000 ;
	f[352].interval = MONTH_INTERVAL ;

	f[353].name     = "strmnr_1_2" ;
	f[353].descrip  = "strmnr_1_2" ;
	f[353].units    = "" ;
	f[353].type     = NC_FLOAT ;
	f[353].scale    = 1.000000 ;
	f[353].interval = MONTH_INTERVAL ;

	f[354].name     = "strmnr_2_2" ;
	f[354].descrip  = "strmnr_2_2" ;
	f[354].units    = "" ;
	f[354].type     = NC_FLOAT ;
	f[354].scale    = 1.000000 ;
	f[354].interval = MONTH_INTERVAL ;

	f[355].name     = "strmnr_1_3" ;
	f[355].descrip  = "strmnr_1_3" ;
	f[355].units    = "" ;
	f[355].type     = NC_FLOAT ;
	f[355].scale    = 1.000000 ;
	f[355].interval = MONTH_INTERVAL ;

	f[356].name     = "strmnr_2_3" ;
	f[356].descrip  = "strmnr_2_3" ;
	f[356].units    = "" ;
	f[356].type     = NC_FLOAT ;
	f[356].scale    = 1.000000 ;
	f[356].interval = MONTH_INTERVAL ;

	f[357].name     = "s1mnr_1_1" ;
	f[357].descrip  = "s1mnr_1_1" ;
	f[357].units    = "" ;
	f[357].type     = NC_FLOAT ;
	f[357].scale    = 1.000000 ;
	f[357].interval = MONTH_INTERVAL ;

	f[358].name     = "s1mnr_2_1" ;
	f[358].descrip  = "s1mnr_2_1" ;
	f[358].units    = "" ;
	f[358].type     = NC_FLOAT ;
	f[358].scale    = 1.000000 ;
	f[358].interval = MONTH_INTERVAL ;

	f[359].name     = "s1mnr_1_2" ;
	f[359].descrip  = "s1mnr_1_2" ;
	f[359].units    = "" ;
	f[359].type     = NC_FLOAT ;
	f[359].scale    = 1.000000 ;
	f[359].interval = MONTH_INTERVAL ;

	f[360].name     = "s1mnr_2_2" ;
	f[360].descrip  = "s1mnr_2_2" ;
	f[360].units    = "" ;
	f[360].type     = NC_FLOAT ;
	f[360].scale    = 1.000000 ;
	f[360].interval = MONTH_INTERVAL ;

	f[361].name     = "s1mnr_1_3" ;
	f[361].descrip  = "s1mnr_1_3" ;
	f[361].units    = "" ;
	f[361].type     = NC_FLOAT ;
	f[361].scale    = 1.000000 ;
	f[361].interval = MONTH_INTERVAL ;

	f[362].name     = "s1mnr_2_3" ;
	f[362].descrip  = "s1mnr_2_3" ;
	f[362].units    = "" ;
	f[362].type     = NC_FLOAT ;
	f[362].scale    = 1.000000 ;
	f[362].interval = MONTH_INTERVAL ;

	f[363].name     = "s2mnr_1" ;
	f[363].descrip  = "s2mnr_1" ;
	f[363].units    = "" ;
	f[363].type     = NC_FLOAT ;
	f[363].scale    = 1.000000 ;
	f[363].interval = MONTH_INTERVAL ;

	f[364].name     = "s2mnr_2" ;
	f[364].descrip  = "s2mnr_2" ;
	f[364].units    = "" ;
	f[364].type     = NC_FLOAT ;
	f[364].scale    = 1.000000 ;
	f[364].interval = MONTH_INTERVAL ;

	f[365].name     = "s2mnr_3" ;
	f[365].descrip  = "s2mnr_3" ;
	f[365].units    = "" ;
	f[365].type     = NC_FLOAT ;
	f[365].scale    = 1.000000 ;
	f[365].interval = MONTH_INTERVAL ;

	f[366].name     = "s3mnr_1" ;
	f[366].descrip  = "s3mnr_1" ;
	f[366].units    = "" ;
	f[366].type     = NC_FLOAT ;
	f[366].scale    = 1.000000 ;
	f[366].interval = MONTH_INTERVAL ;

	f[367].name     = "s3mnr_2" ;
	f[367].descrip  = "s3mnr_2" ;
	f[367].units    = "" ;
	f[367].type     = NC_FLOAT ;
	f[367].scale    = 1.000000 ;
	f[367].interval = MONTH_INTERVAL ;

	f[368].name     = "s3mnr_3" ;
	f[368].descrip  = "s3mnr_3" ;
	f[368].units    = "" ;
	f[368].type     = NC_FLOAT ;
	f[368].scale    = 1.000000 ;
	f[368].interval = MONTH_INTERVAL ;

	f[369].name     = "tcerat_1" ;
	f[369].descrip  = "tcerat_1" ;
	f[369].units    = "" ;
	f[369].type     = NC_FLOAT ;
	f[369].scale    = 1.000000 ;
	f[369].interval = MONTH_INTERVAL ;

	f[370].name     = "tcerat_2" ;
	f[370].descrip  = "tcerat_2" ;
	f[370].units    = "" ;
	f[370].type     = NC_FLOAT ;
	f[370].scale    = 1.000000 ;
	f[370].interval = MONTH_INTERVAL ;

	f[371].name     = "tcerat_3" ;
	f[371].descrip  = "tcerat_3" ;
	f[371].units    = "" ;
	f[371].type     = NC_FLOAT ;
	f[371].scale    = 1.000000 ;
	f[371].interval = MONTH_INTERVAL ;

	f[372].name     = "tcnpro" ;
	f[372].descrip  = "tcnpro" ;
	f[372].units    = "" ;
	f[372].type     = NC_FLOAT ;
	f[372].scale    = 1.000000 ;
	f[372].interval = MONTH_INTERVAL ;

	f[373].name     = "tomres_1" ;
	f[373].descrip  = "tomres_1" ;
	f[373].units    = "" ;
	f[373].type     = NC_FLOAT ;
	f[373].scale    = 1.000000 ;
	f[373].interval = MONTH_INTERVAL ;

	f[374].name     = "tomres_2" ;
	f[374].descrip  = "tomres_2" ;
	f[374].units    = "" ;
	f[374].type     = NC_FLOAT ;
	f[374].scale    = 1.000000 ;
	f[374].interval = MONTH_INTERVAL ;

	f[375].name     = "totalc" ;
	f[375].descrip  = "totalc" ;
	f[375].units    = "" ;
	f[375].type     = NC_FLOAT ;
	f[375].scale    = 1.000000 ;
	f[375].interval = MONTH_INTERVAL ;

	f[376].name     = "totale_1" ;
	f[376].descrip  = "totale_1" ;
	f[376].units    = "" ;
	f[376].type     = NC_FLOAT ;
	f[376].scale    = 1.000000 ;
	f[376].interval = MONTH_INTERVAL ;

	f[377].name     = "totale_2" ;
	f[377].descrip  = "totale_2" ;
	f[377].units    = "" ;
	f[377].type     = NC_FLOAT ;
	f[377].scale    = 1.000000 ;
	f[377].interval = MONTH_INTERVAL ;

	f[378].name     = "totale_3" ;
	f[378].descrip  = "totale_3" ;
	f[378].units    = "" ;
	f[378].type     = NC_FLOAT ;
	f[378].scale    = 1.000000 ;
	f[378].interval = MONTH_INTERVAL ;

	f[379].name     = "volex" ;
	f[379].descrip  = "volex" ;
	f[379].units    = "" ;
	f[379].type     = NC_FLOAT ;
	f[379].scale    = 1.000000 ;
	f[379].interval = MONTH_INTERVAL ;

	f[380].name     = "volgm" ;
	f[380].descrip  = "volgm" ;
	f[380].units    = "" ;
	f[380].type     = NC_FLOAT ;
	f[380].scale    = 1.000000 ;
	f[380].interval = MONTH_INTERVAL ;

	f[381].name     = "volpl" ;
	f[381].descrip  = "volpl" ;
	f[381].units    = "" ;
	f[381].type     = NC_FLOAT ;
	f[381].scale    = 1.000000 ;
	f[381].interval = MONTH_INTERVAL ;

	f[382].name     = "wdfx" ;
	f[382].descrip  = "wdfx" ;
	f[382].units    = "" ;
	f[382].type     = NC_FLOAT ;
	f[382].scale    = 1.000000 ;
	f[382].interval = MONTH_INTERVAL ;

	f[383].name     = "wdfxa" ;
	f[383].descrip  = "wdfxa" ;
	f[383].units    = "" ;
	f[383].type     = NC_FLOAT ;
	f[383].scale    = 1.000000 ;
	f[383].interval = MONTH_INTERVAL ;

	f[384].name     = "wdfxma" ;
	f[384].descrip  = "wdfxma" ;
	f[384].units    = "" ;
	f[384].type     = NC_FLOAT ;
	f[384].scale    = 1.000000 ;
	f[384].interval = MONTH_INTERVAL ;

	f[385].name     = "wdfxms" ;
	f[385].descrip  = "wdfxms" ;
	f[385].units    = "" ;
	f[385].type     = NC_FLOAT ;
	f[385].scale    = 1.000000 ;
	f[385].interval = MONTH_INTERVAL ;

	f[386].name     = "wdfxs" ;
	f[386].descrip  = "wdfxs" ;
	f[386].units    = "" ;
	f[386].type     = NC_FLOAT ;
	f[386].scale    = 1.000000 ;
	f[386].interval = MONTH_INTERVAL ;

	f[387].name     = "acrcis_1" ;
	f[387].descrip  = "acrcis_1" ;
	f[387].units    = "" ;
	f[387].type     = NC_FLOAT ;
	f[387].scale    = 1.000000 ;
	f[387].interval = MONTH_INTERVAL ;

	f[388].name     = "acrcis_2" ;
	f[388].descrip  = "acrcis_2" ;
	f[388].units    = "" ;
	f[388].type     = NC_FLOAT ;
	f[388].scale    = 1.000000 ;
	f[388].interval = MONTH_INTERVAL ;

	f[389].name     = "afbcis_1" ;
	f[389].descrip  = "afbcis_1" ;
	f[389].units    = "" ;
	f[389].type     = NC_FLOAT ;
	f[389].scale    = 1.000000 ;
	f[389].interval = MONTH_INTERVAL ;

	f[390].name     = "afbcis_2" ;
	f[390].descrip  = "afbcis_2" ;
	f[390].units    = "" ;
	f[390].type     = NC_FLOAT ;
	f[390].scale    = 1.000000 ;
	f[390].interval = MONTH_INTERVAL ;

	f[391].name     = "afrcis_1" ;
	f[391].descrip  = "afrcis_1" ;
	f[391].units    = "" ;
	f[391].type     = NC_FLOAT ;
	f[391].scale    = 1.000000 ;
	f[391].interval = MONTH_INTERVAL ;

	f[392].name     = "afrcis_2" ;
	f[392].descrip  = "afrcis_2" ;
	f[392].units    = "" ;
	f[392].type     = NC_FLOAT ;
	f[392].scale    = 1.000000 ;
	f[392].interval = MONTH_INTERVAL ;

	f[393].name     = "alvcis_1" ;
	f[393].descrip  = "alvcis_1" ;
	f[393].units    = "" ;
	f[393].type     = NC_FLOAT ;
	f[393].scale    = 1.000000 ;
	f[393].interval = MONTH_INTERVAL ;

	f[394].name     = "alvcis_2" ;
	f[394].descrip  = "alvcis_2" ;
	f[394].units    = "" ;
	f[394].type     = NC_FLOAT ;
	f[394].scale    = 1.000000 ;
	f[394].interval = MONTH_INTERVAL ;

	f[395].name     = "alwcis_1" ;
	f[395].descrip  = "alwcis_1" ;
	f[395].units    = "" ;
	f[395].type     = NC_FLOAT ;
	f[395].scale    = 1.000000 ;
	f[395].interval = MONTH_INTERVAL ;

	f[396].name     = "alwcis_2" ;
	f[396].descrip  = "alwcis_2" ;
	f[396].units    = "" ;
	f[396].type     = NC_FLOAT ;
	f[396].scale    = 1.000000 ;
	f[396].interval = MONTH_INTERVAL ;

	f[397].name     = "crootc" ;
	f[397].descrip  = "crootc" ;
	f[397].units    = "" ;
	f[397].type     = NC_FLOAT ;
	f[397].scale    = 1.000000 ;
	f[397].interval = MONTH_INTERVAL ;

	f[398].name     = "croote_1" ;
	f[398].descrip  = "croote_1" ;
	f[398].units    = "" ;
	f[398].type     = NC_FLOAT ;
	f[398].scale    = 1.000000 ;
	f[398].interval = MONTH_INTERVAL ;

	f[399].name     = "croote_2" ;
	f[399].descrip  = "croote_2" ;
	f[399].units    = "" ;
	f[399].type     = NC_FLOAT ;
	f[399].scale    = 1.000000 ;
	f[399].interval = MONTH_INTERVAL ;

	f[400].name     = "croote_3" ;
	f[400].descrip  = "croote_3" ;
	f[400].units    = "" ;
	f[400].type     = NC_FLOAT ;
	f[400].scale    = 1.000000 ;
	f[400].interval = MONTH_INTERVAL ;

	f[401].name     = "crtacc" ;
	f[401].descrip  = "crtacc" ;
	f[401].units    = "" ;
	f[401].type     = NC_FLOAT ;
	f[401].scale    = 1.000000 ;
	f[401].interval = MONTH_INTERVAL ;

	f[402].name     = "crtcis_1" ;
	f[402].descrip  = "crtcis_1" ;
	f[402].units    = "" ;
	f[402].type     = NC_FLOAT ;
	f[402].scale    = 1.000000 ;
	f[402].interval = MONTH_INTERVAL ;

	f[403].name     = "crtcis_2" ;
	f[403].descrip  = "crtcis_2" ;
	f[403].units    = "" ;
	f[403].type     = NC_FLOAT ;
	f[403].scale    = 1.000000 ;
	f[403].interval = MONTH_INTERVAL ;

	f[404].name     = "eupprt_1_1" ;
	f[404].descrip  = "eupprt_1_1" ;
	f[404].units    = "" ;
	f[404].type     = NC_FLOAT ;
	f[404].scale    = 1.000000 ;
	f[404].interval = MONTH_INTERVAL ;

	f[405].name     = "eupprt_2_1" ;
	f[405].descrip  = "eupprt_2_1" ;
	f[405].units    = "" ;
	f[405].type     = NC_FLOAT ;
	f[405].scale    = 1.000000 ;
	f[405].interval = MONTH_INTERVAL ;

	f[406].name     = "eupprt_3_1" ;
	f[406].descrip  = "eupprt_3_1" ;
	f[406].units    = "" ;
	f[406].type     = NC_FLOAT ;
	f[406].scale    = 1.000000 ;
	f[406].interval = MONTH_INTERVAL ;

	f[407].name     = "eupprt_4_1" ;
	f[407].descrip  = "eupprt_4_1" ;
	f[407].units    = "" ;
	f[407].type     = NC_FLOAT ;
	f[407].scale    = 1.000000 ;
	f[407].interval = MONTH_INTERVAL ;

	f[408].name     = "eupprt_5_1" ;
	f[408].descrip  = "eupprt_5_1" ;
	f[408].units    = "" ;
	f[408].type     = NC_FLOAT ;
	f[408].scale    = 1.000000 ;
	f[408].interval = MONTH_INTERVAL ;

	f[409].name     = "eupprt_1_2" ;
	f[409].descrip  = "eupprt_1_2" ;
	f[409].units    = "" ;
	f[409].type     = NC_FLOAT ;
	f[409].scale    = 1.000000 ;
	f[409].interval = MONTH_INTERVAL ;

	f[410].name     = "eupprt_2_2" ;
	f[410].descrip  = "eupprt_2_2" ;
	f[410].units    = "" ;
	f[410].type     = NC_FLOAT ;
	f[410].scale    = 1.000000 ;
	f[410].interval = MONTH_INTERVAL ;

	f[411].name     = "eupprt_3_2" ;
	f[411].descrip  = "eupprt_3_2" ;
	f[411].units    = "" ;
	f[411].type     = NC_FLOAT ;
	f[411].scale    = 1.000000 ;
	f[411].interval = MONTH_INTERVAL ;

	f[412].name     = "eupprt_4_2" ;
	f[412].descrip  = "eupprt_4_2" ;
	f[412].units    = "" ;
	f[412].type     = NC_FLOAT ;
	f[412].scale    = 1.000000 ;
	f[412].interval = MONTH_INTERVAL ;

	f[413].name     = "eupprt_5_2" ;
	f[413].descrip  = "eupprt_5_2" ;
	f[413].units    = "" ;
	f[413].type     = NC_FLOAT ;
	f[413].scale    = 1.000000 ;
	f[413].interval = MONTH_INTERVAL ;

	f[414].name     = "eupprt_1_3" ;
	f[414].descrip  = "eupprt_1_3" ;
	f[414].units    = "" ;
	f[414].type     = NC_FLOAT ;
	f[414].scale    = 1.000000 ;
	f[414].interval = MONTH_INTERVAL ;

	f[415].name     = "eupprt_2_3" ;
	f[415].descrip  = "eupprt_2_3" ;
	f[415].units    = "" ;
	f[415].type     = NC_FLOAT ;
	f[415].scale    = 1.000000 ;
	f[415].interval = MONTH_INTERVAL ;

	f[416].name     = "eupprt_3_3" ;
	f[416].descrip  = "eupprt_3_3" ;
	f[416].units    = "" ;
	f[416].type     = NC_FLOAT ;
	f[416].scale    = 1.000000 ;
	f[416].interval = MONTH_INTERVAL ;

	f[417].name     = "eupprt_4_3" ;
	f[417].descrip  = "eupprt_4_3" ;
	f[417].units    = "" ;
	f[417].type     = NC_FLOAT ;
	f[417].scale    = 1.000000 ;
	f[417].interval = MONTH_INTERVAL ;

	f[418].name     = "eupprt_5_3" ;
	f[418].descrip  = "eupprt_5_3" ;
	f[418].units    = "" ;
	f[418].type     = NC_FLOAT ;
	f[418].scale    = 1.000000 ;
	f[418].interval = MONTH_INTERVAL ;

	f[419].name     = "fbrchc" ;
	f[419].descrip  = "fbrchc" ;
	f[419].units    = "" ;
	f[419].type     = NC_FLOAT ;
	f[419].scale    = 1.000000 ;
	f[419].interval = MONTH_INTERVAL ;

	f[420].name     = "fbrche_1" ;
	f[420].descrip  = "fbrche_1" ;
	f[420].units    = "" ;
	f[420].type     = NC_FLOAT ;
	f[420].scale    = 1.000000 ;
	f[420].interval = MONTH_INTERVAL ;

	f[421].name     = "fbrche_2" ;
	f[421].descrip  = "fbrche_2" ;
	f[421].units    = "" ;
	f[421].type     = NC_FLOAT ;
	f[421].scale    = 1.000000 ;
	f[421].interval = MONTH_INTERVAL ;

	f[422].name     = "fbrche_3" ;
	f[422].descrip  = "fbrche_3" ;
	f[422].units    = "" ;
	f[422].type     = NC_FLOAT ;
	f[422].scale    = 1.000000 ;
	f[422].interval = MONTH_INTERVAL ;

	f[423].name     = "fbracc" ;
	f[423].descrip  = "fbracc" ;
	f[423].units    = "" ;
	f[423].type     = NC_FLOAT ;
	f[423].scale    = 1.000000 ;
	f[423].interval = MONTH_INTERVAL ;

	f[424].name     = "fbrcis_1" ;
	f[424].descrip  = "fbrcis_1" ;
	f[424].units    = "" ;
	f[424].type     = NC_FLOAT ;
	f[424].scale    = 1.000000 ;
	f[424].interval = MONTH_INTERVAL ;

	f[425].name     = "fbrcis_2" ;
	f[425].descrip  = "fbrcis_2" ;
	f[425].units    = "" ;
	f[425].type     = NC_FLOAT ;
	f[425].scale    = 1.000000 ;
	f[425].interval = MONTH_INTERVAL ;

	f[426].name     = "fcacc" ;
	f[426].descrip  = "fcacc" ;
	f[426].units    = "" ;
	f[426].type     = NC_FLOAT ;
	f[426].scale    = 1.000000 ;
	f[426].interval = MONTH_INTERVAL ;

	f[427].name     = "forstg_1" ;
	f[427].descrip  = "forstg_1" ;
	f[427].units    = "" ;
	f[427].type     = NC_FLOAT ;
	f[427].scale    = 1.000000 ;
	f[427].interval = MONTH_INTERVAL ;

	f[428].name     = "forstg_2" ;
	f[428].descrip  = "forstg_2" ;
	f[428].units    = "" ;
	f[428].type     = NC_FLOAT ;
	f[428].scale    = 1.000000 ;
	f[428].interval = MONTH_INTERVAL ;

	f[429].name     = "forstg_3" ;
	f[429].descrip  = "forstg_3" ;
	f[429].units    = "" ;
	f[429].type     = NC_FLOAT ;
	f[429].scale    = 1.000000 ;
	f[429].interval = MONTH_INTERVAL ;

	f[430].name     = "frootc" ;
	f[430].descrip  = "frootc" ;
	f[430].units    = "" ;
	f[430].type     = NC_FLOAT ;
	f[430].scale    = 1.000000 ;
	f[430].interval = MONTH_INTERVAL ;

	f[431].name     = "froote_1" ;
	f[431].descrip  = "froote_1" ;
	f[431].units    = "" ;
	f[431].type     = NC_FLOAT ;
	f[431].scale    = 1.000000 ;
	f[431].interval = MONTH_INTERVAL ;

	f[432].name     = "froote_2" ;
	f[432].descrip  = "froote_2" ;
	f[432].units    = "" ;
	f[432].type     = NC_FLOAT ;
	f[432].scale    = 1.000000 ;
	f[432].interval = MONTH_INTERVAL ;

	f[433].name     = "froote_3" ;
	f[433].descrip  = "froote_3" ;
	f[433].units    = "" ;
	f[433].type     = NC_FLOAT ;
	f[433].scale    = 1.000000 ;
	f[433].interval = MONTH_INTERVAL ;

	f[434].name     = "frtacc" ;
	f[434].descrip  = "frtacc" ;
	f[434].units    = "" ;
	f[434].type     = NC_FLOAT ;
	f[434].scale    = 1.000000 ;
	f[434].interval = MONTH_INTERVAL ;

	f[435].name     = "frtcis_1" ;
	f[435].descrip  = "frtcis_1" ;
	f[435].units    = "" ;
	f[435].type     = NC_FLOAT ;
	f[435].scale    = 1.000000 ;
	f[435].interval = MONTH_INTERVAL ;

	f[436].name     = "frtcis_2" ;
	f[436].descrip  = "frtcis_2" ;
	f[436].units    = "" ;
	f[436].type     = NC_FLOAT ;
	f[436].scale    = 1.000000 ;
	f[436].interval = MONTH_INTERVAL ;

	f[437].name     = "frstc" ;
	f[437].descrip  = "frstc" ;
	f[437].units    = "" ;
	f[437].type     = NC_FLOAT ;
	f[437].scale    = 1.000000 ;
	f[437].interval = MONTH_INTERVAL ;

	f[438].name     = "frste_1" ;
	f[438].descrip  = "frste_1" ;
	f[438].units    = "" ;
	f[438].type     = NC_FLOAT ;
	f[438].scale    = 1.000000 ;
	f[438].interval = MONTH_INTERVAL ;

	f[439].name     = "frste_2" ;
	f[439].descrip  = "frste_2" ;
	f[439].units    = "" ;
	f[439].type     = NC_FLOAT ;
	f[439].scale    = 1.000000 ;
	f[439].interval = MONTH_INTERVAL ;

	f[440].name     = "frste_3" ;
	f[440].descrip  = "frste_3" ;
	f[440].units    = "" ;
	f[440].type     = NC_FLOAT ;
	f[440].scale    = 1.000000 ;
	f[440].interval = MONTH_INTERVAL ;

	f[441].name     = "fsysc" ;
	f[441].descrip  = "fsysc" ;
	f[441].units    = "" ;
	f[441].type     = NC_FLOAT ;
	f[441].scale    = 1.000000 ;
	f[441].interval = MONTH_INTERVAL ;

	f[442].name     = "fsyse_1" ;
	f[442].descrip  = "fsyse_1" ;
	f[442].units    = "" ;
	f[442].type     = NC_FLOAT ;
	f[442].scale    = 1.000000 ;
	f[442].interval = MONTH_INTERVAL ;

	f[443].name     = "fsyse_2" ;
	f[443].descrip  = "fsyse_2" ;
	f[443].units    = "" ;
	f[443].type     = NC_FLOAT ;
	f[443].scale    = 1.000000 ;
	f[443].interval = MONTH_INTERVAL ;

	f[444].name     = "fsyse_3" ;
	f[444].descrip  = "fsyse_3" ;
	f[444].units    = "" ;
	f[444].type     = NC_FLOAT ;
	f[444].scale    = 1.000000 ;
	f[444].interval = MONTH_INTERVAL ;

	f[445].name     = "rleavc" ;
	f[445].descrip  = "rleavc" ;
	f[445].units    = "" ;
	f[445].type     = NC_FLOAT ;
	f[445].scale    = 1.000000 ;
	f[445].interval = MONTH_INTERVAL ;

	f[446].name     = "rleave_1" ;
	f[446].descrip  = "rleave_1" ;
	f[446].units    = "" ;
	f[446].type     = NC_FLOAT ;
	f[446].scale    = 1.000000 ;
	f[446].interval = MONTH_INTERVAL ;

	f[447].name     = "rleave_2" ;
	f[447].descrip  = "rleave_2" ;
	f[447].units    = "" ;
	f[447].type     = NC_FLOAT ;
	f[447].scale    = 1.000000 ;
	f[447].interval = MONTH_INTERVAL ;

	f[448].name     = "rleave_3" ;
	f[448].descrip  = "rleave_3" ;
	f[448].units    = "" ;
	f[448].type     = NC_FLOAT ;
	f[448].scale    = 1.000000 ;
	f[448].interval = MONTH_INTERVAL ;

	f[449].name     = "rlvacc" ;
	f[449].descrip  = "rlvacc" ;
	f[449].units    = "" ;
	f[449].type     = NC_FLOAT ;
	f[449].scale    = 1.000000 ;
	f[449].interval = MONTH_INTERVAL ;

	f[450].name     = "rlvcis_1" ;
	f[450].descrip  = "rlvcis_1" ;
	f[450].units    = "" ;
	f[450].type     = NC_FLOAT ;
	f[450].scale    = 1.000000 ;
	f[450].interval = MONTH_INTERVAL ;

	f[451].name     = "rlvcis_2" ;
	f[451].descrip  = "rlvcis_2" ;
	f[451].units    = "" ;
	f[451].type     = NC_FLOAT ;
	f[451].scale    = 1.000000 ;
	f[451].interval = MONTH_INTERVAL ;

	f[452].name     = "rlwodc" ;
	f[452].descrip  = "rlwodc" ;
	f[452].units    = "" ;
	f[452].type     = NC_FLOAT ;
	f[452].scale    = 1.000000 ;
	f[452].interval = MONTH_INTERVAL ;

	f[453].name     = "rlwode_1" ;
	f[453].descrip  = "rlwode_1" ;
	f[453].units    = "" ;
	f[453].type     = NC_FLOAT ;
	f[453].scale    = 1.000000 ;
	f[453].interval = MONTH_INTERVAL ;

	f[454].name     = "rlwode_2" ;
	f[454].descrip  = "rlwode_2" ;
	f[454].units    = "" ;
	f[454].type     = NC_FLOAT ;
	f[454].scale    = 1.000000 ;
	f[454].interval = MONTH_INTERVAL ;

	f[455].name     = "rlwode_3" ;
	f[455].descrip  = "rlwode_3" ;
	f[455].units    = "" ;
	f[455].type     = NC_FLOAT ;
	f[455].scale    = 1.000000 ;
	f[455].interval = MONTH_INTERVAL ;

	f[456].name     = "rlwacc" ;
	f[456].descrip  = "rlwacc" ;
	f[456].units    = "" ;
	f[456].type     = NC_FLOAT ;
	f[456].scale    = 1.000000 ;
	f[456].interval = MONTH_INTERVAL ;

	f[457].name     = "rlwcis_1" ;
	f[457].descrip  = "rlwcis_1" ;
	f[457].units    = "" ;
	f[457].type     = NC_FLOAT ;
	f[457].scale    = 1.000000 ;
	f[457].interval = MONTH_INTERVAL ;

	f[458].name     = "rlwcis_2" ;
	f[458].descrip  = "rlwcis_2" ;
	f[458].units    = "" ;
	f[458].type     = NC_FLOAT ;
	f[458].scale    = 1.000000 ;
	f[458].interval = MONTH_INTERVAL ;

	f[459].name     = "sumrsp" ;
	f[459].descrip  = "sumrsp" ;
	f[459].units    = "" ;
	f[459].type     = NC_FLOAT ;
	f[459].scale    = 1.000000 ;
	f[459].interval = MONTH_INTERVAL ;

	f[460].name     = "tcrem" ;
	f[460].descrip  = "tcrem" ;
	f[460].units    = "" ;
	f[460].type     = NC_FLOAT ;
	f[460].scale    = 1.000000 ;
	f[460].interval = MONTH_INTERVAL ;

	f[461].name     = "terem_1" ;
	f[461].descrip  = "terem_1" ;
	f[461].units    = "" ;
	f[461].type     = NC_FLOAT ;
	f[461].scale    = 1.000000 ;
	f[461].interval = MONTH_INTERVAL ;

	f[462].name     = "terem_2" ;
	f[462].descrip  = "terem_2" ;
	f[462].units    = "" ;
	f[462].type     = NC_FLOAT ;
	f[462].scale    = 1.000000 ;
	f[462].interval = MONTH_INTERVAL ;

	f[463].name     = "terem_3" ;
	f[463].descrip  = "terem_3" ;
	f[463].units    = "" ;
	f[463].type     = NC_FLOAT ;
	f[463].scale    = 1.000000 ;
	f[463].interval = MONTH_INTERVAL ;

	f[464].name     = "w1lig" ;
	f[464].descrip  = "w1lig" ;
	f[464].units    = "" ;
	f[464].type     = NC_FLOAT ;
	f[464].scale    = 1.000000 ;
	f[464].interval = MONTH_INTERVAL ;

	f[465].name     = "w2lig" ;
	f[465].descrip  = "w2lig" ;
	f[465].units    = "" ;
	f[465].type     = NC_FLOAT ;
	f[465].scale    = 1.000000 ;
	f[465].interval = MONTH_INTERVAL ;

	f[466].name     = "w3lig" ;
	f[466].descrip  = "w3lig" ;
	f[466].units    = "" ;
	f[466].type     = NC_FLOAT ;
	f[466].scale    = 1.000000 ;
	f[466].interval = MONTH_INTERVAL ;

	f[467].name     = "w1mnr_1" ;
	f[467].descrip  = "w1mnr_1" ;
	f[467].units    = "" ;
	f[467].type     = NC_FLOAT ;
	f[467].scale    = 1.000000 ;
	f[467].interval = MONTH_INTERVAL ;

	f[468].name     = "w1mnr_2" ;
	f[468].descrip  = "w1mnr_2" ;
	f[468].units    = "" ;
	f[468].type     = NC_FLOAT ;
	f[468].scale    = 1.000000 ;
	f[468].interval = MONTH_INTERVAL ;

	f[469].name     = "w1mnr_3" ;
	f[469].descrip  = "w1mnr_3" ;
	f[469].units    = "" ;
	f[469].type     = NC_FLOAT ;
	f[469].scale    = 1.000000 ;
	f[469].interval = MONTH_INTERVAL ;

	f[470].name     = "w2mnr_1" ;
	f[470].descrip  = "w2mnr_1" ;
	f[470].units    = "" ;
	f[470].type     = NC_FLOAT ;
	f[470].scale    = 1.000000 ;
	f[470].interval = MONTH_INTERVAL ;

	f[471].name     = "w2mnr_2" ;
	f[471].descrip  = "w2mnr_2" ;
	f[471].units    = "" ;
	f[471].type     = NC_FLOAT ;
	f[471].scale    = 1.000000 ;
	f[471].interval = MONTH_INTERVAL ;

	f[472].name     = "w2mnr_3" ;
	f[472].descrip  = "w2mnr_3" ;
	f[472].units    = "" ;
	f[472].type     = NC_FLOAT ;
	f[472].scale    = 1.000000 ;
	f[472].interval = MONTH_INTERVAL ;

	f[473].name     = "w3mnr_1" ;
	f[473].descrip  = "w3mnr_1" ;
	f[473].units    = "" ;
	f[473].type     = NC_FLOAT ;
	f[473].scale    = 1.000000 ;
	f[473].interval = MONTH_INTERVAL ;

	f[474].name     = "w3mnr_2" ;
	f[474].descrip  = "w3mnr_2" ;
	f[474].units    = "" ;
	f[474].type     = NC_FLOAT ;
	f[474].scale    = 1.000000 ;
	f[474].interval = MONTH_INTERVAL ;

	f[475].name     = "w3mnr_3" ;
	f[475].descrip  = "w3mnr_3" ;
	f[475].units    = "" ;
	f[475].type     = NC_FLOAT ;
	f[475].scale    = 1.000000 ;
	f[475].interval = MONTH_INTERVAL ;

	f[476].name     = "wd1cis_1" ;
	f[476].descrip  = "wd1cis_1" ;
	f[476].units    = "" ;
	f[476].type     = NC_FLOAT ;
	f[476].scale    = 1.000000 ;
	f[476].interval = MONTH_INTERVAL ;

	f[477].name     = "wd1cis_2" ;
	f[477].descrip  = "wd1cis_2" ;
	f[477].units    = "" ;
	f[477].type     = NC_FLOAT ;
	f[477].scale    = 1.000000 ;
	f[477].interval = MONTH_INTERVAL ;

	f[478].name     = "wd2cis_1" ;
	f[478].descrip  = "wd2cis_1" ;
	f[478].units    = "" ;
	f[478].type     = NC_FLOAT ;
	f[478].scale    = 1.000000 ;
	f[478].interval = MONTH_INTERVAL ;

	f[479].name     = "wd2cis_2" ;
	f[479].descrip  = "wd2cis_2" ;
	f[479].units    = "" ;
	f[479].type     = NC_FLOAT ;
	f[479].scale    = 1.000000 ;
	f[479].interval = MONTH_INTERVAL ;

	f[480].name     = "wd3cis_1" ;
	f[480].descrip  = "wd3cis_1" ;
	f[480].units    = "" ;
	f[480].type     = NC_FLOAT ;
	f[480].scale    = 1.000000 ;
	f[480].interval = MONTH_INTERVAL ;

	f[481].name     = "wd3cis_2" ;
	f[481].descrip  = "wd3cis_2" ;
	f[481].units    = "" ;
	f[481].type     = NC_FLOAT ;
	f[481].scale    = 1.000000 ;
	f[481].interval = MONTH_INTERVAL ;

	f[482].name     = "wood1c" ;
	f[482].descrip  = "wood1c" ;
	f[482].units    = "" ;
	f[482].type     = NC_FLOAT ;
	f[482].scale    = 1.000000 ;
	f[482].interval = MONTH_INTERVAL ;

	f[483].name     = "wood2c" ;
	f[483].descrip  = "wood2c" ;
	f[483].units    = "" ;
	f[483].type     = NC_FLOAT ;
	f[483].scale    = 1.000000 ;
	f[483].interval = MONTH_INTERVAL ;

	f[484].name     = "wood3c" ;
	f[484].descrip  = "wood3c" ;
	f[484].units    = "" ;
	f[484].type     = NC_FLOAT ;
	f[484].scale    = 1.000000 ;
	f[484].interval = MONTH_INTERVAL ;

	f[485].name     = "woodc" ;
	f[485].descrip  = "woodc" ;
	f[485].units    = "" ;
	f[485].type     = NC_FLOAT ;
	f[485].scale    = 1.000000 ;
	f[485].interval = MONTH_INTERVAL ;

	f[486].name     = "wood1e_1" ;
	f[486].descrip  = "wood1e_1" ;
	f[486].units    = "" ;
	f[486].type     = NC_FLOAT ;
	f[486].scale    = 1.000000 ;
	f[486].interval = MONTH_INTERVAL ;

	f[487].name     = "wood1e_2" ;
	f[487].descrip  = "wood1e_2" ;
	f[487].units    = "" ;
	f[487].type     = NC_FLOAT ;
	f[487].scale    = 1.000000 ;
	f[487].interval = MONTH_INTERVAL ;

	f[488].name     = "wood1e_3" ;
	f[488].descrip  = "wood1e_3" ;
	f[488].units    = "" ;
	f[488].type     = NC_FLOAT ;
	f[488].scale    = 1.000000 ;
	f[488].interval = MONTH_INTERVAL ;

	f[489].name     = "wood2e_1" ;
	f[489].descrip  = "wood2e_1" ;
	f[489].units    = "" ;
	f[489].type     = NC_FLOAT ;
	f[489].scale    = 1.000000 ;
	f[489].interval = MONTH_INTERVAL ;

	f[490].name     = "wood2e_2" ;
	f[490].descrip  = "wood2e_2" ;
	f[490].units    = "" ;
	f[490].type     = NC_FLOAT ;
	f[490].scale    = 1.000000 ;
	f[490].interval = MONTH_INTERVAL ;

	f[491].name     = "wood2e_3" ;
	f[491].descrip  = "wood2e_3" ;
	f[491].units    = "" ;
	f[491].type     = NC_FLOAT ;
	f[491].scale    = 1.000000 ;
	f[491].interval = MONTH_INTERVAL ;

	f[492].name     = "wood3e_1" ;
	f[492].descrip  = "wood3e_1" ;
	f[492].units    = "" ;
	f[492].type     = NC_FLOAT ;
	f[492].scale    = 1.000000 ;
	f[492].interval = MONTH_INTERVAL ;

	f[493].name     = "wood3e_2" ;
	f[493].descrip  = "wood3e_2" ;
	f[493].units    = "" ;
	f[493].type     = NC_FLOAT ;
	f[493].scale    = 1.000000 ;
	f[493].interval = MONTH_INTERVAL ;

	f[494].name     = "wood3e_3" ;
	f[494].descrip  = "wood3e_3" ;
	f[494].units    = "" ;
	f[494].type     = NC_FLOAT ;
	f[494].scale    = 1.000000 ;
	f[494].interval = MONTH_INTERVAL ;

	f[495].name     = "woode_1" ;
	f[495].descrip  = "woode_1" ;
	f[495].units    = "" ;
	f[495].type     = NC_FLOAT ;
	f[495].scale    = 1.000000 ;
	f[495].interval = MONTH_INTERVAL ;

	f[496].name     = "woode_2" ;
	f[496].descrip  = "woode_2" ;
	f[496].units    = "" ;
	f[496].type     = NC_FLOAT ;
	f[496].scale    = 1.000000 ;
	f[496].interval = MONTH_INTERVAL ;

	f[497].name     = "woode_3" ;
	f[497].descrip  = "woode_3" ;
	f[497].units    = "" ;
	f[497].type     = NC_FLOAT ;
	f[497].scale    = 1.000000 ;
	f[497].interval = MONTH_INTERVAL ;

	f[498].name     = "mx_index" ;
	f[498].descrip  = "mx_index" ;
	f[498].units    = "" ;
	f[498].type     = NC_FLOAT ;
	f[498].scale    = 1.000000 ;
	f[498].interval = MONTH_INTERVAL ;

	f[499].name     = "c3c4_index" ;
	f[499].descrip  = "C3 production as a % of total production" ;
	f[499].units    = "%" ;
	f[499].type     = NC_FLOAT ;
	f[499].scale    = 1.000000 ;
	f[499].interval = MONTH_INTERVAL ;

	f[500].name     = "burn" ;
	f[500].descrip  = "burn" ;
	f[500].units    = "" ;
	f[500].type     = NC_FLOAT ;
	f[500].scale    = 1.000000 ;
	f[500].interval = MONTH_INTERVAL ;

	f[501].name     = "CENsurface_runoff" ;
	f[501].descrip  = "surface runoff" ;
	f[501].units    = "mmH2O" ;
	f[501].type     = NC_FLOAT ;
	f[501].scale    = 1.000000 ;
	f[501].interval = MONTH_INTERVAL ;

	f[502].name     = "tmp_index" ;
	f[502].descrip  = "tmp_index" ;
	f[502].units    = "" ;
	f[502].type     = NC_FLOAT ;
	f[502].scale    = 1.000000 ;
	f[502].interval = MONTH_INTERVAL ;

	f[503].name     = "ppt_index" ;
	f[503].descrip  = "ppt_index" ;
	f[503].units    = "" ;
	f[503].type     = NC_FLOAT ;
	f[503].scale    = 1.000000 ;
	f[503].interval = MONTH_INTERVAL ;

	f[504].name     = "symb" ;
	f[504].descrip  = "symb" ;
	f[504].units    = "" ;
	f[504].type     = NC_FLOAT ;
	f[504].scale    = 1.000000 ;
	f[504].interval = MONTH_INTERVAL ;

	f[505].name     = "afcacc" ;
	f[505].descrip  = "afcacc" ;
	f[505].units    = "" ;
	f[505].type     = NC_FLOAT ;
	f[505].scale    = 1.000000 ;
	f[505].interval = MONTH_INTERVAL ;

	f[506].name     = "bgrema" ;
	f[506].descrip  = "bgrema" ;
	f[506].units    = "" ;
	f[506].type     = NC_FLOAT ;
	f[506].scale    = 1.000000 ;
	f[506].interval = MONTH_INTERVAL ;

	f[507].name     = "bgrmai_1" ;
	f[507].descrip  = "bgrmai(1)" ;
	f[507].units    = "" ;
	f[507].type     = NC_FLOAT ;
	f[507].scale    = 1.000000 ;
	f[507].interval = MONTH_INTERVAL ;

	f[508].name     = "bgrmai_2" ;
	f[508].descrip  = "bgrmai(2)" ;
	f[508].units    = "" ;
	f[508].type     = NC_FLOAT ;
	f[508].scale    = 1.000000 ;
	f[508].interval = MONTH_INTERVAL ;

	f[509].name     = "bgrmae_1" ;
	f[509].descrip  = "bgrmae(1)" ;
	f[509].units    = "" ;
	f[509].type     = NC_FLOAT ;
	f[509].scale    = 1.000000 ;
	f[509].interval = MONTH_INTERVAL ;

	f[510].name     = "bgrmae_2" ;
	f[510].descrip  = "bgrmae(2)" ;
	f[510].units    = "" ;
	f[510].type     = NC_FLOAT ;
	f[510].scale    = 1.000000 ;
	f[510].interval = MONTH_INTERVAL ;

	f[511].name     = "bgrmae_3" ;
	f[511].descrip  = "bgrmae(3)" ;
	f[511].units    = "" ;
	f[511].type     = NC_FLOAT ;
	f[511].scale    = 1.000000 ;
	f[511].interval = MONTH_INTERVAL ;

	f[512].name     = "asmos_subsoil" ;
	f[512].descrip  = "soil water content below asmos(nlayer)" ;
	f[512].units    = "cm H2O" ;
	f[512].type     = NC_FLOAT ;
	f[512].scale    = 1.000000 ;
	f[512].interval = MONTH_INTERVAL ;

	// 513-521 spare

	f[522].name     = "tmp_prev_0" ;
	f[522].descrip  = "Jan efolded tmp" ;
	f[522].units    = "" ;
	f[522].type     = NC_FLOAT ;
	f[522].scale    = 1.000000 ;
	f[522].interval = MONTH_INTERVAL ;

	f[523].name     = "tmp_prev_1" ;
	f[523].descrip  = "Feb efolded tmp" ;
	f[523].units    = "" ;
	f[523].type     = NC_FLOAT ;
	f[523].scale    = 1.000000 ;
	f[523].interval = MONTH_INTERVAL ;

	f[524].name     = "tmp_prev_2" ;
	f[524].descrip  = "Mar efolded tmp" ;
	f[524].units    = "" ;
	f[524].type     = NC_FLOAT ;
	f[524].scale    = 1.000000 ;
	f[524].interval = MONTH_INTERVAL ;

	f[525].name     = "tmp_prev_3" ;
	f[525].descrip  = "Apr efolded tmp" ;
	f[525].units    = "" ;
	f[525].type     = NC_FLOAT ;
	f[525].scale    = 1.000000 ;
	f[525].interval = MONTH_INTERVAL ;

	f[526].name     = "tmp_prev_4" ;
	f[526].descrip  = "May efolded tmp" ;
	f[526].units    = "" ;
	f[526].type     = NC_FLOAT ;
	f[526].scale    = 1.000000 ;
	f[526].interval = MONTH_INTERVAL ;

	f[527].name     = "tmp_prev_5" ;
	f[527].descrip  = "Jun efolded tmp" ;
	f[527].units    = "" ;
	f[527].type     = NC_FLOAT ;
	f[527].scale    = 1.000000 ;
	f[527].interval = MONTH_INTERVAL ;

	f[528].name     = "tmp_prev_6" ;
	f[528].descrip  = "Jul efolded tmp" ;
	f[528].units    = "" ;
	f[528].type     = NC_FLOAT ;
	f[528].scale    = 1.000000 ;
	f[528].interval = MONTH_INTERVAL ;

	f[529].name     = "tmp_prev_7" ;
	f[529].descrip  = "Aug efolded tmp" ;
	f[529].units    = "" ;
	f[529].type     = NC_FLOAT ;
	f[529].scale    = 1.000000 ;
	f[529].interval = MONTH_INTERVAL ;

	f[530].name     = "tmp_prev_8" ;
	f[530].descrip  = "Sep efolded tmp" ;
	f[530].units    = "" ;
	f[530].type     = NC_FLOAT ;
	f[530].scale    = 1.000000 ;
	f[530].interval = MONTH_INTERVAL ;

	f[531].name     = "tmp_prev_9" ;
	f[531].descrip  = "Oct efolded tmp" ;
	f[531].units    = "" ;
	f[531].type     = NC_FLOAT ;
	f[531].scale    = 1.000000 ;
	f[531].interval = MONTH_INTERVAL ;

	f[532].name     = "tmp_prev_10" ;
	f[532].descrip  = "Nov efolded tmp" ;
	f[532].units    = "" ;
	f[532].type     = NC_FLOAT ;
	f[532].scale    = 1.000000 ;
	f[532].interval = MONTH_INTERVAL ;

	f[533].name     = "tmp_prev_11" ;
	f[533].descrip  = "Dec efolded tmp" ;
	f[533].units    = "" ;
	f[533].type     = NC_FLOAT ;
	f[533].scale    = 1.000000 ;
	f[533].interval = MONTH_INTERVAL ;

	f[534].name     = "ppt_prev_0" ;
	f[534].descrip  = "Jan efolded ppt" ;
	f[534].units    = "" ;
	f[534].type     = NC_FLOAT ;
	f[534].scale    = 1.000000 ;
	f[534].interval = MONTH_INTERVAL ;

	f[535].name     = "ppt_prev_1" ;
	f[535].descrip  = "Feb efolded ppt" ;
	f[535].units    = "" ;
	f[535].type     = NC_FLOAT ;
	f[535].scale    = 1.000000 ;
	f[535].interval = MONTH_INTERVAL ;

	f[536].name     = "ppt_prev_2" ;
	f[536].descrip  = "Mar efolded ppt" ;
	f[536].units    = "" ;
	f[536].type     = NC_FLOAT ;
	f[536].scale    = 1.000000 ;
	f[536].interval = MONTH_INTERVAL ;

	f[537].name     = "ppt_prev_3" ;
	f[537].descrip  = "Apr efolded ppt" ;
	f[537].units    = "" ;
	f[537].type     = NC_FLOAT ;
	f[537].scale    = 1.000000 ;
	f[537].interval = MONTH_INTERVAL ;

	f[538].name     = "ppt_prev_4" ;
	f[538].descrip  = "May efolded ppt" ;
	f[538].units    = "" ;
	f[538].type     = NC_FLOAT ;
	f[538].scale    = 1.000000 ;
	f[538].interval = MONTH_INTERVAL ;

	f[539].name     = "ppt_prev_5" ;
	f[539].descrip  = "Jun efolded ppt" ;
	f[539].units    = "" ;
	f[539].type     = NC_FLOAT ;
	f[539].scale    = 1.000000 ;
	f[539].interval = MONTH_INTERVAL ;

	f[540].name     = "ppt_prev_6" ;
	f[540].descrip  = "Jul efolded ppt" ;
	f[540].units    = "" ;
	f[540].type     = NC_FLOAT ;
	f[540].scale    = 1.000000 ;
	f[540].interval = MONTH_INTERVAL ;

	f[541].name     = "ppt_prev_7" ;
	f[541].descrip  = "Aug efolded ppt" ;
	f[541].units    = "" ;
	f[541].type     = NC_FLOAT ;
	f[541].scale    = 1.000000 ;
	f[541].interval = MONTH_INTERVAL ;

	f[542].name     = "ppt_prev_8" ;
	f[542].descrip  = "Sep efolded ppt" ;
	f[542].units    = "" ;
	f[542].type     = NC_FLOAT ;
	f[542].scale    = 1.000000 ;
	f[542].interval = MONTH_INTERVAL ;

	f[543].name     = "ppt_prev_9" ;
	f[543].descrip  = "Oct efolded ppt" ;
	f[543].units    = "" ;
	f[543].type     = NC_FLOAT ;
	f[543].scale    = 1.000000 ;
	f[543].interval = MONTH_INTERVAL ;

	f[544].name     = "ppt_prev_10" ;
	f[544].descrip  = "Nov efolded ppt" ;
	f[544].units    = "" ;
	f[544].type     = NC_FLOAT ;
	f[544].scale    = 1.000000 ;
	f[544].interval = MONTH_INTERVAL ;

	f[545].name     = "ppt_prev_11" ;
	f[545].descrip  = "Dec efolded ppt" ;
	f[545].units    = "" ;
	f[545].type     = NC_FLOAT ;
	f[545].scale    = 1.000000 ;
	f[545].interval = MONTH_INTERVAL ;

	f[546].name     = "initial_class_mapss" ;
	f[546].descrip  = "initial class from MAPSS" ;
	f[546].units    = "" ;
	f[546].type     = NC_FLOAT ;
	f[546].scale    = 1.000000 ;
	f[546].interval = MONTH_INTERVAL ;

	f[547].name     = "initial_vclass_mapss" ;
	f[547].descrip  = "initial vclass from MAPSS" ;
	f[547].units    = "" ;
	f[547].type     = NC_FLOAT ;
	f[547].scale    = 1.000000 ;
	f[547].interval = MONTH_INTERVAL ;

	f[548].name     = "efold_t" ;
	f[548].descrip  = " efolding tau" ;
	f[548].units    = "" ;
	f[548].type     = NC_FLOAT ;
	f[548].scale    = 1.000000 ;
	f[548].interval = MONTH_INTERVAL ;

	f[549].name     = "annet" ;
	f[549].descrip  = "year-to-date actual evapotranspiration" ;
	f[549].units    = "cm H2O" ;
	f[549].type     = NC_FLOAT ;
	f[549].scale    = 1.000000 ;
	f[549].interval = MONTH_INTERVAL ;

	f[550].name     = "final_year" ;
	f[550].descrip  = "last year of EQ run" ;
	f[550].units    = "year" ;
	f[550].type     = NC_FLOAT ;
	f[550].scale    = 1.000000 ;
	f[550].interval = MONTH_INTERVAL ;

	for (int i = 0; i<filter_length; i++) f[i].show = strcmp(f[i].name, "")!=0;

} // end of initializeEQoutputList()


