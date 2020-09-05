/* MAPSSbiogeographyModel.cpp */

#define UNIX_MAPSS 1

#ifndef UNIX_MAPSS
#include "StdAfx.h"
#endif

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#ifdef UNIX_MAPSS
#include "netcdf.h"
#else
#include "EnvModel.h"
#endif
#include "ScienceFcns.h"
#include "category_bgc.h"
#include "MAPSSvegClasses.h"

#include "ProcessModel.h"
// #include "ScienceFcns.h"
#include "MAPSSbiogeographyModel.h"

#ifdef UNIX_MAPSS
#include "MC2.h"
bool diagsMAPSS = false;
#else
typedef bool Boolean;
extern bool diagsMAPSS;
extern CString diagMsg;
#endif

#ifndef ASSERT
#define ASSERT assert
#endif

MAPSS_model_parameters::MAPSS_model_parameters()
{
	default_parameters();
} // end of default constructor for MAPSS_model_parameters class


void MAPSS_model_parameters::default_parameters()   
{
	int zone_index, cl, ileaf;

	ice_boundary = 0.f; // growing degree days base 0 C
	frost = 13.0f;	// threshold (C) for beginning, end of growing season
	n_decid_bound = -16.0f; // northern temp threshold, deciduous
	s_decid_bound = 2.25f;  // southern temp threshold, deciduous
	evergreen_productivity_tree_lad = 10.0f; 
	evergreen_productivity_shrub_lad = 10.0f;
	evergreen_gdd = 600.0f; //	growing degree days for evergreen (frost based)
	evergreen_events = 40.f;
	evergreen_events_grow_std = 2.0f; //	(1.5 to 2.0)
	evergreen_events_grow_min = 5.f;
	evergreen_events_grow_cov = 0.30f;
	evergreen_events_grow_rng = 4.0f;	
	evergreen_events_grow = 26.0f;
	evergreen_events_avg = 11.0f;
	evergreen_events_grow_avg = 7.5f; //	(6.0 to 7.5)
	evergreen_gdd_ratio = 100.0f;
	evergreen_gdd_north = 250.0f;
	evergreen_gdd_south = 600.0f;
	forest_threshold = 3.75f; // lai below which forest -> savanna (classif. only)
	LaiUpperBoundsLifeform = 15.f;
	fire_threshold_1 = 2.0f;	/* Threshold for Grass Sum LAI */
	fire_threshold_2 = 2.5f;	/* Threshold for Shrub LAI */
	fire_threshold_5 = 100.f;	/* Threshold for pet */
	fire_threshold_6 = 50.f;	/* Threshold for ppt in High Month */
	max_grass_shrub_threshold = 0.70f;
	short_grass_threshold = 1.15f;
	broad_ppt_min = 40.f; //	minimum monthly ppt for broadleaf, growing season
	spring_grass = 1.5f; // maximum grass lai in 1st month of growing season
	grass_snowpack = 10.f; // max depth of snow expressed as mm of water under which grass will grow 
	pj_xeric_threshold = 2.f;
	pj_max_lai_continental = 2.10f;
	pj_max_lai_maritime = 2.10f;
	north_hard_threshold = 9.00f;
	dry_trop_threshold = 2.f;
	semi_desert_threshold = 0.45f;
	tall_grass_threshold = 2.0f;
	desert_grass_sum_threshold = 1.20f;
	desert_shrub_threshold = 0.175f;
	desert_grass_threshold = 0.10f;
	cool_grass_threshold = 3.00f; 
	max_grass_threshold = 1.50f;
	xeric_savanna_threshold = 0.50f;
	mediterranean_savanna_threshold = 0.75f;
	snow0 = 5.0f; //		temp (C) above which snow fraction equals 0 
	melt_slope = 4.0f; //	temp coefficient for snow melt rate (mm)    
	event_ppt = 0.1f; //	nominal # of rainfall events per unit of rain (i.e. 1/(rain/event))
	event_pet = 50.f; //	pet threshold (mm/mo) for max_events determination
	interc_lai = 3.0f; //	maximum ppt interception per event (mm)
	full_attenuation_lai = 5.0f; // woody lai threshold for light attenuation
	no_attenuation_lai = 0.0f; 
	k_factor_constraint = 170.f;
	tsg_threshold = 0.50f;


	evergreen_productivity_tree = 10.0f; // tree productivity for growing season
	evergreen_productivity_shrub = 10.0f; // used for evergreen/decid decision
	snow1 = -7.0f; //		temp (C) below which snow fraction equals 1
	no_melt = -14.f; //		temp (C) below which no snow melt occurs
	tallgrass_pet_factor = 1.22f;
	tree_pet_factor = 2.55f;

	max_events[0] = 5.0f; //	maximum number of events at pet <= event_pet
	max_events[1] = 10.0f; //   maximum number of events at pet > event_pet

	wp[GRASS] = -1.5f; // permanent wilting point, MPa
	wp[TREE] = -1.5f;
	wp[SHRUB] = -6.0f;

	mapss_subalpine_threshold = 1900.; // growing degree days from 0 C
	wue = 1.0f; // water use efficiency scalar, unitless, nominally 1.0f

	for (zone_index = 0; zone_index<NZONES; zone_index++)
	{
		taiga_tundra_boundary[zone_index] = 1165.f;
		tundra_boundary[zone_index] = 615.0f;
		maritime_boundary[zone_index] = 20.f;
		normal_cond_max[zone_index] = 3.5f;
		LaiUpperBoundsEnergy[zone_index][GRASS] = 6.f;
		LaiUpperBoundsEnergy[zone_index][TREE] = 15.f;
		LaiUpperBoundsEnergy[zone_index][SHRUB] = 3.5f;    
		for (ileaf = NEEDLELEAF; ileaf<=BROADLEAF; ileaf++) 
		{
			min_tree_lai[zone_index][ileaf] = 1.85f; 
			k_factor_slope[zone_index][ileaf] = -2.50f;
		} 
	}  
	taiga_tundra_boundary[mapssBOREAL] = 1330.f;
	tundra_boundary[mapssBOREAL] = 735.f;
	maritime_boundary[mapssTEMPERATE] = 18.f;
	normal_cond_max[mapssTROPICAL] = 7.6f;
	LaiUpperBoundsEnergy[ICE][GRASS] = 0.f;
	LaiUpperBoundsEnergy[ICE][TREE] = 0.f;
	LaiUpperBoundsEnergy[ICE][SHRUB] = 0.f;
	LaiUpperBoundsEnergy[TUNDRA][GRASS] = 2.f;
	LaiUpperBoundsEnergy[TUNDRA][TREE] = 0.f;
	LaiUpperBoundsEnergy[TUNDRA][SHRUB] = 1.5f;
	LaiUpperBoundsEnergy[TAIGA_TUNDRA][GRASS] = 3.f;
	LaiUpperBoundsEnergy[TAIGA_TUNDRA][TREE] = 0.f;
	LaiUpperBoundsEnergy[TAIGA_TUNDRA][SHRUB] = 2.5f;
	for (ileaf = NEEDLELEAF; ileaf<=BROADLEAF; ileaf++) min_tree_lai[mapssTROPICAL][ileaf] = 0.65f;

	cond_max_nominal[mapssBOREAL][GRASS] = 3.5f;  
	cond_max_nominal[mapssBOREAL][TREE] = 2.5f;
	cond_max_nominal[mapssBOREAL][SHRUB] = 1.5f;  

	cond_max_nominal[mapssTEMPERATE][GRASS] = 3.5f;
	cond_max_nominal[mapssTEMPERATE][TREE] = 2.5f;  
	cond_max_nominal[mapssTEMPERATE][SHRUB] = 1.5f;

	cond_max_nominal[mapssSUBTROPICAL][GRASS] = 5.5f;
	cond_max_nominal[mapssSUBTROPICAL][TREE] = 2.5f;  
	cond_max_nominal[mapssSUBTROPICAL][SHRUB] = 1.5f;

	cond_max_nominal[mapssTROPICAL][GRASS] = 5.5f;
	cond_max_nominal[mapssTROPICAL][TREE] = 7.6f;  
	cond_max_nominal[mapssTROPICAL][SHRUB] = 1.5f;

	cond_min_nominal[mapssBOREAL][GRASS] = 1.0f;  
	cond_min_nominal[mapssBOREAL][TREE] = 1.5f;
	cond_min_nominal[mapssBOREAL][SHRUB] = 0.8f;  

	cond_min_nominal[mapssTEMPERATE][GRASS] = 1.0f;
	cond_min_nominal[mapssTEMPERATE][TREE] = 1.5f;  
	cond_min_nominal[mapssTEMPERATE][SHRUB] = 0.8f;

	cond_min_nominal[mapssSUBTROPICAL][GRASS] = 0.2f;
	cond_min_nominal[mapssSUBTROPICAL][TREE] = 1.5f;  
	cond_min_nominal[mapssSUBTROPICAL][SHRUB] = 0.8f;

	cond_min_nominal[mapssTROPICAL][GRASS] = 0.2f;
	cond_min_nominal[mapssTROPICAL][TREE] = 1.5f;  
	cond_min_nominal[mapssTROPICAL][SHRUB] = 0.8f;

	cond_surface_max[mapssBOREAL][GRASS] = 5.f;  
	cond_surface_max[mapssBOREAL][TREE] = 20.6f;
	cond_surface_max[mapssBOREAL][SHRUB] = 9.4f;  

	cond_surface_max[mapssTEMPERATE][GRASS] = 23.0f;
	cond_surface_max[mapssTEMPERATE][TREE] = 20.6f;  
	cond_surface_max[mapssTEMPERATE][SHRUB] = 9.4f;

	cond_surface_max[mapssSUBTROPICAL][GRASS] = 17.1f;
	cond_surface_max[mapssSUBTROPICAL][TREE] = 12.1f;  
	cond_surface_max[mapssSUBTROPICAL][SHRUB] = 9.4f;

	cond_surface_max[mapssTROPICAL][GRASS] = 17.1f;
	cond_surface_max[mapssTROPICAL][TREE] = 16.1f;  
	cond_surface_max[mapssTROPICAL][SHRUB] = 9.4f;

	cond_min_months[GRASS] = 0; //	allowed months of below min conductance
	cond_min_months[TREE] = 1; 
	cond_min_months[SHRUB] = 1;

	k_transp[mapssBOREAL][GRASS][0][NEEDLELEAF] = 4.00f; //	coefficient of transpiration ratio model
	k_transp[mapssBOREAL][GRASS][0][BROADLEAF] = 4.00f;	

	k_transp[mapssBOREAL][TREE][0][NEEDLELEAF] = 1.50f; 	
	k_transp[mapssBOREAL][TREE][0][BROADLEAF] = 1.95f;	

	k_transp[mapssBOREAL][SHRUB][0][NEEDLELEAF] = 8.75f; 	
	k_transp[mapssBOREAL][SHRUB][0][BROADLEAF] = 8.75f;

	k_transp[mapssTEMPERATE][GRASS][0][NEEDLELEAF] = 4.00f;
	k_transp[mapssTEMPERATE][GRASS][0][BROADLEAF] = 4.00f;

	k_transp[mapssTEMPERATE][TREE][0][NEEDLELEAF] = 2.50f;
	k_transp[mapssTEMPERATE][TREE][0][BROADLEAF] = 2.95f;

	k_transp[mapssTEMPERATE][SHRUB][0][NEEDLELEAF] = 8.75f;
	k_transp[mapssTEMPERATE][SHRUB][0][BROADLEAF] = 8.75f;

	k_transp[mapssSUBTROPICAL][GRASS][0][NEEDLELEAF] = 4.00f;
	k_transp[mapssSUBTROPICAL][GRASS][0][BROADLEAF] = 4.00f;

	k_transp[mapssSUBTROPICAL][TREE][0][NEEDLELEAF] = 2.25f;
	k_transp[mapssSUBTROPICAL][TREE][0][BROADLEAF] = 2.70f;

	k_transp[mapssSUBTROPICAL][SHRUB][0][NEEDLELEAF] = 8.75f;
	k_transp[mapssSUBTROPICAL][SHRUB][0][BROADLEAF] = 8.75f;

	k_transp[mapssTROPICAL][GRASS][0][NEEDLELEAF] = 4.00f;
	k_transp[mapssTROPICAL][GRASS][0][BROADLEAF] = 4.00f;

	k_transp[mapssTROPICAL][TREE][0][NEEDLELEAF] = 2.70f; 
	k_transp[mapssTROPICAL][TREE][0][BROADLEAF] = 2.70f;

	k_transp[mapssTROPICAL][SHRUB][0][NEEDLELEAF] = 6.75f;
	k_transp[mapssTROPICAL][SHRUB][0][BROADLEAF] = 6.75f;

	a_slope[GRASS] = 0.002f;
	a_slope[TREE] = 0.1f;
	a_slope[SHRUB] = 0.03f;

	k_factor_pet_boundary[mapssBOREAL][NEEDLELEAF] = -1325.f;
	k_factor_pet_boundary[mapssBOREAL][BROADLEAF] = -1325.f;

	k_factor_pet_boundary[mapssTEMPERATE][NEEDLELEAF] = -1325.f;
	k_factor_pet_boundary[mapssTEMPERATE][BROADLEAF] = 1325.f;

	k_factor_pet_boundary[mapssSUBTROPICAL][NEEDLELEAF] = -1325.f;
	k_factor_pet_boundary[mapssSUBTROPICAL][BROADLEAF] = 1325.f;

	k_factor_pet_boundary[mapssTROPICAL][NEEDLELEAF] = -1325.f;
	k_factor_pet_boundary[mapssTROPICAL][BROADLEAF] = -1325.f;

	for (cl = 0; cl<NCL; cl++) for (ileaf = NEEDLELEAF; ileaf<=BROADLEAF; ileaf++) 
		k_trans_addend[cl][ileaf] = 0.;

	WhichSoils = 2; // SCS soils from soils_scs.nc

} // end of MAPSS_model_parameters::default_parameters()


void MAPSS_model_parameters::adjust_parameters()
{
	int i;

	for (i = mapssBOREAL; i <= mapssTROPICAL; i++) 
	{
		k_transp[i][GRASS][0][C3GRASS] += k_trans_addend[GRASS][C3GRASS];
		k_transp[i][GRASS][0][C4GRASS] += k_trans_addend[GRASS][C4GRASS];

		k_transp[i][TREE][0][NEEDLELEAF] += k_trans_addend[TREE][NEEDLELEAF];
		k_transp[i][TREE][0][BROADLEAF] += k_trans_addend[TREE][BROADLEAF];

		k_transp[i][SHRUB][0][NEEDLELEAF] += k_trans_addend[SHRUB][NEEDLELEAF];
		k_transp[i][SHRUB][0][BROADLEAF] += k_trans_addend[SHRUB][BROADLEAF];
	}

} // end of adjust_parameters()


void MAPSS_model_parameters::ORNL_model_parameters()
{
	// Changes from default values start here.
	min_tree_lai[mapssBOREAL][NEEDLELEAF] = 1.5f;
	min_tree_lai[mapssBOREAL][BROADLEAF] = 1.5f;

	min_tree_lai[mapssTEMPERATE][NEEDLELEAF] = 1.7f;
	min_tree_lai[mapssTEMPERATE][BROADLEAF] = 1.5f;

	min_tree_lai[mapssSUBTROPICAL][NEEDLELEAF] = 1.9f;
	min_tree_lai[mapssSUBTROPICAL][BROADLEAF] = 1.5f;

	k_transp[mapssBOREAL][GRASS][0][NEEDLELEAF] = 4.65f + 0.10f; //	coefficient of transpiration ratio model
	k_transp[mapssBOREAL][GRASS][0][BROADLEAF] = 4.25f + 0.10f;	

	k_transp[mapssBOREAL][TREE][0][NEEDLELEAF] = 3.3f + -0.74f; 	
	k_transp[mapssBOREAL][TREE][0][BROADLEAF] = 2.0f + -0.74f;	

	k_transp[mapssBOREAL][SHRUB][0][NEEDLELEAF] = 8.35f + -0.80f; 	
	k_transp[mapssBOREAL][SHRUB][0][BROADLEAF] = 9.25f + -0.80f;

	k_transp[mapssTEMPERATE][GRASS][0][NEEDLELEAF] = 3.5f + 0.10f;
	k_transp[mapssTEMPERATE][GRASS][0][BROADLEAF] = 3.5f + 0.10f;

	k_transp[mapssTEMPERATE][TREE][0][NEEDLELEAF] = 3.2f + -0.74f;
	k_transp[mapssTEMPERATE][TREE][0][BROADLEAF] = 3.5f + -0.74f;

	k_transp[mapssTEMPERATE][SHRUB][0][NEEDLELEAF] = 6.6f + -0.80f;
	k_transp[mapssTEMPERATE][SHRUB][0][BROADLEAF] = 9.25f + -0.80f;

	k_transp[mapssSUBTROPICAL][GRASS][0][NEEDLELEAF] = 4.3f + 0.10f;
	k_transp[mapssSUBTROPICAL][GRASS][0][BROADLEAF] += 0.10f;

	k_transp[mapssSUBTROPICAL][TREE][0][NEEDLELEAF] = 3.1f + -0.74f;
	k_transp[mapssSUBTROPICAL][TREE][0][BROADLEAF] = 2.6f + -0.74f;

	k_transp[mapssSUBTROPICAL][SHRUB][0][NEEDLELEAF] = 9.f + -0.80f;
	k_transp[mapssSUBTROPICAL][SHRUB][0][BROADLEAF] = 9.f + -0.80f;

	k_transp[mapssTROPICAL][GRASS][0][NEEDLELEAF] = 4.25f + 0.10f;
	k_transp[mapssTROPICAL][GRASS][0][BROADLEAF] = 3.35f + 0.10f;

	k_transp[mapssTROPICAL][TREE][0][NEEDLELEAF] = 3.5f + -0.74f; 
	k_transp[mapssTROPICAL][TREE][0][BROADLEAF] = 3.8f + -0.74f;

	k_transp[mapssTROPICAL][SHRUB][0][NEEDLELEAF] = 7.25f + -0.80f;
	k_transp[mapssTROPICAL][SHRUB][0][BROADLEAF] = 7.25f + -0.80f;

	taiga_tundra_boundary[mapssTEMPERATE] = 1345.f;

	k_factor_slope[mapssBOREAL][NEEDLELEAF] = -3.5f;
	k_factor_slope[mapssBOREAL][BROADLEAF] = -3.5f;
	k_factor_slope[mapssTEMPERATE][NEEDLELEAF] = -3.5f;
	k_factor_slope[mapssTEMPERATE][BROADLEAF] = -2.9f;
	k_factor_slope[mapssSUBTROPICAL][NEEDLELEAF] = -3.5f;
	k_factor_slope[mapssSUBTROPICAL][BROADLEAF] = -3.5f;
	k_factor_slope[mapssTROPICAL][NEEDLELEAF] = -3.5f;
	k_factor_slope[mapssSUBTROPICAL][BROADLEAF] = -3.5f;

	// The version of MAPSS in the ORNL repository does not distinguish subalpine forest (114) from 
	// other classes.  Setting subalpine_threshold to 0 here effectively disables the logic to 
	// distinguish subalpine forest.
	mapss_subalpine_threshold = 0.f; 

} // end of ORNL_model_parameters()


void MAPSS_model_parameters::US10km_model_parameters(MAPSSparameterSetName paramSet)
{
	// Changes from default values start here.
	min_tree_lai[mapssBOREAL][NEEDLELEAF] = 1.5f;
	min_tree_lai[mapssBOREAL][BROADLEAF] = 1.5f;

	min_tree_lai[mapssTEMPERATE][NEEDLELEAF] = 1.5f;
	min_tree_lai[mapssTEMPERATE][BROADLEAF] = 1.5f;

	min_tree_lai[mapssSUBTROPICAL][NEEDLELEAF] = 1.5f;
	min_tree_lai[mapssSUBTROPICAL][BROADLEAF] = 1.5f;

	k_transp[mapssBOREAL][GRASS][0][NEEDLELEAF] = 4.25f; 
	k_transp[mapssBOREAL][GRASS][0][BROADLEAF] = 4.25f; 

	k_transp[mapssBOREAL][TREE][0][NEEDLELEAF] = 2.75f; 	
	k_transp[mapssBOREAL][TREE][0][BROADLEAF] = 2.75f;	

	k_transp[mapssBOREAL][SHRUB][0][NEEDLELEAF] = 9.25f; 	
	k_transp[mapssBOREAL][SHRUB][0][BROADLEAF] = 9.25f; 	

	k_transp[mapssTEMPERATE][GRASS][0][NEEDLELEAF] = 4.25f;
	k_transp[mapssTEMPERATE][GRASS][0][BROADLEAF] = 4.25f;

	k_transp[mapssTEMPERATE][TREE][0][NEEDLELEAF] = 3.75f;
	k_transp[mapssTEMPERATE][TREE][0][BROADLEAF] = 3.75f;

	k_transp[mapssTEMPERATE][SHRUB][0][NEEDLELEAF] = 9.25f;
	k_transp[mapssTEMPERATE][SHRUB][0][BROADLEAF] = 9.25f;

	k_transp[mapssSUBTROPICAL][GRASS][0][NEEDLELEAF] = 4.25f;
	k_transp[mapssSUBTROPICAL][GRASS][0][BROADLEAF] = 4.25f;

	k_transp[mapssSUBTROPICAL][TREE][0][NEEDLELEAF] = 3.5f;
	k_transp[mapssSUBTROPICAL][TREE][0][BROADLEAF] = 3.5f;

	k_transp[mapssSUBTROPICAL][SHRUB][0][NEEDLELEAF] = 9.25f;
	k_transp[mapssSUBTROPICAL][SHRUB][0][BROADLEAF] = 9.25f;

	k_transp[mapssTROPICAL][GRASS][0][NEEDLELEAF] = 4.25f;
	k_transp[mapssTROPICAL][GRASS][0][BROADLEAF] = 4.25f;

	k_transp[mapssTROPICAL][TREE][0][NEEDLELEAF] = 3.5f;
	k_transp[mapssTROPICAL][TREE][0][BROADLEAF] = 3.5f;

	k_transp[mapssTROPICAL][SHRUB][0][NEEDLELEAF] = 7.25f;
	k_transp[mapssTROPICAL][SHRUB][0][BROADLEAF] = 7.25f;

	tundra_boundary[mapssBOREAL] = 367.5f;
	tundra_boundary[mapssTEMPERATE] = 307.5f;
	tundra_boundary[mapssSUBTROPICAL] = 307.5f;
	tundra_boundary[mapssTROPICAL] = 307.5f;

	taiga_tundra_boundary[mapssBOREAL] = 665.f;
	taiga_tundra_boundary[mapssTEMPERATE] = 582.5f;
	taiga_tundra_boundary[mapssSUBTROPICAL] = 582.5f;
	taiga_tundra_boundary[mapssTROPICAL] = 582.5f;

	k_factor_slope[mapssBOREAL][NEEDLELEAF] = -3.5;
	k_factor_slope[mapssBOREAL][BROADLEAF] = -3.5;
	k_factor_slope[mapssTEMPERATE][NEEDLELEAF] = -3.5;
	k_factor_slope[mapssTEMPERATE][BROADLEAF] = -3.5;
	k_factor_slope[mapssSUBTROPICAL][NEEDLELEAF] = -3.5;
	k_factor_slope[mapssSUBTROPICAL][BROADLEAF] = -3.5;
	k_factor_slope[mapssTROPICAL][NEEDLELEAF] = -3.5;
	k_factor_slope[mapssTROPICAL][BROADLEAF] = -3.5;

	if (paramSet==US10kmFAOparams) WhichSoils = 1; // FAO soils from soils_fao.nc
	else if (paramSet==US10kmSCSparams) WhichSoils = 2; // SCS soils from soils_scs.nc
	else ASSERT(0);

} // end of US10km_model_parameters()


void MAPSS_model_parameters::Wies_Yose_model_parameters()
{
	mapss_subalpine_threshold = 2350; // degree days referenced to 0 C
	s_decid_bound = 0.9f; // deg C
	forest_threshold = 1.5f; // LAI
	for (int zone = 0; zone<NZONES; zone++) min_tree_lai[zone][NEEDLELEAF] = min_tree_lai[zone][BROADLEAF] = 1.0;

} // end of Wies_Yose_model_parameters()


MAPSS_model_parameters::MAPSS_model_parameters(MAPSSparameterSetName paramSet, float in_wue)
{
	default_parameters();

	wue = in_wue;

	// Set transpiration constants.
	for (int zone_index = mapssBOREAL; zone_index <= mapssTROPICAL; zone_index++) {
		for (int cl = GRASS; cl <= SHRUB; cl++) 
		{ // note when cl==GRASS, NEEDLELEAF means C3GRASS, and BROADLEAF means C4GRASS
			// cond_lai_max[zone_index][cl][NEEDLELEAF] = normal_cond_max[zone_index]*LaiUpperBoundsLifeform;
			cond_max[zone_index][cl][NEEDLELEAF] = cond_max_nominal[zone_index][cl]*wue;
			cond_min[zone_index][cl][NEEDLELEAF] = cond_min_nominal[zone_index][cl]*wue;
			c[zone_index][cl][NEEDLELEAF] = 4.0f*cond_max[zone_index][cl][NEEDLELEAF]*cond_max[zone_index][cl][NEEDLELEAF];
			b[zone_index][cl][NEEDLELEAF] = c[zone_index][cl][NEEDLELEAF] / (wp[cl] * wp[cl]);
			k_transp[zone_index][cl][1][NEEDLELEAF] = k_transp[zone_index][cl][0][NEEDLELEAF] * 0.5f;

			// cond_lai_max[zone_index][cl][BROADLEAF] = normal_cond_max[zone_index]*LaiUpperBoundsLifeform;
			cond_max[zone_index][cl][BROADLEAF] = cond_max_nominal[zone_index][cl]*wue;
			cond_min[zone_index][cl][BROADLEAF] = cond_min_nominal[zone_index][cl]*wue;
			c[zone_index][cl][BROADLEAF] = 4.0f*cond_max[zone_index][cl][BROADLEAF]*cond_max[zone_index][cl][BROADLEAF];
			b[zone_index][cl][BROADLEAF] = c[zone_index][cl][BROADLEAF] / (wp[cl] * wp[cl]);
			k_transp[zone_index][cl][1][BROADLEAF] = k_transp[zone_index][cl][0][BROADLEAF] * 0.5f;
		}
	}

	if (paramSet==ORNLparams) 
		ORNL_model_parameters();
	else if (paramSet==US10kmFAOparams || paramSet==US10kmSCSparams)
		US10km_model_parameters(paramSet);

} // end of MAPSS_model_parameters(MAPSSparameterSetName paramSet, float wue) constructor


#ifdef UNIX_MAPSS
MAPSS_BiogeogModel::MAPSS_BiogeogModel(Simulation * pRun)
{
	ProcessModelInit(pRun);
} // end of constructor for class CENTURY_BiogeochemModel
#endif  


MAPSS_BiogeogModel::MAPSS_BiogeogModel(MAPSSparameterSetName paramSet, SoilDataStruct in_soil_data, float * pptP, float * tmpP, float * vprP, float * wndP,
		float in_elev, float in_wue)
{
	int month, zone_index, cl;

	// Set transpiration constants.
	wue = in_wue;
	for (zone_index = mapssBOREAL; zone_index <= mapssTROPICAL; zone_index++) {
		for (cl = GRASS; cl <= SHRUB; cl++) 
		{ // note when cl==GRASS, NEEDLELEAF means C3GRASS, and BROADLEAF means C4GRASS
			// cond_lai_max[zone_index][cl][NEEDLELEAF] = normal_cond_max[zone_index]*LaiUpperBoundsLifeform;
			cond_max[zone_index][cl][NEEDLELEAF] = cond_max_nominal[zone_index][cl]*wue;
			cond_min[zone_index][cl][NEEDLELEAF] = cond_min_nominal[zone_index][cl]*wue;
			c[zone_index][cl][NEEDLELEAF] = 4.0f*cond_max[zone_index][cl][NEEDLELEAF]*cond_max[zone_index][cl][NEEDLELEAF];
			b[zone_index][cl][NEEDLELEAF] = c[zone_index][cl][NEEDLELEAF] / (wp[cl] * wp[cl]);
			k_transp[zone_index][cl][1][NEEDLELEAF] = k_transp[zone_index][cl][0][NEEDLELEAF] * 0.5f;

			// cond_lai_max[zone_index][cl][BROADLEAF] = normal_cond_max[zone_index]*LaiUpperBoundsLifeform;
			cond_max[zone_index][cl][BROADLEAF] = cond_max_nominal[zone_index][cl]*wue;
			cond_min[zone_index][cl][BROADLEAF] = cond_min_nominal[zone_index][cl]*wue;
			c[zone_index][cl][BROADLEAF] = 4.0f*cond_max[zone_index][cl][BROADLEAF]*cond_max[zone_index][cl][BROADLEAF];
			b[zone_index][cl][BROADLEAF] = c[zone_index][cl][BROADLEAF] / (wp[cl] * wp[cl]);
			k_transp[zone_index][cl][1][BROADLEAF] = k_transp[zone_index][cl][0][BROADLEAF] * 0.5f;
		}
	}

	if (paramSet==ORNLparams) 
		ORNL_model_parameters();
	else if (paramSet==US10kmFAOparams || paramSet==US10kmSCSparams)
		US10km_model_parameters(paramSet);

	/* Set initial values of internal variables. */
	m_cats.broadleaf = m_cats.deciduous = m_cats.grassfire = m_cats.mixed = FALSE;
	m_cats.mclass = Unknown;
	m_max_grass = 0.0;
	m_equilibrium_factor = 0.0;
	sciFn.FindSlopeYIntercept(evergreen_gdd_north, 100.0, evergreen_gdd_south, 0.0, &m_mix_slope, &m_mix_y_intercept);


	/* Set initial values specific to the current point simulation. */
	for (month = JAN; month<=DEC; month++)
	{
		m_wnd[month] = wndP!=NULL ? *(wndP + month) : 3.5f;
		m_ppt[month] = *(pptP + month);
		m_tmp[month] = *(tmpP + month);
		m_vpr[month] = *(vprP + month);
	}
	m_elevation = in_elev;

	m_soil = Soil_MAPSS(in_soil_data);

	/* set controlling boundaries for lai iteration */
	m_excess_h2o[SURFACE] = m_soil.awc[SURFACE] * EXCESS_H2O;
	m_excess_h2o[INTERMEDIATE] = (m_soil.awc[SURFACE] + m_soil.awc[INTERMEDIATE]) * EXCESS_H2O;
	m_inadeq_h2o[SURFACE] = m_soil.awc[SURFACE] * INADEQ_H2O;
	m_inadeq_h2o[INTERMEDIATE] = (m_soil.awc[SURFACE] + m_soil.awc[INTERMEDIATE]) * INADEQ_H2O;
	m_unavail = m_soil.wilt_pt[SURFACE] + m_soil.wilt_pt[INTERMEDIATE];	
	m_total_saturation[GRASS] = m_soil.swhc[SURFACE];
	m_total_saturation[WOODY] = m_soil.swhc[SURFACE] + m_soil.swhc[INTERMEDIATE];

} // end of MAPSS_BiogeogModel(ParameterSetName paramSet, SoilDataStruct in_soil_data, float * pptP, float * tmpP, float * vprP, float * wndP, float in_elev, float in_wue)


#ifdef UNIX_MAPSS
void MAPSS_BiogeogModel::MAPSSmodelInit(MAPSSparameterSetName paramSet, SoilDataStruct in_soil_data, float * pptP, float * tmpP, float * vprP, float * wndP,
		float in_elev)
{
	int month, zone_index, cl;

	MAPSS_model_parameters();

	// Set transpiration constants.
	for (zone_index = mapssBOREAL; zone_index <= mapssTROPICAL; zone_index++) {
		for (cl = GRASS; cl <= SHRUB; cl++) 
		{ // note when cl==GRASS, NEEDLELEAF means C3GRASS, and BROADLEAF means C4GRASS
			cond_lai_max[zone_index][cl][NEEDLELEAF] = normal_cond_max[zone_index]*LaiUpperBoundsLifeform;
			cond_max[zone_index][cl][NEEDLELEAF] = cond_max_nominal[zone_index][cl];
			cond_min[zone_index][cl][NEEDLELEAF] = cond_min_nominal[zone_index][cl];
			c[zone_index][cl][NEEDLELEAF] = 4.0f*cond_max[zone_index][cl][NEEDLELEAF]*cond_max[zone_index][cl][NEEDLELEAF];
			b[zone_index][cl][NEEDLELEAF] = c[zone_index][cl][NEEDLELEAF] / (wp[cl] * wp[cl]);
			k_transp[zone_index][cl][1][NEEDLELEAF] = k_transp[zone_index][cl][0][NEEDLELEAF] * 0.5f;

			cond_lai_max[zone_index][cl][BROADLEAF] = normal_cond_max[zone_index]*LaiUpperBoundsLifeform;
			cond_max[zone_index][cl][BROADLEAF] = cond_max_nominal[zone_index][cl];
			cond_min[zone_index][cl][BROADLEAF] = cond_min_nominal[zone_index][cl];
			c[zone_index][cl][BROADLEAF] = 4.0f*cond_max[zone_index][cl][BROADLEAF]*cond_max[zone_index][cl][BROADLEAF];
			b[zone_index][cl][BROADLEAF] = c[zone_index][cl][BROADLEAF] / (wp[cl] * wp[cl]);
			k_transp[zone_index][cl][1][BROADLEAF] = k_transp[zone_index][cl][0][BROADLEAF] * 0.5f;
		}
	}

	if (paramSet==ORNLparams) 
		ORNL_model_parameters();
	else if (paramSet==US10kmFAOparams || paramSet==US10kmSCSparams)
		US10km_model_parameters(paramSet);
	else if (paramSet==WiesYoseParams) Wies_Yose_model_parameters();

	if (modelParamsP->MAPSS_subalpine_threshold!=-9999.) mapss_subalpine_threshold = modelParamsP->MAPSS_subalpine_threshold;

	adjust_parameters();

	/* Set initial values of internal variables. */
	m_cats.broadleaf = m_cats.deciduous = m_cats.grassfire = m_cats.mixed = FALSE;
	m_cats.mclass = Unknown;
	m_max_grass = 0.0;
	m_equilibrium_factor = 0.0;
	sciFn.FindSlopeYIntercept(evergreen_gdd_north, 100.0, evergreen_gdd_south, 0.0, &m_mix_slope, &m_mix_y_intercept);


	/* Set initial values specific to the current point simulation. */
	for (month = JAN; month<=DEC; month++)
	{
		m_wnd[month] = wndP!=NULL ? *(wndP + month) : 3.5f;
		m_ppt[month] = *(pptP + month);
		m_tmp[month] = *(tmpP + month);
		m_vpr[month] = *(vprP + month);
	}
	m_elevation = in_elev;

	m_soil = Soil_MAPSS(in_soil_data);

	/* set controlling boundaries for lai iteration */
	m_excess_h2o[SURFACE] = m_soil.awc[SURFACE] * EXCESS_H2O;
	m_excess_h2o[INTERMEDIATE] = (m_soil.awc[SURFACE] + m_soil.awc[INTERMEDIATE]) * EXCESS_H2O;
	m_inadeq_h2o[SURFACE] = m_soil.awc[SURFACE] * INADEQ_H2O;
	m_inadeq_h2o[INTERMEDIATE] = (m_soil.awc[SURFACE] + m_soil.awc[INTERMEDIATE]) * INADEQ_H2O;
	m_unavail = m_soil.wilt_pt[SURFACE] + m_soil.wilt_pt[INTERMEDIATE];	
	m_total_saturation[GRASS] = m_soil.swhc[SURFACE];
	m_total_saturation[WOODY] = m_soil.swhc[SURFACE] + m_soil.swhc[INTERMEDIATE];

} // end of MAPSS_BiogeogModel::MAPSSmodelInit(MAPSSparameterSetName paramSet, SoilDataStruct in_soil_data, float * pptP, float * tmpP, float * vprP, float * wndP, float in_elev)


bool MAPSS_BiogeogModel::runModelAndCollectData(const int year_index, const int row_index, const int col_index)
{
	bool rtn_flag;  
	size_t coords[3];
	int rtn_val;
	MAPSSvegClass mclass;

	bool ckLatLonFlag = row_index<=1 && col_index<=1;
	int grid_row = row_index + runParamsP->row_begin;
	int grid_col = col_index + runParamsP->col_begin;

	/* Get the input data for this point. */
	InputDataClass input_data;
	inputDataP = &input_data; // inputDataP is member data of class ProcessModel
	rtn_flag = pS->getInputData(0, grid_row, grid_col, ckLatLonFlag, &input_data);
	if (!rtn_flag) 
	{
		printf("getInputData() failed\n");
		return(false);
	}

	if (runParamsP->latlonFlag) modelParamsP->southernHemisphereFlag |= input_data.lat<0.;

	MAPSSmodelInit(modelParamsP->MAPSSparameterSet, input_data.soilData,
			input_data.pptP, input_data.tmpP, input_data.vprP, input_data.wndP, input_data.elev);

	/* Run the MAPSS model. */
	MAPSSdataOutputs mapssOutput;
	mclass = mapss_model(&mapssOutput); // Run the MAPSS model on one point.
	mapssOutput.mapssOutvars[MAPSSmclass] = (float)mclass;
	mapssOutput.mapssOutvars[MAPSSvclass] = (float)convertMAPSStoVEMAP(mclass);
	mapssOutput.mapssOutvars[MAPSScanopy] = (float)m_canopy_type;
	mapssOutput.mapssOutvars[MAPSSzone] = (float)m_zone;
	// mapssOutput.mapssOutvars[MAPSS ] =  ;

	rtn_flag = mclass>Unknown && mclass<=Ice; 
	if (!rtn_flag) 
	{
		printf("MAPSS model failed to return a vegetation type: mclass, msg = %d, %s\n",
				mclass, mapssOutput.msgP==NULL ? "no message" : mapssOutput.msgP);
		return(false);
	}
	/* Collect the output data for this point. */  
	coords[pS->MAPSSoutFile_ncid.rowNdxPos] = grid_row - pS->MAPSSoutFile_ncid.row_offset;
	coords[pS->MAPSSoutFile_ncid.colNdxPos] = grid_col - pS->MAPSSoutFile_ncid.col_offset;
	assert(pS->MAPSSoutFile_ncid.timeNdxPos>=0);
	coords[pS->MAPSSoutFile_ncid.timeNdxPos] = 0;

	for (int i = 0; i<NUM_MAPSS_OUTVARS; i++) if (pS->MAPSSoutputList[i].show) switch (pS->MAPSSoutputList[i].type)
	{ short shortval; float floatval;
		case NC_SHORT:
			shortval = (short)mapssOutput.mapssOutvars[i];
			rtn_val = nc_put_var1_short(pS->MAPSSoutFile_ncid.fileid, pS->MAPSSoutputList[i].var_id, coords, &shortval); 
			if (!(chk_nc(rtn_val)))
				assert(false);
			break;
		case NC_FLOAT:
			floatval = (float)mapssOutput.mapssOutvars[i];
			rtn_val = nc_put_var1_float(pS->MAPSSoutFile_ncid.fileid, pS->MAPSSoutputList[i].var_id, coords, &floatval); assert(chk_nc(rtn_val));
			break;
		default: assert(0); break;
	}

	return(true);

} // end of MAPSS_BiogeogModel::runModelAndCollectData()
#endif

Soil_MAPSS::Soil_MAPSS() { }
Soil_MAPSS::~Soil_MAPSS() { }

Soil_MAPSS::Soil_MAPSS(SoilDataStruct soil_data)
{
	{ // This code was originally in the default constructor for Soil_MAPSS
		/************
		  from Table 2, Saxton et al., pg 1033,
		  Soil Society of America Journal, Vol. 50, 1986
		 *************/
		a = -4.396;
		b = -0.0715;
		c = -4.880e-4;
		d = -4.285e-5;
		e = -3.140;
		f = -2.22e-3;
		g = -3.484e-5;
		h = 0.332;
		j = -7.251e-4;
		k = 0.1276;
		m = -0.108;
		n = 0.341;
		rock_frag_max = 0.50;
		k_surfrun = 1.7f; //	coeff of surface runoff (increase to reduce runoff) 
	} 

	float thick_top2layers;
	double thickA, thickB; // These have to be doubles to get the exact same answers as the old MC1 code.          
	int sl;
	const float thicknessI[NSL] = {500., 1000., 1800.}; // mm
	const double ppI[NSL] = {21.012, 21.012, 19.012};
	const double qqI[NSL] = {-8.0e-2, -8.0e-2, -8.0e-2};
	const double rrI[NSL] = {-3.895, -3.895, -4.3};
	const double ttI[NSL] = {4.55e-2, 3.92e-2, 3.42e-2};
	const double uuI[NSL] = {-0.03, -0.04, -0.05};
	const double vvI[NSL] = {8.7546e-5, 8.7546e-5, 8.7546e-5};
	const static float KK_I[2][NSL] = {
		{0.5f, 0.8f, 0.8f}, // saturated
		{0.8f, 0.5f, 0.2f}}; // unsaturated
	const static float exp_percI[2][NSL] = {
		{1.0f, 2.5f, 10.0f}, // saturated
		{2.5f, 3.0f, 10.0f}}; // unsaturated 
	const float matrix_potI = -1.5f; // MPa
	const float field_potI = -0.033f;

	for (sl=0; sl<NSL; sl++)
	{
		thickness[sl] = thicknessI[sl];
		pp[sl] = ppI[sl];
		qq[sl] = qqI[sl];
		rr[sl] = rrI[sl];
		tt[sl] = ttI[sl];
		uu[sl] = uuI[sl];
		vv[sl] = vvI[sl];
		KK[SATURATED][sl] = KK_I[SATURATED][sl];
		KK[UNSATURATED][sl] = KK_I[UNSATURATED][sl];
		exp_perc[SATURATED][sl] = exp_percI[SATURATED][sl];
		exp_perc[UNSATURATED][sl] = exp_percI[UNSATURATED][sl];
	}
	matrix_pot = matrix_potI;
	field_pot = field_potI;

	thick_top2layers = thickness[SURFACE] + thickness[INTERMEDIATE];
	mineral_depth = soil_data.soils_record[0];
	for (sl = 0; sl<NUMBER_SOIL_LAYERS; sl++)
	{
		sand[sl] = soil_data.soils_record[1 + sl];
		clay[sl] = soil_data.soils_record[4 + sl];
		rock_frag_mineral[sl] = soil_data.soils_record[7 + sl];
		// printf("MC2 sl, sand, clay, rock_frag_mineral = %d, %f, %f, %f\n", sl, sand[sl], clay[sl], rock_frag_mineral[sl]);
	}
	bulk_density = soil_data.bd;
	// In the lower limit check for bulk density below, we use a limit of 0.01 g/cm3, relative to a lowest bulk
	// density in nature of 0.03 g/cm3, for histosols (e.g. peat soils), as reported in 
	//     Batjes NH (1996).  Total carbon and nitrogen in soils of the world.
	//     European Journal of Soil Science 47:151-163.
	// Thanks to Wendy Peterman for providing this information.

	valid = (bulk_density==NO_DATA_TOKEN || (bulk_density>=0.01 && bulk_density<=3.0)) 
		&& sand[SURFACE]>0.0 && clay[SURFACE]>0.0 && rock_frag_mineral[SURFACE]>=0.0 && mineral_depth>10.0;
	if (!valid) return;

	if (sand[INTERMEDIATE] <= 0.0) sand[INTERMEDIATE] = sand[SURFACE];
	if (clay[INTERMEDIATE] <= 0.0) clay[INTERMEDIATE] = clay[SURFACE];
	if (rock_frag_mineral[INTERMEDIATE] < 0.0) rock_frag_mineral[INTERMEDIATE] = rock_frag_mineral[SURFACE];

	if (sand[DEEP] <= 0.0) sand[DEEP] = sand[INTERMEDIATE];
	if (clay[DEEP] <= 0.0) clay[DEEP] = clay[INTERMEDIATE];
	if (rock_frag_mineral[DEEP] < 0.0) rock_frag_mineral[DEEP] = rock_frag_mineral[INTERMEDIATE];

	// Set effective thickness of soil layers.
	if (mineral_depth >= thick_top2layers) for (sl = SURFACE; sl <= DEEP; sl++) 
		eff_thickness[sl] = thickness[sl] * (1.0 - (rock_frag_mineral[sl] / 100.0));
	else if ((mineral_depth > thickness[SURFACE]) && (mineral_depth <= thick_top2layers)) 
	{
		eff_thickness[SURFACE] = thickness[SURFACE] * (1.0 - (rock_frag_mineral[SURFACE] / 100.0));
		thickA = (mineral_depth - thickness[SURFACE]) * (1.0 - (rock_frag_mineral[INTERMEDIATE] / 100.0));
		thickB = thickness[INTERMEDIATE] * (1.0 - rock_frag_max);
		eff_thickness[INTERMEDIATE] = MAX(thickA, thickB);
		eff_thickness[DEEP] = thickness[DEEP] * (1.0 - rock_frag_max);
	}
	else if (mineral_depth <= thickness[SURFACE]) 
	{
		thickA = mineral_depth * (1.0 - (rock_frag_mineral[SURFACE] / 100.0));
		thickB = thickness[SURFACE] * (1.0 - rock_frag_max);
		eff_thickness[SURFACE] = MAX(thickA, thickB);
		eff_thickness[INTERMEDIATE] = thickness[INTERMEDIATE] * (1.0 - rock_frag_max);
		eff_thickness[DEEP] = thickness[DEEP] * (1.0 - rock_frag_max);
	}
	else ASSERT(0);

	/* set global drainage variables */
	// Originally these calculations were all done in KPa.  MAPSS works in
	// MPa.  The numbers taken as parameters are in MPa and the numbers
	// given as results are in mm.
	// wilt_pt and theta_m return the same values.
	for (sl = SURFACE; sl <= DEEP; sl++) 
	{
		double bb_inv;
		double spct = sand[sl];
		double cpct = clay[sl];
		double thick = eff_thickness[sl];

		aa[sl] = exp(a + b*cpct + c*spct*spct + d*spct*spct*cpct)*100.0;
		bb[sl] = e + f*cpct*cpct + g*spct*spct + g*spct*spct*cpct;

		bb_inv = 1./bb[sl];
		field_cap[sl] = (float)(pow((-field_pot*1000.)/aa[sl], bb_inv)*thick);
		wilt_pt[sl] = (float)(pow((-matrix_pot*1000.)/aa[sl], bb_inv)*thick);
		theta_m[sl] = (float)(pow((-matrix_pot*1000.)/aa[sl], bb_inv)*thick);

		h2o_sat[sl] = h + (j * spct) + (k * log10(cpct));
		swhc[sl] = (float)(h2o_sat[sl]*thick);
		awc[sl] = field_cap[sl] - wilt_pt[sl];
		if (diagsMAPSS)
		{
#ifndef UNIX_MAPSS
			diagMsg.Format("MC2 sl, aa, bb_inv, thick = %d, %.8lf, %.8lf, %.8lf\n", sl, aa[sl], bb_inv, thick);
			Report::LogMsg(diagMsg);
			diagMsg.Format("MC2 sl, field_cap, wilt_pt, awc = %d, %.8f, %.8f, %.8f\n", sl, field_cap[sl], wilt_pt[sl], awc[sl]); 
			Report::LogMsg(diagMsg);
#endif
		} 
		sat2mat[sl] = swhc[sl] - wilt_pt[sl];
	}

} // end of Soil_MAPSS(SoilDataStruct soil_data)


/*
 * 
 * ----------------------------  MAPSS biogeography model
 * 
 * Returns MAPSS vegetation type directly, or an error code<=0. 
 * If returning a veg type, loads up the mapss_output structure 
 * as a side effect.  Can abort due to an error.
 *
 * This model is described in Neilson R.P. (1995). A model for predicting
 * continental-scale vegetation distribution and water balance. Ecological
 * Applications 5(2): 362-385.
 *
 * The original MAPSS model C code was incorporated into the MC1 model.
 * This C++ version of the MAPSS model was adapted from the MAPSS C code
 * in the MC1 model by David Conklin at Conservation Biology Institute in 2010-11.
 */
MAPSSvegClass MAPSS_BiogeogModel::mapss_model(MAPSSdataOutputs * in_mapss_outputP)
{
	int mo, veg, cl;
	State site[MONTHS*2];
	float max_ppt;
	State	* curr_year;
	State * prev_year;

	m_mapss_outputP = in_mapss_outputP;
	memset(m_mapss_outputP, 0, sizeof(MAPSSdataOutputs));
	m_mapss_outputP->msgP = "";

	// Prepare the climate data
	m_min_tmp = m_tmp[JAN];
	m_max_tmp = m_tmp[JAN];
	m_max_tmp_mo = JAN;
	for (mo = FEB; mo <= DEC; mo++) 
	{
		if (m_tmp[mo] < m_min_tmp) m_min_tmp = m_tmp[mo];
		if (m_tmp[mo] > m_max_tmp) 
		{
			m_max_tmp = m_tmp[mo];
			m_max_tmp_mo = mo;
		}
	}

	{ // Calculate m_thaw and m_length_winter
		// thaw will be set to 1=Jan,...,12=Dec or FROST_ANY_MONTH or NEVER_FREEZES.
		// All we have for input is mean monthly temperature, not minimum temperature.  So there 
		// can be freezing temperatures in a given month even if the tmp[mo] value is above 0.
		// The "frost" parameter is the mean monthly temperature below which we assume there are
		// freezing temperatures. As of the date of these comments, "frost" has a value of 13.0 C,
		// the "s_decid_bound" parameter is +2.25 C, and "n_decid_bound" is -16.0 C. (12/30/10)
		int	mo, thaw_mo = -1;
		bool	does_freeze = FALSE;

		for (mo = JAN, m_length_winter = 0; mo <= DEC; mo++) 
		{
			if (m_tmp[mo] < frost) 
			{
				does_freeze = TRUE;
				m_length_winter++;
			}
			else if (thaw_mo == -1) thaw_mo = mo;
		}
		if (thaw_mo == -1) m_thaw = FROST_ANY_MONTH;
		else if (!does_freeze) m_thaw = NEVER_FREEZES;
		else m_thaw = thaw_mo + 1; // m_thaw = 1 for Jan, ..., 12 for Dec
	} // end of block to calculate m_thaw and m_length_winter

	// Determine the climate zone.
	m_zone = sciFn.ClimateZone4MAPSS(m_min_tmp, n_decid_bound, s_decid_bound, frost);

	if (!m_soil.valid) // Soil data was prepared by the constructor.
	{
		m_mapss_outputP->msgP = INVALID_SOIL_DATA_MSG;
		return(Unknown);
	}

	/* Calc monthly vp_sat[] and monthly potential evapotranspiration for shrub and tree layers */
	for (mo = JAN; mo <= DEC; mo++) 
	{
		m_vp_sat[mo] = (float)sciFn.satvp((double)m_tmp[mo]);
		if (m_vp_sat[mo] < m_vpr[mo] ) m_vp_sat[mo] = m_vpr[mo] + 1.0f; // Constrain m_vp_sat to >= m_vpr + 1.
		for (cl = 0; cl<NCL; cl++) m_pet_array[cl][mo] = 0.0f;
	}
	CalcTrbxfr(TREE);
	CalcTrbxfr(SHRUB);
	CalcTrbxfr(GRASS);
	m_petP = m_pet_array[TREE]; // m_petP points to the 12 monthly PET values for TREEs

	/* calculate storm event frequency, melt rate, etc */
	{ // Snow and melt - For each month, calculates amount of precip which falls as snow and the melt rate. 
		// Snow is in units of snow water equivalent.
		float	b, a, a_melt;

		b = -1.0f / (snow0 - snow1);			/* slope */
		a = 1.0f + b * (- snow1);			/* slope */
		a_melt = melt_slope * (- no_melt);	/* melt y-intercept */

		for (mo = JAN; mo <= DEC; mo++) 
		{
			if (m_tmp[mo] <= snow1) m_months[mo].snow = m_ppt[mo];
			else if (m_tmp[mo] >= snow0) m_months[mo].snow = 0.0;
			else m_months[mo].snow = (a + (b * m_tmp[mo])) * m_ppt[mo];

			if (m_tmp[mo] <= no_melt) m_months[mo].melt_rate = 0.0;
			else m_months[mo].melt_rate = a_melt + (melt_slope * m_tmp[mo]);

		}	
	} // end of snow and melt block  
	{ // Rainfall events
		int	season, length_summer;
		float	rain;
		float events_summer;

		m_ppt_winter = 0.0;
		events_summer = 0.0;

		for (mo = JAN; mo < MONTHS; mo++) 
		{
			rain = m_ppt[mo] - m_months[mo].snow;
			season = (*(m_petP + mo) > event_pet) ? 1 : 0; // if PET for trees is > event_pet (e.g. 50), then season==1
			m_months[mo].events = MIN(event_ppt * rain, max_events[season]);

			if (m_tmp[mo] >= frost) events_summer += m_months[mo].events;
			if (m_tmp[mo] < frost) m_ppt_winter += m_ppt[mo];
		}

		length_summer = 12 - m_length_winter;
		m_events_summer_avg = length_summer ? events_summer/length_summer : 0.0f;
	} // end of Rainfall events block

	m_growing_deg_days_zero = sciFn.GrowingDegreeDays(m_tmp, BASE_ZERO);
	m_growing_deg_days_frost = sciFn.GrowingDegreeDays(m_tmp, frost);
	m_mix_ratio = sciFn.FindRatio(m_growing_deg_days_frost, m_mix_slope, m_mix_y_intercept);
	BroadDecid(&m_cats, -1.0, FALSE, FALSE); // Distinguish between EB, DB, and EN

	max_ppt = -1.0;
	for (mo = JAN; mo <= DEC; mo ++) if (m_ppt[mo] > max_ppt) max_ppt = m_ppt[mo];
	if (max_ppt <= PPT_CUTOFF) 
	{		
		if (m_growing_deg_days_zero <= ice_boundary) return(Ice);
		else return( DesertExtreme );
	}

	// Set m_LaiUpperBoundEnergy for each veg type from the LaiUpperBoundsEnergy parameters for the appropriate m_zone.
	for (veg = GRASS; veg <= SHRUB; veg++)
	{
		if (m_growing_deg_days_zero > taiga_tundra_boundary[m_zone]) m_LaiUpperBoundEnergy[veg] = 
			LaiUpperBoundsEnergy[m_zone][veg];
		else if (m_growing_deg_days_zero <= ice_boundary)  m_LaiUpperBoundEnergy[veg] = LaiUpperBoundsEnergy[ICE][veg];
		else if (m_growing_deg_days_zero <= tundra_boundary[m_zone]) m_LaiUpperBoundEnergy[veg] = 
			LaiUpperBoundsEnergy[TUNDRA][veg];
		else m_LaiUpperBoundEnergy[veg] = LaiUpperBoundsEnergy[TAIGA_TUNDRA][veg];
	}

	// Initialize monthly lai and pet values.
	m_maxlai = initialize_lai();

	// Initialize site, the array of 24 months of State structures.
	for (mo = 0; mo<MONTHS*2; mo++)
	{
		site[mo].transpire[GRASS] = site[mo].transpire[WOODY] = 0.;
		site[mo].evap = site[mo].infiltrate = 0.;
		site[mo].soil_h2o[SURFACE] = site[mo].soil_h2o[INTERMEDIATE] = site[mo].soil_h2o[DEEP] = 0.;
		site[mo].snow = site[mo].melt = site[mo].surf_run = site[mo].base_run = 0.;
		site[mo].swp[SURFACE] = site[mo].swp[INTERMEDIATE] = -1000.;
		site[mo].pct_soil_h2o[SURFACE] = site[mo].pct_soil_h2o[INTERMEDIATE] = 0.;
	}

	LaiCycle(&m_offset, site, m_offset, &curr_year, &prev_year); /* Iterate through lai cycle */

	GrassWoodyPart1(curr_year, prev_year, site, m_offset);
	if (!m_cats.deciduous || !CheckEvergreen(curr_year, site, m_offset, m_lai_values, &curr_year, &prev_year))
	{ // GrassWoodyPart2
		bool did_burn, shrubs;
		PsPath biomass;

		shrubs = !(m_pet_adj == tree_pet_factor);

		{ // Calculate fire_condition
			float shrub_lai;
			unsigned		mo;

			did_burn = FALSE;

			shrub_lai = (m_pet_adj == tree_pet_factor) ? 0.0f : m_maxlai;
			if (((m_lai_values[GRASS_SUM] > fire_threshold_1) || (shrub_lai > fire_threshold_2))) 
			{
				for (mo = JAN; mo <= DEC; mo++) if ((*(m_petP + mo) > fire_threshold_5) && (m_ppt[mo] > fire_threshold_6)) 
				{
					did_burn = TRUE;
					break;
				}
			}
		} // end of fire_condition block

		if ((!(((m_zone == mapssSUBTROPICAL) || (m_zone == mapssTROPICAL)) && (!m_cats.broadleaf))) && (shrubs)) 
			m_cats.grassfire = did_burn;

		if (m_cats.grassfire) 
		{
			GrassAlone(prev_year, m_conductance);
			if (0 && diagsMAPSS) for (mo = 0; mo<12; mo++) printf("MC2 after GrassAlone m_lai[%d][] = %f, %f\n", mo, m_lai[mo][0], m_lai[mo][1]);
			m_lai_values[GRASS] = FindMeanLai(&m_lai_values[GRASS_SUM], GRASS);
			m_lai_values[TREE] = m_lai_values[SHRUB] = m_lai_values[TREE_SUM] = m_lai_values[SHRUB_SUM] = 0.0;
		}
		else 
		{
			if (m_pet_adj == tree_pet_factor) 
			{
				FindMeanLai(&m_lai_values[TREE_SUM], WOODY);
				m_lai_values[TREE] = m_maxlai;
			}
			else 
			{
				FindMeanLai(&m_lai_values[SHRUB_SUM], WOODY);
				m_lai_values[SHRUB] = m_maxlai;
			}
		}

		m_max_grass = 0.0;
		for (mo = JAN; mo <= DEC; mo++) if (m_lai[mo][GRASS] > m_max_grass) m_max_grass = m_lai[mo][GRASS];

		PS_CenturyNew(m_lai, &biomass);
		m_canopy_type = biomass.canopy;
		// m_c3pct = biomass.ratio;

		m_cats.mclass = ClassifyStation(m_lai_values, m_canopy_type);

		for (mo = JAN; mo <= DEC; mo++) *(m_petP + mo) *= m_pet_adj;
		m_capacity = m_soil.awc[SURFACE] + m_soil.awc[INTERMEDIATE];

	} // end of GrassWoodyPart2 

	/*Load output structure.*/
	m_mapss_outputP->mclass = m_cats.mclass;
	for (mo = JAN; mo <= DEC; mo++) 
	{
		m_mapss_outputP->lai[mo][WOODY] = m_lai[mo][WOODY];
		m_mapss_outputP->lai[mo][GRASS] = m_lai[mo][GRASS];
		m_mapss_outputP->conductance[mo][WOODY] = m_conductance[mo][WOODY];
		m_mapss_outputP->conductance[mo][GRASS] = m_conductance[mo][GRASS];
	}
	m_mapss_outputP->lai_values[GRASS] = m_lai_values[GRASS];
	m_mapss_outputP->lai_values[TREE] = m_lai_values[TREE];
	m_mapss_outputP->lai_values[SHRUB] = m_lai_values[SHRUB];
	m_mapss_outputP->lai_values[GRASS_SUM] = m_lai_values[GRASS_SUM];
	m_mapss_outputP->lai_values[TREE_SUM] = m_lai_values[TREE_SUM];
	m_mapss_outputP->lai_values[SHRUB_SUM] = m_lai_values[SHRUB_SUM];  
	m_mapss_outputP->at_woody_ppt = m_at_woody_ppt;
	m_mapss_outputP->at_woody_ppt_norm = m_at_woody_ppt_norm;
	m_mapss_outputP->old_tree_lai = m_old_tree_lai;
	m_mapss_outputP->new_tree_lai = m_new_tree_lai;
	m_mapss_outputP->woody_alone_at = m_woody_alone_at;
	m_mapss_outputP->k_factor = m_k_factor;
	m_mapss_outputP->canopy_type = m_canopy_type;
	m_mapss_outputP->h2o_capacity = m_capacity;    
	m_mapss_outputP->zone = m_zone;  
	// m_mapss_outputP->c3c4_ratio = m_c3pct;

	return(m_cats.mclass);

} // end of mapss_model()


