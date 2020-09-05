/*
 *  CENTURY.cpp
 *  mc2
 *
 *  Created by Dave Conklin on 3/14/12.
 *
 */

#include <stdio.h> // for printf()
#include <math.h> 
#include <cstring> 

#include "assert.h"

#include "netcdf.h"
#include "category_bgc.h"
#include "MAPSSvegClasses.h"

#include "ScienceFcns.h"
#include "ProcessModel.h"
// #include "ScienceFcns.h"
#include "CENTURY.h"
#include "MAPSSbiogeographyModel.h"
#include "MCfire.h"
#include "MC2.h"
#include "MCbiogeog.h"


CENTURY_BiogeochemModel::CENTURY_BiogeochemModel(Simulation * pRun)
{
	ProcessModelInit(pRun);
	m_fire_last_month = false;
} // end of constructor for class CENTURY_BiogeochemModel


bool CENTURY_BiogeochemModel::runModelAndCollectData(const int year_index, const int row_index, const int col_index)
{
	bool rtn_flag;  
	size_t coords[3];
	int rtn_val;

	bool ckLatLonFlag = row_index<=1 || col_index<=1 || pS->m_num_active_cells<10;
	int grid_row = row_index + runParamsP->row_begin;
	int grid_col = col_index + runParamsP->col_begin;
	int month = 12*(year_index + runParamsP->years_offset_into_input_data);

	InputDataClass input_data;
	inputDataP = &input_data; // inputDataP is member data of class ProcessModel
	rtn_flag = pS->getInputData(month, grid_row, grid_col, ckLatLonFlag, &input_data);
	if (!rtn_flag) return(false);

	// Wind data, when available at all, comes as mean wind for each month - no interannual variation.
	float dummy_wnd[12] = {3.5, 3.5, 3.5,  3.5, 3.5, 3.5,  3.5, 3.5, 3.5,  3.5, 3.5, 3.5}; // m s-1 
	for (int mo = 0; mo<12; mo++)
	{
		int source_mo = modelParamsP->southernHemisphereFlag ? (mo + 6)%12 : mo;
		if (inputDataP->wndP!=NULL)  m_wnd[mo] = *(inputDataP->wndP + source_mo);
		else
		{ 
			assert(runParamsP->dummyWindFlag);
			m_wnd[mo] = dummy_wnd[source_mo];
		}
	} // end of month loop

	switch (runParamsP->runMode)
	{
		case CENTURY_EQ:
			mapssOutput.mclass = (int)*(input_data.inFromMAPSS + MAPSSmclass); 
			mapssOutput.zone = (int)*(input_data.inFromMAPSS + MAPSSzone); 

			// Execute the EQ model      
			rtn_flag = century_eq_model(); 
			if (!rtn_flag) return(false);

			// Collect the EQ output data for this point. 
			EQdata.cenOutvars[CENmapss_canopy] = *(input_data.inFromMAPSS + MAPSScanopy);
			EQdata.cenOutvars[CENmapss_zone] = (int)*(input_data.inFromMAPSS + MAPSSzone);

			coords[pS->EQoutFile_ncid.rowNdxPos] = grid_row - pS->EQoutFile_ncid.row_offset;
			coords[pS->EQoutFile_ncid.colNdxPos] = grid_col - pS->EQoutFile_ncid.col_offset;
			if (pS->EQoutFile_ncid.timeNdxPos>=0) coords[pS->EQoutFile_ncid.timeNdxPos] = 0;

			if (modelParamsP->code_flags[FAST_WARMSTART_OUT_FLAG])
			{
				size_t ws_coords[4], ws_counts[4];
				ws_coords[pS->EQoutFile_ncid.timeNdxPos] = 0;
				ws_coords[pS->EQoutFile_ncid.rowNdxPos] = coords[pS->EQoutFile_ncid.rowNdxPos];
				ws_coords[pS->EQoutFile_ncid.colNdxPos] = coords[pS->EQoutFile_ncid.colNdxPos];
				ws_counts[pS->EQoutFile_ncid.timeNdxPos] = ws_counts[pS->EQoutFile_ncid.rowNdxPos] = ws_counts[pS->EQoutFile_ncid.colNdxPos] = 1;
				ws_coords[pS->EQoutFile_ncid.wsNdxPos] = 0;
				ws_counts[pS->EQoutFile_ncid.wsNdxPos] = NUM_EQ_OUTVARS;
				rtn_val = nc_put_vara(pS->EQoutFile_ncid.fileid, pS->EQoutFile_ncid.ws_varid, ws_coords,
						ws_counts, EQdata.cenOutvars);
				assert(chk_nc(rtn_val));
			}
			else for (int i = 0; i<NUM_EQ_OUTVARS; i++) if (pS->EQandWSoutputList[i].show) 
			{ float floatval;
				floatval = (float)EQdata.cenOutvars[i];
				rtn_val = nc_put_var1_float(pS->EQoutFile_ncid.fileid, pS->EQandWSoutputList[i].var_id, coords, &floatval); 
				if (!chk_nc(rtn_val))
					assert(false);
			}

			return(true);
			break;

		case SPINUP:
			for (int i = 0; i<NUM_WS_OUTVARS; i++) EQdata.cenOutvars[i] = *(input_data.inFromEQ + i);

			// Run spinup for 1 point.
			rtn_flag = (EQdata.cenOutvars[CENinitial_class_mapss]!=NC_FILL_FLOAT) &&
				mc2_spinup(); 
			if (!rtn_flag) return(false);

			// Collect the spinup output data for this point. 
			coords[pS->WSoutFile_ncid.rowNdxPos] = grid_row - pS->WSoutFile_ncid.row_offset;
			coords[pS->WSoutFile_ncid.colNdxPos] = grid_col - pS->WSoutFile_ncid.col_offset;
			if (pS->WSoutFile_ncid.timeNdxPos>=0) coords[pS->WSoutFile_ncid.timeNdxPos] = 0;

			stateData.copyFStoCENoutvars();

			if (modelParamsP->code_flags[FAST_WARMSTART_OUT_FLAG])
			{
				size_t ws_coords[4], ws_counts[4];
				ws_coords[pS->WSoutFile_ncid.timeNdxPos] = 0;
				ws_coords[pS->WSoutFile_ncid.rowNdxPos] = coords[pS->WSoutFile_ncid.rowNdxPos];
				ws_coords[pS->WSoutFile_ncid.colNdxPos] = coords[pS->WSoutFile_ncid.colNdxPos];
				ws_counts[pS->WSoutFile_ncid.timeNdxPos] = ws_counts[pS->WSoutFile_ncid.rowNdxPos] = ws_counts[pS->WSoutFile_ncid.colNdxPos] = 1;
				ws_coords[pS->WSoutFile_ncid.wsNdxPos] = 0;
				ws_counts[pS->WSoutFile_ncid.wsNdxPos] = NUM_WS_OUTVARS;
				rtn_val = nc_put_vara(pS->WSoutFile_ncid.fileid, pS->WSoutFile_ncid.ws_varid, ws_coords,
						ws_counts, stateData.eqState.cenOutvars);
				assert(chk_nc(rtn_val));
			}
			else for (int i = 0; i<NUM_WS_OUTVARS; i++) if (pS->EQandWSoutputList[i].show) 
			{ float floatval;
				floatval = (float)stateData.eqState.cenOutvars[i];
				rtn_val = nc_put_var1_float(pS->WSoutFile_ncid.fileid, pS->EQandWSoutputList[i].var_id, coords, &floatval); 
				assert(chk_nc(rtn_val));
			}

			// Collect the discretionary output data for this point. 
			if (runParamsP->monthlyOutputFlag && pS->m_num_monthly_vars>0)
			{
				coords[pS->m_moOutFile_ncid.rowNdxPos] = grid_row - pS->m_moOutFile_ncid.row_offset;
				coords[pS->m_moOutFile_ncid.colNdxPos] = grid_col - pS->m_moOutFile_ncid.col_offset;
				pS->writeMonthlyData(coords);
			}
			if (runParamsP->yearlyOutputFlag && pS->m_num_yearly_vars>0)
			{
				coords[pS->m_yrOutFile_ncid.rowNdxPos] = grid_row - pS->m_yrOutFile_ncid.row_offset;
				coords[pS->m_yrOutFile_ncid.colNdxPos] = grid_col - pS->m_yrOutFile_ncid.col_offset;
				pS->writeYearlyData(coords);
			}
			if (runParamsP->multiyrOutputFlag && pS->m_num_multiyr_vars>0)
			{
				coords[pS->m_multiyrOutFile_ncid.rowNdxPos] = grid_row - pS->m_multiyrOutFile_ncid.row_offset;
				coords[pS->m_multiyrOutFile_ncid.colNdxPos] = grid_col - pS->m_multiyrOutFile_ncid.col_offset;
				pS->writeMultiyrData(coords);
			}
			if (runParamsP->timeInvariantOutputFlag && pS->m_num_single_vars>0)
			{
				coords[pS->m_singleOutFile_ncid.rowNdxPos] = grid_row - pS->m_singleOutFile_ncid.row_offset;
				coords[pS->m_singleOutFile_ncid.colNdxPos] = grid_col - pS->m_singleOutFile_ncid.col_offset;
				pS->writeSingleData(coords);
			}

			return(true);  
			break;  

		case TRANSIENT:
			for (int i = 0; i<NUM_WS_OUTVARS; i++) EQdata.cenOutvars[i] = *(input_data.inFromWS + i);

			rtn_flag = (EQdata.cenOutvars[CENinitial_class_mapss]!=NC_FILL_FLOAT) &&
				mc2_transient();
			if (!rtn_flag) return(false);

			if (runParamsP->warmstartOutputFlag)
			{ // Collect the warmstart output data for this point. 
				coords[pS->WSoutFile_ncid.rowNdxPos] = grid_row - pS->WSoutFile_ncid.row_offset;
				coords[pS->WSoutFile_ncid.colNdxPos] = grid_col - pS->WSoutFile_ncid.col_offset;
				if (pS->WSoutFile_ncid.timeNdxPos>=0) coords[pS->WSoutFile_ncid.timeNdxPos] = 0;

				stateData.copyFStoCENoutvars();

				if (modelParamsP->code_flags[FAST_WARMSTART_OUT_FLAG])
				{
					size_t ws_coords[4], ws_counts[4];
					ws_coords[pS->WSoutFile_ncid.timeNdxPos] = 0;
					ws_coords[pS->WSoutFile_ncid.rowNdxPos] = coords[pS->WSoutFile_ncid.rowNdxPos];
					ws_coords[pS->WSoutFile_ncid.colNdxPos] = coords[pS->WSoutFile_ncid.colNdxPos];
					ws_counts[pS->WSoutFile_ncid.timeNdxPos] = ws_counts[pS->WSoutFile_ncid.rowNdxPos] = ws_counts[pS->WSoutFile_ncid.colNdxPos] = 1;
					ws_coords[pS->WSoutFile_ncid.wsNdxPos] = 0;
					ws_counts[pS->WSoutFile_ncid.wsNdxPos] = NUM_WS_OUTVARS;
					rtn_val = nc_put_vara(pS->WSoutFile_ncid.fileid, pS->WSoutFile_ncid.ws_varid, ws_coords,
							ws_counts, stateData.eqState.cenOutvars);
					assert(chk_nc(rtn_val));
				}
				else for (int i = 0; i<NUM_WS_OUTVARS; i++) if (pS->EQandWSoutputList[i].show) 
				{ float floatval;
					floatval = (float)stateData.eqState.cenOutvars[i];
					rtn_val = nc_put_var1_float(pS->WSoutFile_ncid.fileid, pS->EQandWSoutputList[i].var_id, coords, &floatval); 
					if (!chk_nc(rtn_val))
						assert(false);
				}

			} // end of block for warmstart output data

			// Collect the discretionary output data for this point. 
			if (runParamsP->monthlyOutputFlag && pS->m_num_monthly_vars>0)
			{
				coords[pS->m_moOutFile_ncid.rowNdxPos] = grid_row - pS->m_moOutFile_ncid.row_offset;
				coords[pS->m_moOutFile_ncid.colNdxPos] = grid_col - pS->m_moOutFile_ncid.col_offset;
				pS->writeMonthlyData(coords);
			}
			if (runParamsP->yearlyOutputFlag && pS->m_num_yearly_vars>0)
			{
				coords[pS->m_yrOutFile_ncid.rowNdxPos] = grid_row - pS->m_yrOutFile_ncid.row_offset;
				coords[pS->m_yrOutFile_ncid.colNdxPos] = grid_col - pS->m_yrOutFile_ncid.col_offset;
				pS->writeYearlyData(coords);
			}
			if (runParamsP->multiyrOutputFlag && pS->m_num_multiyr_vars>0)
			{
				coords[pS->m_multiyrOutFile_ncid.rowNdxPos] = grid_row - pS->m_multiyrOutFile_ncid.row_offset;
				coords[pS->m_multiyrOutFile_ncid.colNdxPos] = grid_col - pS->m_multiyrOutFile_ncid.col_offset;
				pS->writeMultiyrData(coords);
			}
			if (runParamsP->timeInvariantOutputFlag && pS->m_num_single_vars>0)
			{
				coords[pS->m_singleOutFile_ncid.rowNdxPos] = grid_row - pS->m_singleOutFile_ncid.row_offset;
				coords[pS->m_singleOutFile_ncid.colNdxPos] = grid_col - pS->m_singleOutFile_ncid.col_offset;
				pS->writeSingleData(coords);
			}

			return(true);  
			break;  

		default:
			assert(false);
			break;
	} // end of switch (runParamsP->runMode)

	assert(false); // should never be reached  
	return(false);
} // end of CENTURY_BiogeochemModel::runModelAndCollectData()


bool CENTURY_BiogeochemModel::simulate1pt1yrInner(
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
		float tmax_smoothed[],
		float ppt_raw[],
		float tmp_raw[],
		float tmax_raw[],
		float tmin_raw[],
		float tdmean_raw[],
		float vpr_raw[],
		float rh_raw[])
		// Simulate one point for one year.
{
	bool rtnFlag = true;
	float max_tree_lai, max_grass_lai;
	// float century_carbon_woody[12]; // g C m-2
	// float century_carbon_herb[12]; // g C m-2
	float geochem_lai_woody[12];
	float geochem_lai_herb[12];
	int fire_doy = NO_FIRE;
	float part_burn = 0.;
	float everg, needl;
	bool possibleFireFlag = false;
	logical diagsLogical = runParamsP->diags ? TRUE : FALSE;  
	logical onestep = TRUE;
	float century_c3c4_ratio;
	float century_tmp_index, century_ppt_index;
	int month_of_fire;
	float n_volatil_mo[12], c_ecosys_mo[12], tslc_mo[12], c_live_abovegr_mo[12], c_live_belowgr_mo[12], c_all_abovegr_mo[12];
	float rleavc_mo[12], fbrchc_mo[12], rlwodc_mo[12], crootc_mo[12], frootc_mo[12], adeadc_mo[12], bdeadc_mo[12], metabc_mo[12];
	float strucc_mo[12], litterc_mo[12], wood2c_mo[12], stdedc_mo[12], somsc_mo[12], somtc_mo[12], crmvst_mo[12];
	float gfrac_mo[12], surface_runoff_mo[12], h2o_streamflow_mo[12], pet_mo[12];

	// Set values assuming no fire will be simulated this year.
	// If a fire is simulated, these values will be overwritten.
	pS->VarDict[FIRE].save(0.f, YEAR_INTERVAL);
	pS->VarDict[FIRE_UNSUPPRESSED].save(0.f, YEAR_INTERVAL);
	pS->VarDict[CONSUMED].save(0., YEAR_INTERVAL);
	pS->VarDict[FIRE_KILLED].save(0., YEAR_INTERVAL);
	pS->VarDict[CONSUMED_LIVE].save(0., YEAR_INTERVAL);
	pS->VarDict[CONSUMED_DEAD].save(0., YEAR_INTERVAL);
	pS->VarDict[PART_BURN].save(0., YEAR_INTERVAL);
	pS->VarDict[FIRE_FLI].save(0.f, YEAR_INTERVAL);
	pS->VarDict[FIRE_ROS].save(0.f, YEAR_INTERVAL);
	pS->VarDict[FIRE_ERC].save(0.f, YEAR_INTERVAL);

	if (runParamsP->fireModelFlag) 
	{
		assert(0<=*vtypeP && *vtypeP<=MAX_VTYPE);

		pF->m_ffmc_threshold = modelParamsP->ffmc_threshold_by_vtype[*vtypeP];
		if (pF->m_ffmc_threshold<0.f) 
			pF->m_ffmc_threshold = modelParamsP->ffmc_threshold[m_zone4biogeog][*tree_typP];
		assert(pF->m_ffmc_threshold>0.f);

		pF->m_bui_threshold = modelParamsP->bui_threshold_by_vtype[*vtypeP];
		if (pF->m_bui_threshold<0.f) 
			pF->m_bui_threshold = modelParamsP->bui_threshold[m_zone4biogeog][*tree_typP];
		assert(pF->m_bui_threshold>0.f);

		pF->DFuelMC();
		pS->VarDict[FFMC_ANN_MAX].save(pF->m_ffmc_ann_max, YEAR_INTERVAL);
		pS->VarDict[BUI_ANN_MAX].save(pF->m_bui_ann_max, YEAR_INTERVAL);

		float fire_margin_ann_min = 1.;
		for (int doy = 0; doy<365; doy++)
		{
			float bui_day = pF->m_bui_doy[doy]; if (bui_day<0.f) bui_day = 0.f;
			float bui_frac = bui_day<pF->m_bui_threshold ? (bui_day/pF->m_bui_threshold) : 1.f;
			float ffmc_day = pF->m_ffmc_doy[doy]; if (ffmc_day<0.f) ffmc_day = 0.f;
			float ffmc_frac = ffmc_day<pF->m_ffmc_threshold ? (ffmc_day/pF->m_ffmc_threshold) : 1.f;
			float fire_margin = ffmc_frac<bui_frac ? (1.f - ffmc_frac) : (1.f - bui_frac);
			if (fire_margin<fire_margin_ann_min) fire_margin_ann_min = fire_margin;
			pF->m_fire_margin_doy[doy] = fire_margin;
		}
		pS->VarDict[FIRE_MARGIN_ANN_MIN].save(fire_margin_ann_min, YEAR_INTERVAL);

		// Is a fire possible this year?
		bool long_enough_since_last_fireFlag = true; // Assume yes.
		if (!modelParamsP->part_burnFlag)
		{ // partial cell burn mechanism is turned off; Get min fire return interval for current veg class
			float min_fri, max_fri; 
			pF->fire_return_interval(*vtypeP, &min_fri, &max_fri);
			long_enough_since_last_fireFlag = m_burn_count>=min_fri;
		}
		possibleFireFlag = (runParamsP->fire_set_interval>0 && runParamsP->fire_set_interval<=m_burn_count)
			|| (fire_margin_ann_min==0.f /* pF->m_ffmc_ann_max>=pF->m_ffmc_threshold && pF->m_bui_ann_max>=pF->m_bui_threshold */
					&& long_enough_since_last_fireFlag);               
	}
	else
	{
		pS->VarDict[FFMC_ANN_MAX].save(NC_FILL_FLOAT, YEAR_INTERVAL);
		pS->VarDict[BUI_ANN_MAX].save(NC_FILL_FLOAT, YEAR_INTERVAL);
		pS->VarDict[FIRE_MARGIN_ANN_MIN].save(NC_FILL_FLOAT, YEAR_INTERVAL);
	} // end of if (runParamsP->fireModelFlag) ... else ...

	month_of_fire = 0;
	int days_per_month[] = DAYS_PER_MONTH;
	int first_doy_of_next_month = 0;

	float npp_ytd = 0;
	for (int mo = JAN; mo <= DEC; mo++) 
	{ 
		int first_doy_of_month = first_doy_of_next_month;
		int n_days = days_per_month[modelParamsP->southernHemisphereFlag ? (mo + 6)%12 : mo]; 
		first_doy_of_next_month = first_doy_of_month + n_days;

		int woody_begin_grow, woody_end_grow, grass_begin_grow, grass_end_grow, grass_senescence;
		float century_carbon_tree, century_carbon_grass; // g C m-2
		float century_mix_index;
		char century_woody_fire_name[6] = "";
		char century_grass_fire_name[6] = "";
		char unusedCharStr[6] = "";
		logical warmstart_init_flag;
		logical burn_year;
		float global_burn_count;
		float time_in = 0.9999;
		logical unlimited_N_flag = modelParamsP->unlimitedNflag ? TRUE : FALSE;
		float clai;

		// The next line reflects that warmstart_init_flag is a (FORTRAN) logical, not a bool.
		warmstart_init_flag = (runParamsP->runMode==TRANSIENT && yrNdx==0 && mo==0) ? TRUE : FALSE; 
		burn_year = FALSE;
		global_burn_count = m_burn_count;           
		cen_state[187] = modelParamsP->m_century_runoff_x_intercept;
		cen_state[188] = modelParamsP->m_century_runoff_slope;
		cen_state[189] = modelParamsP->m_lait_lower_limit;

		modelParamsP->code_flags[ALT_FUEL_LOAD_CODE_FLAG] = modelParamsP->altFuelLoadFlag;
		modelParamsP->code_flags[ALT_TREE_ALLOMETRY_CODE_FLAG] = modelParamsP->altTreeAllometryFlag;
		for (int i=190; i<200; i++) cen_state[i] = modelParamsP->code_flags[i - 190] ? 1 : 0;

		if (m_fire_last_month)
		{
			strcpy(century_woody_fire_name, "BURN ");
			strcpy(century_grass_fire_name, "H    ");
		}
		cen_step_(
				prev_vclassP,
				&woody_begin_grow, &woody_end_grow,		/* result */
				&grass_begin_grow, &grass_end_grow, &grass_senescence,	/* result */
				century_woody_fire_name,
				century_grass_fire_name,
				&century_carbon_tree, 	/* result */
				&century_carbon_grass, /* result */
				cen_state,
				EQdata.cenOutvars,
				&mix_index,
				&century_mix_index,		/* result */
				&m_c3pct,
				&century_c3c4_ratio, 	/* result */
				unusedCharStr, unusedCharStr,
				&diagsLogical, &onestep,
				&burn_year,				/* result */
				&global_burn_count,		/* result */
				&tmp_index,
				&century_tmp_index,		/* result */
				&ppt_index,
				&century_ppt_index,		/* result */
				&mapss_climate_zone,	
				&warmstart_init_flag,			/* input */
				&time_in,				/* input */
				&runParamsP->actual_years_to_run,				/* input */
				&unlimited_N_flag,
				&m_frost_index_this_year);

		if (m_fire_last_month)
		{
			strcpy(century_woody_fire_name, "");
			strcpy(century_grass_fire_name, "");
			m_fire_last_month = false;
		}
		m_burn_count = global_burn_count;
		// century_carbon_woody[mo] = century_carbon_tree;
		// century_carbon_herb[mo] = century_carbon_grass; 
		geochem_lai_herb[mo] = century_carbon_grass*2.2*0.01; // 0.01 is specific leaf area "from lit."; 2.2 is "grass factor"

		m_npp_mo[mo] = cen_state[23] - npp_ytd;
		npp_ytd = cen_state[23];
		pS->VarDict[NPP].save(m_npp_mo[mo], MONTH_INTERVAL);
		m_aet_mo[mo] = cen_state[26]; 
		pS->VarDict[AET].save(m_aet_mo[mo], MONTH_INTERVAL);
		m_vegc_mo[mo] = cen_state[28];
		clai = geochem_lai_woody[mo] = cen_state[35];
		m_frstc_mo[mo] = cen_state[39];
		frootc_mo[mo] = cen_state[45];
		crootc_mo[mo] = cen_state[46];
		m_bflivc_mo[mo] = frootc_mo[mo] + crootc_mo[mo];
		m_aflivc_mo[mo] = m_frstc_mo[mo] - m_bflivc_mo[mo];
		m_aglivc_mo[mo] = cen_state[40];
		m_bglivc_mo[mo] = cen_state[41];
		c_live_abovegr_mo[mo] = m_aflivc_mo[mo] + m_aglivc_mo[mo];
		c_live_belowgr_mo[mo] = m_bflivc_mo[mo] + m_bglivc_mo[mo];
		c_all_abovegr_mo[mo] = m_vegc_mo[mo] - (m_bglivc_mo[mo] + frootc_mo[mo] + crootc_mo[mo]);
		gfrac_mo[mo] = m_vegc_mo[mo]>0. ? (m_aglivc_mo[mo] + m_bglivc_mo[mo])/m_vegc_mo[mo] : 0.0;
		m_fprd_ppt_mo[mo] = cen_state[50];
		pS->VarDict[FPRD_PPT].save(m_fprd_ppt_mo[mo], MONTH_INTERVAL);
		m_fprd_tmp_mo[mo] = cen_state[51];
		pS->VarDict[FPRD_TMP].save(m_fprd_tmp_mo[mo], MONTH_INTERVAL);
		m_gprd_ppt_mo[mo] = cen_state[52];
		pS->VarDict[GPRD_PPT].save(m_gprd_ppt_mo[mo], MONTH_INTERVAL);
		m_gprd_tmp_mo[mo] = cen_state[53];
		pS->VarDict[GPRD_TMP].save(m_gprd_tmp_mo[mo], MONTH_INTERVAL);
		m_ddecid_mo[mo] = cen_state[58];
		n_volatil_mo[mo] = EQdata.cenOutvars[CENvolgma] + EQdata.cenOutvars[CENvolexa] + EQdata.cenOutvars[CENvolpla];
		pS->VarDict[N_VOLATIL].save(n_volatil_mo[mo], MONTH_INTERVAL);
		tslc_mo[mo] = cen_state[29];
		c_ecosys_mo[mo] = m_vegc_mo[mo] + tslc_mo[mo];
		rleavc_mo[mo] = EQdata.cenOutvars[CENrleavc];
		fbrchc_mo[mo] = EQdata.cenOutvars[CENfbrchc];
		rlwodc_mo[mo] = EQdata.cenOutvars[CENrlwodc];
		crootc_mo[mo] = EQdata.cenOutvars[CENcrootc];
		frootc_mo[mo] = EQdata.cenOutvars[CENfrootc];
		adeadc_mo[mo] = cen_state[47];
		bdeadc_mo[mo] = cen_state[48];
		metabc_mo[mo] = EQdata.cenOutvars[CENmetabc_1];
		strucc_mo[mo] = EQdata.cenOutvars[CENstrucc_1];
		litterc_mo[mo] = metabc_mo[mo] + strucc_mo[mo];
		wood2c_mo[mo] = EQdata.cenOutvars[CENwood2c];
		stdedc_mo[mo] = EQdata.cenOutvars[CENstdedc];
		somsc_mo[mo] = EQdata.cenOutvars[CENsomsc];
		somtc_mo[mo] = EQdata.cenOutvars[CENsomtc];
		crmvst_mo[mo] = EQdata.cenOutvars[CENcrmvst];
		pS->VarDict[C_HARVEST].save(crmvst_mo[mo], MONTH_INTERVAL);
		surface_runoff_mo[mo] = EQdata.cenOutvars[CENsurface_runoff]; // ??? x 10 to convert cmH2O to mmH2O?
		pS->VarDict[SFC_RUNOFF].save(surface_runoff_mo[mo], MONTH_INTERVAL);
		h2o_streamflow_mo[mo] = EQdata.cenOutvars[CENstream_1];  // ??? x 10 to convert cmH2O to mmH2O?
		pS->VarDict[H2O_STREAM_FLOW].save(h2o_streamflow_mo[mo], MONTH_INTERVAL);
		pet_mo[mo] = EQdata.cenOutvars[CENpet];  // ??? x 10 to convert cmH2O to mmH2O?
		pS->VarDict[PET].save(pet_mo[mo], MONTH_INTERVAL);

		float fire_margin_month = NC_FILL_FLOAT;
		if (runParamsP->fireModelFlag) 
		{ 
			float tree_stress, grass_stress;
			int fire_dom, month_length;

			fire_margin_month = 1.f;     
			for (int doy = first_doy_of_month; doy < first_doy_of_next_month; doy++)
				if (fire_margin_month>pF->m_fire_margin_doy[doy]) fire_margin_month = pF->m_fire_margin_doy[doy];

			pF->m_mc_duff = cen_state[20] * 100.;  
			grass_stress = 1. -  cen_state[20];
			pF->m_mc_grass = pF->liveFuelMoistureContent(grass_stress, pF->fireParams.mc_grass_min, pF->fireParams.mc_grass_max);
			tree_stress = 1. - cen_state[19];
			pF->m_mc_tree = pF->liveFuelMoistureContent(tree_stress, pF->fireParams.mc_tree_min, pF->fireParams.mc_tree_max);

			/* make var assigments from CENTURY output and convert from g C m-2 to g DM m-2 */
			pF->m_lgras   = EQdata.cenOutvars[3]*2.5f;
			pF->m_lleaf   = EQdata.cenOutvars[445]*2.5f;
			pF->m_lwod1   = EQdata.cenOutvars[419]*2.0f;
			pF->m_lwod100 = EQdata.cenOutvars[452]*2.0f;
			pF->m_dstnd   = EQdata.cenOutvars[262]*2.0f;
			pF->m_mlittr  = EQdata.cenOutvars[104]*2.0f;
			pF->m_slittr  = EQdata.cenOutvars[251]*2.0f;
			pF->m_dwod1   = EQdata.cenOutvars[482]*2.0f;
			pF->m_dwod100 = EQdata.cenOutvars[483]*2.0f;

			month_length = sciFn.days_per_mo[modelParamsP->southernHemisphereFlag ? (mo + 6)%12 : mo];
			pF->FuelLoad(*vtypeP, clai, month_length);
			// pF->erc_G(month_length);

			if (mo==6) pF->stateData.fireState.clai_in_midseason = clai;

			if (mo==pF->m_month_of_min_mc) 
			{
				pF->m_fine_dead_fuel_yr = pF->m_d1hr + pF->m_d10hr;
				pF->m_coarse_dead_fuel_yr = pF->m_d100hr + pF->m_d1000hr;
			}

			if (fire_doy==NO_FIRE)
			{ // no fire so far this year; should a fire be simulated this month?
				bool not_enough_C = modelParamsP->code_flags[FIRE_CODE_FLAG] ?
					(m_vegc_mo[mo]<=105.) : (c_all_abovegr_mo[mo]<60.);
				if (!possibleFireFlag || fire_margin_month>0.f || not_enough_C) fire_dom = NO_FIRE;
				else pF->FireOccur(mo, m_burn_count, *vtypeP, &fire_dom, &fire_doy);
				if (fire_dom!=NO_FIRE) 
				{ // Simulate a fire.
					bool local_fire_suppression_flag = runParamsP->fireSuppressionFlag &&
						(runParamsP->first_calendar_year_of_run + yrNdx)>=runParamsP->fire_suppression_first_year;
					month_of_fire = mo + 1;
					part_burn = pF->part_burn(fire_doy, fire_dom, m_burn_count, &local_fire_suppression_flag, *vtypeP);
					pS->VarDict[PART_BURN].save(part_burn, YEAR_INTERVAL);
					if (!local_fire_suppression_flag) 
					{
						pS->VarDict[FIRE_UNSUPPRESSED].save(1.f, YEAR_INTERVAL);
						pS->VarDict[FIRE_FLI].save(pF->m_fire_fli, YEAR_INTERVAL);
						pS->VarDict[FIRE_ROS].save(pF->m_fire_ros, YEAR_INTERVAL);
						pS->VarDict[FIRE_ERC].save(pF->m_fire_erc, YEAR_INTERVAL);
					}
					m_fire_last_month = true; // cen_step() will simulate a fire next month

					pF->FireEffect(fire_dom, fire_doy, *vtypeP);

					cen_state[0]  = pF->m_consumed[1]/100.f * part_burn; // frac live leaf consumed REMF(1)
					cen_state[1]  = pF->m_consumed[2]/100.f * part_burn; // frac live fine branch consumed REMF(2)
					cen_state[2]  = 0.; // live stems are not consumed REMF(3)
					cen_state[3]  = pF->m_consumed[3]/100.f * part_burn; // frac dead fine branch consumed REMF(4)
					cen_state[4]  = pF->m_consumed[5]/100.f * part_burn; // frac dead large wood consumed REMF(5)
					cen_state[5]  = pF->m_consumed[0]/100.f * part_burn; // frac live grass consumed FLFREM
					cen_state[6]  = pF->m_consumed[3]/100.f * part_burn; // frac dead standing fuel consumed FDREM(1)
					cen_state[7]  = pF->m_consumed[3]/100.f * part_burn; // frac litter consumed FDREM(2)

					cen_state[8]  = pF->m_killed[0]/100.f * part_burn; // fire-induced live leaf turnover LEAFDR(MO)
					cen_state[9]  = pF->m_killed[1]/100.f * part_burn; // fire-induced live fine branch turnover WOODDR(3)
					cen_state[10] = pF->m_killed[2]/100.f * part_burn; // fire-induced coarse wood turnover WOODDR(4)
					cen_state[11] = pF->m_killed[3]/100.f * part_burn; // frac live fine roots killed FD(1)
					cen_state[12] = pF->m_killed[3]/100.f * part_burn; // frac live coarse roots killed FD(2)

					// Limit effects to 99.9%
					for (int i = 0; i<13; i++) if (cen_state[i]>0.999) cen_state[i] = 0.999; 

					cen_state[13] = TRUE; // fire signal DID_BURN

					pF->ann_effect(EQdata.cenOutvars, cen_state, part_burn);

					pS->VarDict[FIRE].save(1.f, YEAR_INTERVAL);
					pS->VarDict[CONSUMED].save(pF->m_consume_totbio, YEAR_INTERVAL);
					pS->VarDict[FIRE_KILLED].save(pF->m_death_totbio, YEAR_INTERVAL);
					pS->VarDict[CONSUMED_LIVE].save(pF->m_consume_live, YEAR_INTERVAL);
					pS->VarDict[CONSUMED_DEAD].save(pF->m_consume_dead, YEAR_INTERVAL);
				} // end of block to simulate a fire (if (fire_dom!=NO_FIRE))
			} // end of block to decide whether a fire occurs this month (if (fire_doy==NO_FIRE))

			pF->stateData.fireState.prev_l1hr   = pF->m_lgras;
			pF->stateData.fireState.prev_1hr    = pF->m_d1hr;
			pF->stateData.fireState.prev_10hr   = pF->m_d10hr;
			pF->stateData.fireState.prev_100hr  = pF->m_d100hr;
			pF->stateData.fireState.prev_1000hr = pF->m_d1000hr;
			pF->stateData.fireState.prev_depth  = pF->m_fuel_depth;
		} // end of if (runParamsP->fireModelFlag)

		pS->VarDict[PPT].save(ppt_raw[mo], MONTH_INTERVAL);
		pS->VarDict[TMP].save(tmp_raw[mo], MONTH_INTERVAL);
		pS->VarDict[TMAX].save(tmax_raw[mo], MONTH_INTERVAL);
		pS->VarDict[TMAX_SMOOTHED].save(tmax_smoothed[mo], MONTH_INTERVAL);
		pS->VarDict[TMIN].save(tmin_raw[mo], MONTH_INTERVAL);
		pS->VarDict[TDMEAN].save(tdmean_raw[mo], MONTH_INTERVAL);
		pS->VarDict[VPR].save(vpr_raw[mo], MONTH_INTERVAL);
		pS->VarDict[RH].save(rh_raw[mo], MONTH_INTERVAL);      
		pS->VarDict[PPT_SMOOTHED].save(ppt_smoothed[mo], MONTH_INTERVAL);      
		pS->VarDict[FIRE_MARGIN].save(fire_margin_month, MONTH_INTERVAL);

		pS->saveOneMonthOneCellOutputData(12*yrNdx + mo);

	} // end of month loop

	m_npp_yr = cen_state[23];
	m_rsp_yr = cen_state[24];
	m_nep_yr = cen_state[25];
	m_nidx = cen_state[33];
	m_eidx = cen_state[34];
	m_max_gfrac = sciFn.annual_max(gfrac_mo);
	pS->VarDict[NPP_TREE].save(cen_state[36], YEAR_INTERVAL);
	pS->VarDict[NPP_GRASS].save(cen_state[37] + cen_state[38], YEAR_INTERVAL); // agcaccx + bgcaccx
	m_bio_consume_century = cen_state[56];
	m_nbp_yr = m_nep_yr - m_bio_consume_century;
	m_nlayer = cen_state[59];
	pS->VarDict[AET].save(sciFn.annual_sum(m_aet_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[PET].save(sciFn.annual_sum(pet_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[N_VOLATIL].save(sciFn.annual_sum(n_volatil_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_ECOSYS].save(sciFn.annual_average(c_ecosys_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_ECOSYS_DEC].save(c_ecosys_mo[11], YEAR_INTERVAL);
	pS->VarDict[C_LIVE_ABOVEGR].save(sciFn.annual_average(c_live_abovegr_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_LIVE_BELOWGR].save(sciFn.annual_average(c_live_belowgr_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_MAX_FOREST_LEAF].save(sciFn.annual_max(rleavc_mo), YEAR_INTERVAL);
	pS->VarDict[C_FINE_BRANCH].save(sciFn.annual_average(fbrchc_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_BOLE].save(sciFn.annual_average(rlwodc_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_MAX_COARSE_ROOT].save(sciFn.annual_max(crootc_mo), YEAR_INTERVAL);
	pS->VarDict[C_MAX_FINE_ROOT].save(sciFn.annual_max(frootc_mo), YEAR_INTERVAL);
	pS->VarDict[C_MAX_LIVE_GRASS_BELOWGR].save(sciFn.annual_max(m_bglivc_mo), YEAR_INTERVAL);
	pS->VarDict[C_DEAD_ABOVEGR].save(sciFn.annual_average(adeadc_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_DEAD_BELOWGR].save(sciFn.annual_average(bdeadc_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_LITTER_METAB].save(sciFn.annual_average(metabc_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_LITTER_STRUC].save(sciFn.annual_average(strucc_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_LITTER].save(sciFn.annual_average(litterc_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_DEAD_WOOD].save(sciFn.annual_average(wood2c_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_MAX_STANDING_DEAD].save(sciFn.annual_max(stdedc_mo), YEAR_INTERVAL);
	pS->VarDict[C_SOIL_AND_LITTER].save(sciFn.annual_average(tslc_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_SOM_X_STRUC_METAB].save(sciFn.annual_average(somsc_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_SOM].save(sciFn.annual_average(somtc_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[TREE_HT].save(pF->m_tree_ht_m, YEAR_INTERVAL);
	pS->VarDict[D1HR].save(pF->m_d1hr, YEAR_INTERVAL);
	pS->VarDict[D10HR].save(pF->m_d10hr, YEAR_INTERVAL);
	pS->VarDict[D100HR].save(pF->m_d100hr, YEAR_INTERVAL);
	pS->VarDict[D1000HR].save(pF->m_d1000hr, YEAR_INTERVAL);
	pS->VarDict[EM_CO].save(pF->m_em_co, YEAR_INTERVAL);
	pS->VarDict[EM_CO2].save(pF->m_em_co2, YEAR_INTERVAL);
	pS->VarDict[EM_CH4].save(pF->m_em_ch4, YEAR_INTERVAL);
	pS->VarDict[EM_NMHC].save(pF->m_em_nmhc, YEAR_INTERVAL);
	pS->VarDict[EM_PM].save(pF->m_em_pm, YEAR_INTERVAL);
	pS->VarDict[C_GRAIN].save(EQdata.cenOutvars[CENcgrain], YEAR_INTERVAL);
	pS->VarDict[C_HARVEST].save(sciFn.annual_average(crmvst_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[TMP_INDEX].save(century_tmp_index, YEAR_INTERVAL);
	pS->VarDict[PPT_INDEX].save(century_ppt_index, YEAR_INTERVAL);
	pS->VarDict[NEEDLE_INDEX].save(m_nidx, YEAR_INTERVAL);
	pS->VarDict[EVERGREEN_INDEX].save(m_eidx, YEAR_INTERVAL);
	pS->VarDict[FPRD_PPT].save(sciFn.annual_average(m_fprd_ppt_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[FPRD_TMP].save(sciFn.annual_average(m_fprd_tmp_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[GPRD_PPT].save(sciFn.annual_average(m_gprd_ppt_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[GPRD_TMP].save(sciFn.annual_average(m_gprd_tmp_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[SFC_RUNOFF].save(sciFn.annual_sum(surface_runoff_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[MAX_SFC_RUNOFF].save(sciFn.annual_max(surface_runoff_mo), YEAR_INTERVAL);
	pS->VarDict[H2O_STREAM_FLOW].save(sciFn.annual_sum(h2o_streamflow_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	// pS->VarDict[].save(sciFn.annual_average(, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);

	sciFn.PhytomorphologicalIndices(century_tmp_index, century_ppt_index, &everg, &needl);
	assert(sciFn.close_enough(needl, m_nidx, 0.0001, 0.001)
			&& sciFn.close_enough(everg, m_eidx, 0.0001, 0.001));
	m_c3pct = century_c3c4_ratio;               

	max_tree_lai = sciFn.annual_max(geochem_lai_woody);
	max_grass_lai = sciFn.annual_max(geochem_lai_herb);
	pS->VarDict[MAX_TREE_LAI].save(max_tree_lai, YEAR_INTERVAL);
	pS->VarDict[MAX_GRASS_LAI].save(max_grass_lai, YEAR_INTERVAL);
	if (runParamsP->runMode==SPINUP && yrNdx==0)
	{
		m_smoothed_max_tree_lai = max_tree_lai;
		m_smoothed_max_grass_lai = max_grass_lai;
	}
	else 
	{ float tau_bio;
		tau_bio = MIN(m_burn_count, modelParamsP->efold_t_max);
		m_smoothed_max_tree_lai = sciFn.efold(tau_bio, max_tree_lai, m_smoothed_max_tree_lai);
		m_smoothed_max_grass_lai = sciFn.efold(tau_bio, max_grass_lai, m_smoothed_max_grass_lai);
	}

	// biogeography calculations
	BiogeographyInputData inVals;
	bool broadleafFlag, deciduousFlag;
	float grassc_mo[12], treec_mo[12];

	for (int mo = 0; mo < 12; mo++)
	{ 
		grassc_mo[mo] = m_aglivc_mo[mo] + m_bglivc_mo[mo];
		treec_mo[mo] = m_aflivc_mo[mo] + m_bflivc_mo[mo];
	}  

	treeType(everg, needl, sciFn.annual_sum(ppt_smoothed, modelParamsP->southernHemisphereFlag), 
			runParamsP->baseCalibration, tree_typP, &broadleafFlag, &deciduousFlag);
	inVals.tree_typ = *tree_typP;
	inVals.npp_yr = m_npp_yr;
	inVals.mean_treec = sciFn.annual_average(treec_mo, modelParamsP->southernHemisphereFlag);
	inVals.max_grassc = sciFn.annual_max(grassc_mo); 
	inVals.mean_vegc = sciFn.annual_average(m_vegc_mo, modelParamsP->southernHemisphereFlag); 
	inVals.max_grass_frac = m_max_gfrac;
	inVals.zone = m_zone4biogeog;
	inVals.c3pct = m_c3pct;
	inVals.gdd0 = m_gdd_zero;
	inVals.cont_index = m_cont_index;
	inVals.min_smoothed_tmp = m_min_smoothed_tmp;
	inVals.ppt_yr = sciFn.annual_sum(ppt_smoothed, modelParamsP->southernHemisphereFlag);
	inVals.tmp_yr = sciFn.annual_average(tmp_smoothed, modelParamsP->southernHemisphereFlag);
	inVals.elev = inputDataP->elev;
	inVals.tot_summer_ppt = ppt_smoothed[5] + ppt_smoothed[6] + ppt_smoothed[7]; // JJA in N hemi, 92 days
	if (modelParamsP->southernHemisphereFlag) inVals.tot_summer_ppt *= 92.f/90.f; // DJF in S hemi, 90 days

	MC2VegType mc2vtype;
	switch (runParamsP->baseCalibration)
	{
		case mc2W_WA:
			inVals.fog = inputDataP->fog; // if (inVals.fog<-100.f) return(false);
			inVals.topomoist = inputDataP->topomoist; // if (inVals.topomoist<-100.f) return(false);
			inVals.deltaTsl = inputDataP->deltaTsl; // if (inVals.deltaTsl<-100.f) return(false);
			inVals.aspect = inputDataP->aspect; // if (inVals.aspect<-100.f) return(false);
			inVals.sw = inputDataP->sw; // if (inVals.sw<-100.f) return(false);
			inVals.cad = inputDataP->cad; // if (inVals.cad<-100.f) return(false);
			inVals.subregion = pS->m_maskVal;
			// fall into mc2ConUS case
		case mc2ConUS:
		case mc2GLOBAL:
		case mc2California:
		case mc2BlueMtns:
			{ // Don't change veg type in the year of a fire, or within regrowth_years of the last fire.
				bool vtype_change_allowed_flag = !runParamsP->fireModelFlag ||
					modelParamsP->code_flags[REGRESSION_TEST_2B101_FLAG] ||
					m_burn_count>modelParamsP->regrowth_years;
				if (vtype_change_allowed_flag) 
				{ 
					rtnFlag = pB->BiogeogMC2(inVals, &m_biome, &m_physiognomic_class, &mc2vtype);
					*vtypeP = (int)mc2vtype;
				}
			}
			break;
		case mc2ConUS_LC:
			inVals.lulcType = Lynx_Undefined;
			inVals.ddecid_mo = m_ddecid_mo;
			inVals.tree_lai = max_tree_lai;
			inVals.grass_lai = max_grass_lai;
			inVals.mean_grass_frac = sciFn.annual_average(gfrac_mo, modelParamsP->southernHemisphereFlag);
			inVals.tmmin = m_tmmin;
			inVals.nlayer = m_nlayer;
			inVals.aflivc = sciFn.annual_average(m_aflivc_mo, modelParamsP->southernHemisphereFlag);
			// Don't change veg type in the year of a fire.
			if (fire_doy==NO_FIRE || modelParamsP->code_flags[REGRESSION_TEST_2B101_FLAG]) 
			{ 
				mc2vtype = pB->BiogeogLC(inVals);
				*vtypeP = (int)mc2vtype;
			}
			m_biome = UNKNOWNbiome;
			m_physiognomic_class = UNKNOWNpclass;
			break;
		default: assert(0); break;
	} // end of switch (runParamsP->baseCalibration)
	// end of biogeography calculations

	if (!rtnFlag) return(false);

	if (pS->VarDict[GROUSE_HABITAT].show) 
	{
		int grouse_habitat_type; 
		float augmaxt = tmax_smoothed[7]; // August tmax in N hemisphere; Feb tmax in S hemisphere
		float smrpre = ppt_smoothed[5] + ppt_smoothed[6] + ppt_smoothed[7]; // JJA precip in N hemi; DJF in S hemi
		float anntmp = sciFn.annual_average(tmp_smoothed, modelParamsP->southernHemisphereFlag);
		grouse_habitat_type = pB->GrouseModel(*vtypeP, smrpre, augmaxt, anntmp);
		pS->VarDict[GROUSE_HABITAT].save(grouse_habitat_type, YEAR_INTERVAL);
		pS->VarDict[GROUSE_SMRPRE].save(smrpre, YEAR_INTERVAL);
		pS->VarDict[GROUSE_AUGMAXT].save(augmaxt, YEAR_INTERVAL);
		pS->VarDict[GROUSE_ANNTMP].save(anntmp, YEAR_INTERVAL);
	}

	if (runParamsP->fireModelFlag)
	{
		pF->stateData.fireState.prev_mc_grass = pF->m_mc_grass;
		pF->stateData.fireState.prev_mc_tree = pF->m_mc_tree;
	}
	pS->VarDict[MONTH_OF_FIRE].save(month_of_fire, YEAR_INTERVAL);
	pS->VarDict[MINERL_5].save(EQdata.cenOutvars[120], YEAR_INTERVAL);
	pS->VarDict[C_MAX_LIVE_GRASS].save(inVals.max_grassc, YEAR_INTERVAL);         
	return(rtnFlag);
} // end of CENTURY_BiogeochemModel::simulate1pt1yrInner()


///////////////////////////////// mc2_transient() //////////////////////////////////////////////////

bool CENTURY_BiogeochemModel::mc2_transient()
{
	Soil_MAPSS soils_data(inputDataP->soilData);   
	SoilsType	WhichSoils = ScsSoilsData;
	bool rtnFlag = true; 
	VCategory initial_vclass;
	MAPSSvegClass initial_class;
	int initial_vtype, vtype;
	float initial_mix_index, initial_c3c4_ratio, initial_tmp_index, initial_ppt_index;
	float initial_everg, initial_needl;
	TreeType initial_tree_typ, tree_typ;
	float cen_state[CEN_STATE_SIZE];
	VCategory prev_vclass;
	float ppt_smoothed[12], tmp_smoothed[12], tmax_smoothed[12], tmin_smoothed[12], tdmean_smoothed[12];
	logical diagsLogical = runParamsP->diags ? TRUE : FALSE;  

	pS->VarDict[ELEV].save(inputDataP->elev);
	pS->VarDict[SWHC_TOP].save(soils_data.swhc[0]);
	pS->VarDict[SWHC_MID].save(soils_data.swhc[1]);
	pS->VarDict[SWHC_DEEP].save(soils_data.swhc[2]);
	pS->VarDict[SWHC_ALL].save(soils_data.swhc[0] + soils_data.swhc[1] + soils_data.swhc[2]);
	pS->VarDict[SOIL_DEPTH].save(soils_data.mineral_depth);

	/* Note that the EQdata.cenOutVars has been initialized from a warmstart file */
	// Get initial conditions from EQdata.  
	for (int mo = 0; mo<12; mo++)
	{
		ppt_smoothed[mo] = EQdata.cenOutvars[CENppt_prev_0 + mo];
		tmp_smoothed[mo] = EQdata.cenOutvars[CENtmp_prev_0 + mo];

		int source_mo = modelParamsP->southernHemisphereFlag ? (mo + 6) : mo;
		tmax_smoothed[mo] = EQdata.cenOutvars[CENtmax_prev_0 + mo]; 
		if (tmax_smoothed[mo]==NC_FILL_FLOAT) tmax_smoothed[mo] = *(inputDataP->tmaxP + source_mo);
		tmin_smoothed[mo] = EQdata.cenOutvars[CENtmin_prev_0 + mo];
		if (tmin_smoothed[mo]==NC_FILL_FLOAT) tmin_smoothed[mo] = *(inputDataP->tminP + source_mo);      
		tdmean_smoothed[mo] = EQdata.cenOutvars[CENtdmean_prev_0 + mo];
		if (tdmean_smoothed[mo]==NC_FILL_FLOAT) tdmean_smoothed[mo] = *(inputDataP->tdmeanP + source_mo);
	}    

	m_frost_index_this_year = EQdata.cenOutvars[CENfrost_index_this_year]; // maybe this should be the frost index of the current year
	initial_tmp_index = EQdata.cenOutvars[CENtmp_index];
	initial_ppt_index = EQdata.cenOutvars[CENppt_index];
	initial_vtype = (int)EQdata.cenOutvars[BSvtype];
	m_fire_last_month = EQdata.cenOutvars[CENfire_last_month]!=0.;

	sciFn.PhytomorphologicalIndices(initial_tmp_index, initial_ppt_index, &initial_everg, &initial_needl);
	treeType(initial_everg, initial_needl, 
			sciFn.annual_sum(ppt_smoothed, modelParamsP->southernHemisphereFlag),
			runParamsP->baseCalibration, &initial_tree_typ, &stateData.deciduous, &stateData.broadleaf);

	cen_zero_(); /* Zero out the labeled commons in century. */

	/* Pass soils info to Century. */
	float sand, clay, depth, rock[NSL];
	sand = (float) soils_data.sand[SURFACE];
	clay = (float) soils_data.clay[SURFACE];
	depth = (float) soils_data.mineral_depth;
	if (depth<=10.) return(false); // If depth is <= 10 mm, quit.
	for (int i = SURFACE; i <= DEEP; i++) rock[i] = (float) soils_data.rock_frag_mineral[i];  
	float hdeposition = 0.;
	float cdeposition = 0.;
	cen_init_soils_(
			&soils_data.bulk_density,
			&sand, 
			&clay,
			&depth,
			rock,
			&WhichSoils,
			&hdeposition,
			&cdeposition);

	/* Pass latitude and elevation to Century. */
	float latitude = (float)inputDataP->lat;
	cen_init_lat_(&latitude, &inputDataP->elev);

	/* Zero out annet for cen_init_() in order to get wdfxs right in eachyr.F */
	// EQdata.cenOutvars[549] = 0.0; only in spinup, not in transient

	int fire_mode = runParamsP->fireModelFlag ? 1 : 0;
	int unlimitedN = modelParamsP->unlimitedNflag ? 1 : 0;

	MC_BiogeographyModel biogeogModelInstance(pS, runParamsP, modelParamsP);
	// MC_FireModel fmi(m_pRun, &stateData.fireState, inputDataP); // "fmi" is "fireModelInstance"
	MC_FireModel fmi(pS, EQdata.cenOutvars, inputDataP); // "fmi" is "fireModelInstance"
	if (runParamsP->fireModelFlag) fmi.portable_srand(fmi.stateData.fireState.rand_seed);

	m_smoothed_max_tree_lai = EQdata.cenOutvars[BSsmoothed_max_tree_lai];
	m_smoothed_max_grass_lai = EQdata.cenOutvars[BSsmoothed_max_grass_lai];
	m_burn_count = EQdata.cenOutvars[FSyrs_since_fire]; 

	pS->VarDict[FFMC_ANN_MAX].save(0., YEAR_INTERVAL);

	// first year initialization and cen_init() call
	initial_mix_index = EQdata.cenOutvars[CENmx_index];
	initial_tmp_index = EQdata.cenOutvars[CENtmp_index];
	initial_ppt_index = EQdata.cenOutvars[CENppt_index];
	initial_c3c4_ratio = EQdata.cenOutvars[CENc3c4_index];
	int initial_mapss_climate_zone = EQdata.cenOutvars[CENmapss_zone];

	initial_class = (MAPSSvegClass)EQdata.cenOutvars[CENinitial_class_mapss];
	initial_vclass = (VCategory)EQdata.cenOutvars[CENinitial_vclass_mapss];
	prev_vclass = initial_vclass;
	if (initial_class<=NaturalBarren || initial_class==Ice || initial_class==DesertExtreme || initial_vclass==VUnknown) 
		return(false);  

	if (runParamsP->diags) 
	{
		printf("class, vclass: %d, %d \n", initial_class, initial_vclass);
		printf("tmp_index: %.2f \n", initial_tmp_index);
		printf("ppt_index: %.2f \n", initial_ppt_index);
		printf("c3c4_ratio: %.2f \n", initial_c3c4_ratio);
	}

	for (int i=0; i<187; i++) cen_state[i] = 0.0;
	cen_state[187] = modelParamsP->m_century_runoff_x_intercept;
	cen_state[188] = modelParamsP->m_century_runoff_slope;
	cen_state[189] = modelParamsP->m_lait_lower_limit;

	modelParamsP->code_flags[ALT_FUEL_LOAD_CODE_FLAG] = modelParamsP->altFuelLoadFlag;
	modelParamsP->code_flags[ALT_TREE_ALLOMETRY_CODE_FLAG] = modelParamsP->altTreeAllometryFlag;
	for (int i=190; i<200; i++) cen_state[i] = modelParamsP->code_flags[i - 190] ? 1 : 0;

	logical initflag = TRUE;
	logical onestep = TRUE;

	cen_init_(
			&(runParamsP->actual_years_to_run),
			&initial_vclass,
			&diagsLogical, &onestep,
			cen_state,
			EQdata.cenOutvars,
			modelParamsP->century_path,
			&initial_mix_index,
			&initial_c3c4_ratio,
			&initflag,
			&fire_mode,
			&initial_mapss_climate_zone,
			&initial_tmp_index,
			&initial_ppt_index,
			&unlimitedN,
			&m_frost_index_this_year,
			runParamsP->CO2_file);
	m_fire_last_month = false; // *** this should be restored from the saved state, not just set to false 
	vtype = initial_vtype;
	tree_typ = initial_tree_typ;
	// end of first year initialization and cen_init() call 

	int yrs_offset = runParamsP->years_offset_into_input_data;
	int yrs_needed = runParamsP->actual_years_to_run + (modelParamsP->southernHemisphereFlag ? 1 : 0);
	assert((yrs_offset + yrs_needed)<=runParamsP->years_of_climate_data);
	for (int yrNdx = 0; rtnFlag && yrNdx<runParamsP->actual_years_to_run; yrNdx++) 
		rtnFlag = simulate1pt1yr(&fmi, &biogeogModelInstance, &prev_vclass, &vtype, &tree_typ, yrNdx, yrNdx,
				cen_state, tmp_smoothed, tmax_smoothed, tmin_smoothed, tdmean_smoothed, ppt_smoothed, initial_tmp_index, initial_ppt_index);

	cen_end_(); // Century closes files 

	for (int i = 0; i<NUM_WS_OUTVARS; i++) stateData.eqState.cenOutvars[i] = EQdata.cenOutvars[i];

	for (int mo = 0; mo<12; mo++)
	{
		stateData.eqState.cenOutvars[CENppt_prev_0 + mo] = ppt_smoothed[mo];
		stateData.eqState.cenOutvars[CENtmp_prev_0 + mo] = tmp_smoothed[mo];
		stateData.eqState.cenOutvars[CENtmax_prev_0 + mo] = tmax_smoothed[mo];
		stateData.eqState.cenOutvars[CENtmin_prev_0 + mo] = tmin_smoothed[mo];
		stateData.eqState.cenOutvars[CENtdmean_prev_0 + mo] = tdmean_smoothed[mo];
	}
	stateData.eqState.cenOutvars[CENfire_last_month] = m_fire_last_month ? 1. : 0.; 
	stateData.eqState.cenOutvars[BSvtype] = vtype;
	stateData.eqState.cenOutvars[BSsmoothed_max_tree_lai] = m_smoothed_max_tree_lai;
	stateData.eqState.cenOutvars[BSsmoothed_max_grass_lai] = m_smoothed_max_grass_lai;
	fmi.stateData.fireState.yrs_since_fire = m_burn_count;
	stateData.fireState = fmi.stateData.fireState; // Get the fire state data back into the CENTURY model's
	// state data structure

	if (rtnFlag and runParamsP->multiyrOutputFlag && pS->m_num_multiyr_vars>0) pS->saveMultiyrOutputData();

	return(rtnFlag);
} // end of bool CENTURY_BiogeochemModel::mc2_transient()



///////////////////////////////// mc2_spinup() //////////////////////////////////////////////////////  

bool CENTURY_BiogeochemModel::mc2_spinup()
{
	Soil_MAPSS soils_data(inputDataP->soilData);
	SoilsType	WhichSoils = ScsSoilsData;
	bool rtnFlag = true; 
	VCategory initial_vclass;
	MAPSSvegClass initial_class;
	int initial_mapss_climate_zone;
	int initial_vtype, vtype;
	float initial_mix_index, initial_c3c4_ratio, initial_tmp_index, initial_ppt_index;
	float initial_everg, initial_needl;
	TreeType initial_tree_typ, tree_typ;
	float	cen_state[CEN_STATE_SIZE];
	VCategory prev_vclass;
	float ppt_smoothed[12], tmp_smoothed[12], tmax_smoothed[12], tmin_smoothed[12], tdmean_smoothed[12];
	float tau; // for efolding
	logical diagsLogical = runParamsP->diags ? TRUE : FALSE;  
	int efold_init_years, yr;
	int spinup_yr, climate_yr;
	int save_num_vpr_tmp_issues = 0;

	assert(runParamsP->runMode==SPINUP);

	pS->VarDict[ELEV].save(inputDataP->elev);
	pS->VarDict[SWHC_TOP].save(soils_data.swhc[0]);
	pS->VarDict[SWHC_MID].save(soils_data.swhc[1]);
	pS->VarDict[SWHC_DEEP].save(soils_data.swhc[2]);
	pS->VarDict[SWHC_ALL].save(soils_data.swhc[0] + soils_data.swhc[1] + soils_data.swhc[2]);
	pS->VarDict[SOIL_DEPTH].save(soils_data.mineral_depth);

	/* Note that the EQdata.cenOutVars has been initialized from the netCDF file written out at the end of the CENTURY_EQ phase. */

	/* Initialize climate smoothing for ppt, tmp, tmax, tmin, tdmean. */
	tau = modelParamsP->efold_t_max;
	efold_init_years = MIN(1.5*modelParamsP->efold_t_max, runParamsP->years_of_climate_data);
	assert(runParamsP->years_of_climate_data>=2); // must have at least 2 years of data to do southern hemisphere shift
	for (int mo = 0; mo<12; mo++)
	{ int source_mo;
		source_mo = modelParamsP->southernHemisphereFlag ? (mo + 6) : mo;
		ppt_smoothed[mo] = *(inputDataP->pptP + source_mo);
		tmp_smoothed[mo] = *(inputDataP->tmpP + source_mo);
		tmax_smoothed[mo] = *(inputDataP->tmaxP + source_mo);
		tmin_smoothed[mo] = *(inputDataP->tminP + source_mo);
		tdmean_smoothed[mo] = *(inputDataP->tdmeanP + source_mo);
	}
	yr = 1;
	while (yr<efold_init_years)
	{
		for (int mo = 0; mo<12; mo++)
		{ int source_mo;
			float raw_ppt, raw_tmp, raw_tmax, raw_tmin, raw_tdmean;
			source_mo = yr*12 + (modelParamsP->southernHemisphereFlag ? (mo + 6) : mo);
			raw_ppt = *(inputDataP->pptP + source_mo);
			ppt_smoothed[mo] = sciFn.efold(tau, raw_ppt, ppt_smoothed[mo]);
			raw_tmp = *(inputDataP->tmpP + source_mo);
			tmp_smoothed[mo] = sciFn.efold(tau, raw_tmp, tmp_smoothed[mo]);
			raw_tmax = *(inputDataP->tmaxP + source_mo);
			tmax_smoothed[mo] = sciFn.efold(tau, raw_tmax, tmax_smoothed[mo]);
			raw_tmin = *(inputDataP->tminP + source_mo);
			tmin_smoothed[mo] = sciFn.efold(tau, raw_tmin, tmin_smoothed[mo]);
			raw_tdmean = *(inputDataP->tdmeanP + source_mo);
			tdmean_smoothed[mo] = sciFn.efold(tau, raw_tdmean, tdmean_smoothed[mo]);
		}
		yr++;
	}

	if (modelParamsP->code_flags[MAPSS_IN_SPINUP_CODE_FLAG])
	{ // Run MAPSS now to get various initial conditions.  This is how MC1 does it.
		MAPSS_BiogeogModel mapssInstance(pS);

		mapssInstance.MAPSSmodelInit(modelParamsP->MAPSSparameterSet, inputDataP->soilData,
				inputDataP->pptP, inputDataP->tmpP, inputDataP->vprP, inputDataP->wndP, inputDataP->elev);

		// sciFn.MixIndex(inputDataP->pptP, inputDataP->tmpP, modelParamsP->p_hi_mult, 
		//      &initial_tmp_index, &initial_ppt_index, &initial_mix_index);
		m_frost_index_this_year = sciFn.FrostIndex(inputDataP->tmpP);

		/* Run the MAPSS model. */
		MAPSSdataOutputs mapssOutput;
		initial_class = mapssInstance.mapss_model(&mapssOutput); // Run the MAPSS model on one point.
		// initial_c3c4_ratio = mapssInstance.m_c3pct;
		initial_mapss_climate_zone = mapssInstance.m_zone;      
	}
	else 
	{ // Get initial conditions from EQdata.
		m_frost_index_this_year = EQdata.cenOutvars[CENfrost_index_this_year]; // maybe this should be the frost index of the first year of the spinup climate
		initial_class = (MAPSSvegClass)EQdata.cenOutvars[CENinitial_class_mapss];
		initial_mapss_climate_zone = EQdata.cenOutvars[CENmapss_zone];
	}
	m_fire_last_month = false; // There is no dynamic fire in EQ.

	initial_tmp_index = EQdata.cenOutvars[CENtmp_index];
	initial_ppt_index = EQdata.cenOutvars[CENppt_index];
	initial_mix_index = EQdata.cenOutvars[CENmx_index];
	initial_c3c4_ratio = EQdata.cenOutvars[CENc3c4_index];
	initial_vclass = (VCategory)EQdata.cenOutvars[CENinitial_vclass_mapss];
	initial_vtype = convertVEMAPtoVtype(initial_vclass, runParamsP->baseCalibration);
	sciFn.PhytomorphologicalIndices(initial_tmp_index, initial_ppt_index, &initial_everg, &initial_needl);
	treeType(initial_everg, initial_needl, 
			sciFn.annual_sum(ppt_smoothed, modelParamsP->southernHemisphereFlag),
			runParamsP->baseCalibration, &initial_tree_typ, &stateData.deciduous, &stateData.broadleaf);

	// stateData.initializeSpinupData();
	cen_zero_(); /* Zero out the labeled commons in century. */

	/* Pass soils info to Century. */
	float sand, clay, depth, rock[NSL];
	sand = (float) soils_data.sand[SURFACE];
	clay = (float) soils_data.clay[SURFACE];
	depth = (float) soils_data.mineral_depth;
	if (depth<=10.) return(false); // If depth is <= 10 mm, quit.
	for (int i = SURFACE; i <= DEEP; i++) rock[i] = (float) soils_data.rock_frag_mineral[i];  
	float hdeposition = 0.;
	float cdeposition = 0.;
	cen_init_soils_(
			&soils_data.bulk_density,
			&sand, 
			&clay,
			&depth,
			rock,
			&WhichSoils,
			&hdeposition,
			&cdeposition);

	/* Pass latitude and elevation to Century. */
	float latitude = (float)inputDataP->lat;
	cen_init_lat_(&latitude, &inputDataP->elev);

	/* Zero out annet for cen_init_() in order to get wdfxs right in eachyr.F */
	EQdata.cenOutvars[549] = 0.0; 

	int fire_mode = runParamsP->fireModelFlag ? 1 : 0;
	int unlimitedN = modelParamsP->unlimitedNflag ? 1 : 0;

	MC_BiogeographyModel biogeogModelInstance(pS, runParamsP, modelParamsP);
	MC_FireModel fmi(pS, &stateData.fireState, inputDataP);
	if (runParamsP->fireModelFlag) fmi.portable_srand(1);

	assert(runParamsP->years_offset_into_input_data==0);
	m_smoothed_max_tree_lai = 0.0;
	m_smoothed_max_grass_lai = 0.0;
	m_burn_count = 0; // Initialize number of years since last simulated fire; maybe this should be a large number, for faster spinups

	// first year initialization and cen_init() call
	prev_vclass = initial_vclass;
	if (initial_class<=NaturalBarren || initial_class==Ice || initial_class==DesertExtreme || initial_vclass==VUnknown) 
		return(false);  

	if (runParamsP->diags) 
	{
		printf("class, vclass: %d, %d \n", initial_class, initial_vclass);
		printf("tmp_index: %.2f \n", initial_tmp_index);
		printf("ppt_index: %.2f \n", initial_ppt_index);
		printf("c3c4_ratio: %.2f \n", initial_c3c4_ratio);
	}

	for (int i=0; i<187; i++) cen_state[i] = 0.0;
	cen_state[187] = modelParamsP->m_century_runoff_x_intercept;
	cen_state[188] = modelParamsP->m_century_runoff_slope;
	cen_state[189] = modelParamsP->m_lait_lower_limit;
	for (int i=190; i<200; i++) cen_state[i] = modelParamsP->code_flags[i - 190] ? 1 : 0;
	logical initflag = TRUE;
	logical onestep = TRUE;

	cen_init_(
			&(runParamsP->actual_years_to_run),
			&initial_vclass,
			&diagsLogical, &onestep,
			cen_state,
			EQdata.cenOutvars,
			modelParamsP->century_path,
			&initial_mix_index,
			&initial_c3c4_ratio,
			&initflag,
			&fire_mode,
			&initial_mapss_climate_zone,
			&initial_tmp_index,
			&initial_ppt_index,
			&unlimitedN,
			&m_frost_index_this_year,
			runParamsP->CO2_file);
	m_fire_last_month = false;
	m_c3pct = initial_c3c4_ratio;
	vtype = initial_vtype;
	tree_typ = initial_tree_typ;
	// end of first year initialization and cen_init() call

	int eff_spinup_cycle_length = runParamsP->years_of_climate_data;
	if (modelParamsP->southernHemisphereFlag) eff_spinup_cycle_length--; 
	assert(eff_spinup_cycle_length>0);

	for (spinup_yr = 0; spinup_yr<runParamsP->actual_years_to_run; spinup_yr += eff_spinup_cycle_length) 
	{
		for (climate_yr = 0; rtnFlag && climate_yr<eff_spinup_cycle_length; climate_yr++) 
			rtnFlag = simulate1pt1yr(&fmi, &biogeogModelInstance, &prev_vclass, &vtype, &tree_typ, 
					(spinup_yr + climate_yr), climate_yr, cen_state,
					tmp_smoothed, tmax_smoothed, tmin_smoothed, tdmean_smoothed, ppt_smoothed, initial_tmp_index, initial_ppt_index);

		// Count vpr_tmp_issues only on the 1st pass thru the climate data.
		if (spinup_yr==0) save_num_vpr_tmp_issues = pS->m_num_vpr_tmp_issues; 

	} // end of spinup_yr loop

	pS->m_num_vpr_tmp_issues = save_num_vpr_tmp_issues;

	cen_end_(); /* Century closes files */

	for (int i = 0; i<NUM_WS_OUTVARS; i++) stateData.eqState.cenOutvars[i] = EQdata.cenOutvars[i];

	for (int mo = 0; mo<12; mo++)
	{
		stateData.eqState.cenOutvars[CENppt_prev_0 + mo] = ppt_smoothed[mo];
		stateData.eqState.cenOutvars[CENtmp_prev_0 + mo] = tmp_smoothed[mo];
		stateData.eqState.cenOutvars[CENtmax_prev_0 + mo] = tmax_smoothed[mo];
		stateData.eqState.cenOutvars[CENtmin_prev_0 + mo] = tmin_smoothed[mo];
		stateData.eqState.cenOutvars[CENtdmean_prev_0 + mo] = tdmean_smoothed[mo];
	}
	stateData.eqState.cenOutvars[CENclass_mapss] = initial_class;
	stateData.eqState.cenOutvars[CENvclass_mapss] = -9999;
	stateData.eqState.cenOutvars[CENfire_last_month] = m_fire_last_month ? 1. : 0.; 
	stateData.eqState.cenOutvars[BSvtype] = vtype;
	stateData.eqState.cenOutvars[BSsmoothed_max_tree_lai] = m_smoothed_max_tree_lai;
	stateData.eqState.cenOutvars[BSsmoothed_max_grass_lai] = m_smoothed_max_grass_lai;
	fmi.stateData.fireState.yrs_since_fire = m_burn_count;
	stateData.fireState = fmi.stateData.fireState; // Get the fire state data back into the CENTURY model's
	// state data structure

	if (rtnFlag and runParamsP->multiyrOutputFlag && pS->m_num_multiyr_vars>0) pS->saveMultiyrOutputData();

	return(rtnFlag);
} // end of mc2_spinup()


bool CENTURY_BiogeochemModel::century_eq_model()
{
	VCategory vclass; // VEMAP vegetation class

	Soil_MAPSS soils_data(inputDataP->soilData);
	SoilsType	WhichSoils = ScsSoilsData;
	int			years;
	float			cen_state[CEN_STATE_SIZE];
	int		mo, i;
	float		sand, clay, depth, rock[3];
	float		final_growth = 0.67;
	float ppt[12], tmin[12], tmax[12], tmp[12];
	MAPSSvegClass mclass; 
	int zone;

	mclass = (MAPSSvegClass)mapssOutput.mclass;
	zone = mapssOutput.zone;
	vclass = convertMAPSStoVEMAP(mclass);
	if (vclass==VUnknown) 
	{ 
		printf("century_eq_model(): vclass==VUnknown. mclass, vclass = %d, %d\n", mclass, vclass);
		return(false); 
	}

	for (i=0; i<187; i++) cen_state[i] = 0.0;
	cen_state[187] = modelParamsP->m_century_runoff_x_intercept;
	cen_state[188] = modelParamsP->m_century_runoff_slope;
	cen_state[189] = modelParamsP->m_lait_lower_limit;
	for (i=190; i<200; i++) cen_state[i] = modelParamsP->code_flags[i - 190];

	for (mo = 0; mo<12; mo++)
	{ int source_mo;
		source_mo = modelParamsP->southernHemisphereFlag ? (mo + 6) % 12 : mo;
		ppt[mo] = *(inputDataP->pptP + source_mo);
		tmin[mo] = *(inputDataP->tminP + source_mo);
		tmax[mo] = *(inputDataP->tmaxP + source_mo);
		tmp[mo] = *(inputDataP->tmpP + source_mo);
	}
	cen_init_climate_(ppt, tmax, tmin); 

	/* Pass soils info and latitude to Century. */
	sand = (float) soils_data.sand[SURFACE];
	clay = (float) soils_data.clay[SURFACE];
	depth = (float) soils_data.mineral_depth;
	if (depth<=10.) 
	{ // When depth is <= 10 mm, quit.
		printf("century_eq_model(): soil_depth = %f, less than minimum for simulation of vegetation\n", depth);
		return(false); 
	}
	for ( i = SURFACE; i <= DEEP; i++) rock[i] = (float) soils_data.rock_frag_mineral[i];

	float hdeposition = 0.;
	float cdeposition = 0.;
	cen_init_soils_(
			&soils_data.bulk_density,
			&sand, 
			&clay,
			&depth,
			rock,
			&WhichSoils,
			&hdeposition,
			&cdeposition);

	float latitude = (float)inputDataP->lat;
	cen_init_lat_(&latitude, &inputDataP->elev);

	float tmp_index, ppt_index, mix_index;
	sciFn.MixIndex(ppt, tmp, modelParamsP->p_hi_mult, &tmp_index, &ppt_index, &mix_index);
	m_frost_index_this_year = sciFn.FrostIndex(tmp);
	m_c3pct = sciFn.C3prodPct(tmp);
	int fire_mode = runParamsP->fireModelFlag ? 1 : 0;
	int unlimitedN = modelParamsP->unlimitedNflag ? 1 : 0;
	int diags = runParamsP->diags ? 1 : 0;
	logical onestep = FALSE;
	logical initflag = FALSE;
	cen_init_(
			&(runParamsP->actual_years_to_run),
			&vclass,
			&diags, &onestep,
			cen_state,
			EQdata.cenOutvars, // entire array is initialized to 0.0, then some century variables are copied into the initial portion
			modelParamsP->century_path,
			&mix_index,
			&m_c3pct,
			&initflag,
			&fire_mode,
			&zone,
			&tmp_index,
			&ppt_index,
			&unlimitedN,
			&m_frost_index_this_year,
			runParamsP->CO2_file);        

	/** Step through standard Century to get output and state, until EQ **/
	years = 0;
	logical eq_flag;
	do 
	{ 		    
		for (mo = JAN; mo <= DEC; mo++)
		{ 
			if (!modelParamsP->code_flags[TREE_N_FLAG])
			{
				cen_state[187] = modelParamsP->m_century_runoff_x_intercept;
				cen_state[188] = modelParamsP->m_century_runoff_slope;
				cen_state[189] = modelParamsP->m_lait_lower_limit;
			}
			for (i=190; i<200; i++) cen_state[i] = modelParamsP->code_flags[i - 190];
			stand_step_(cen_state, EQdata.cenOutvars, &eq_flag, &final_growth, &diags, &onestep, &m_frost_index_this_year); 
		}
		years++;
	} while (eq_flag==FALSE && years<runParamsP->actual_years_to_run);

	printf("Exit equilibration loop w/ eq_flag = %d and years = %d\n", eq_flag, years);

	cen_end_(); /* Century closes files */

	EQdata.cenOutvars[CENfrost_index_this_year] = m_frost_index_this_year;      
	EQdata.cenOutvars[CENclass_mapss] = mclass; // "class_mapss"
	EQdata.cenOutvars[CENinitial_class_mapss] = mclass; // "initial_class_mapss"
	EQdata.cenOutvars[CENvclass_mapss] = vclass; // "vclass_mapss"
	EQdata.cenOutvars[CENinitial_vclass_mapss] = vclass; // "initial_vclass_mapss"
	EQdata.cenOutvars[CENefold_t] = modelParamsP->efold_t_max; 
	EQdata.cenOutvars[CENequil_time] = years - 1; // for compatibility with MC1
	EQdata.cenOutvars[CENfinal_year] = years; // "final_year"
	EQdata.cenOutvars[CENmapss_mix_index] = mix_index;
	EQdata.cenOutvars[CENmapss_tmp_index] = tmp_index;
	EQdata.cenOutvars[CENmapss_ppt_index] = ppt_index;
	EQdata.cenOutvars[CENfire_last_month] = 0.; // There is no dynamic fire in EQ.

	return(true);

} // end of century_eq_model()


bool CENTURY_BiogeochemModel::simulate1pt1yr(
		MC_FireModel * pF,
		MC_BiogeographyModel * pB,
		VCategory * prev_vclassP,
		int * vtypeP,
		TreeType * tree_typP,
		int yrNdx, 
		int climateYrNdx,
		float	cen_state[CEN_STATE_SIZE],
		float tmp_smoothed[12],
		float tmax_smoothed[12],
		float tmin_smoothed[12],
		float tdmean_smoothed[12],
		float ppt_smoothed[12],
		float initial_tmp_index,
		float initial_ppt_index)
{
	bool rtnFlag = true;
	float ppt_raw[12], tmax_raw[12], tmin_raw[12], tmp_raw[12], tdmean_raw[12], vpr_raw[12], rh_raw[12], pet[12];
	float tmp_index, ppt_index, mix_index;

	for (int mo = 0; mo<12; mo++)
	{ int source_mo;
		source_mo = climateYrNdx*12 + (modelParamsP->southernHemisphereFlag ? (mo + 6) : mo);
		assert(source_mo<(12*runParamsP->years_of_climate_data));
		ppt_raw[mo] = *(inputDataP->pptP + source_mo);
		tmin_raw[mo] = *(inputDataP->tminP + source_mo);
		tmax_raw[mo] = *(inputDataP->tmaxP + source_mo);
		tmp_raw[mo] = *(inputDataP->tmpP + source_mo);
		tdmean_raw[mo] = *(inputDataP->tdmeanP + source_mo);
		if (tmax_raw[mo]>100. || tmp_raw[mo]>tmax_raw[mo] || tmin_raw[mo]>tmp_raw[mo])
		{
			printf("*** CENTURY.cpp/simulate1pt1yr(): "
					"climateYrNdx, mo, tmax_raw, tmp_raw, tmin_raw = %d, %d, %f, %f, %f\n"
					"(can be caused by use of non-standard 'scaled' attribute instead of 'scale_factor' attribute in .nc input file\n"
					"or by tmin>tmax or tmp>tmax or tmp<tmin)\n",
					climateYrNdx, mo, tmax_raw[mo], tmp_raw[mo], tmin_raw[mo]);
			if (tmax_raw[mo]>100.) err_exit("anomalous temperature data: tmax>100. degC");
			float abs_diff = fabs(tmax_raw[mo] - tmin_raw[mo]);
			float tmp_calculated = (tmax_raw[mo] + tmin_raw[mo])/2.f;
			tmax_raw[mo] = tmp_calculated + abs_diff/2.f;
			tmin_raw[mo] = tmp_calculated - abs_diff/2.f;
			tmp_raw[mo] = tmp_calculated;
			printf("Attempted to repair temperature data by setting tmax_raw = %f, tmp_raw = %f, tmin_raw = %f\n",
					tmax_raw[mo], tmp_raw[mo], tmin_raw[mo]);
		}
		vpr_raw[mo] = *(inputDataP->vprP + source_mo);

		assert((tmp_raw[mo] + FREEZE)>0.);
		rh_raw[mo] = (vpr_raw[mo]/sciFn.satvp(tmp_raw[mo]))*100.;
		if (rh_raw[mo]>100.) rh_raw[mo] = 100.;
	}

	float annual_ppt = 0.; for (int mo = 0; mo<12; mo++) annual_ppt += ppt_raw[mo];
	pS->VarDict[PPT].save(annual_ppt, YEAR_INTERVAL);
	pS->VarDict[TMP].save(sciFn.annual_average(tmp_raw, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[TMAX].save(sciFn.annual_average(tmax_raw, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[TMIN].save(sciFn.annual_average(tmin_raw, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[TDMEAN].save(sciFn.annual_average(tdmean_raw, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[VPR].save(sciFn.annual_average(vpr_raw, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);

	cen_init_climate_(ppt_raw, tmax_raw, tmin_raw); 
	/* Calculate this year's smoothed ppt, tmp, tmax, tmin, and tdmean values. */
	for (int mo = JAN; mo <= DEC; mo++) 
	{
		ppt_smoothed[mo] = sciFn.efold(modelParamsP->efold_t_max, ppt_raw[mo], ppt_smoothed[mo]);
		tmp_smoothed[mo] = sciFn.efold(modelParamsP->efold_t_max, tmp_raw[mo], tmp_smoothed[mo]);
		tmax_smoothed[mo] = sciFn.efold(modelParamsP->efold_t_max, tmax_raw[mo], tmax_smoothed[mo]);
		tmin_smoothed[mo] = sciFn.efold(modelParamsP->efold_t_max, tmin_raw[mo], tmin_smoothed[mo]);
		tdmean_smoothed[mo] = sciFn.efold(modelParamsP->efold_t_max, tdmean_raw[mo], tdmean_smoothed[mo]);
	}
	pS->VarDict[TMAX_SMOOTHED].save(sciFn.annual_average(tmax_smoothed, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[PPT_SMOOTHED].save(sciFn.annual_sum(ppt_smoothed, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);

	/* Use tmp_smoothed and ppt_smoothed to calculate tmp_index, ppt_index, mix_index, and c3c4_index, 
	   which will be passed in to cen_step().  */ 
	sciFn.MixIndex(ppt_smoothed, tmp_smoothed, modelParamsP->p_hi_mult, &tmp_index, &ppt_index, &mix_index);
	if (runParamsP->runMode==SPINUP)
	{
		if (modelParamsP->code_flags[MAPSS_IN_SPINUP_CODE_FLAG] && yrNdx==0) 
		{
			tmp_index = initial_tmp_index;
			ppt_index = initial_ppt_index;
			// m_c3pct = initial_c3c4_ratio;
		}
		else m_c3pct = sciFn.C3prodPct(tmp_smoothed);
	}
	else m_c3pct = sciFn.C3prodPct(tmp_smoothed);
	m_min_smoothed_tmp = sciFn.annual_min(tmp_smoothed);
	m_max_smoothed_tmp = sciFn.annual_max(tmp_smoothed);
	m_cont_index = m_max_smoothed_tmp - m_min_smoothed_tmp;
	m_gdd_zero = sciFn.GrowingDegreeDays(tmp_smoothed, 0., modelParamsP->southernHemisphereFlag);
	m_zone4biogeog = sciFn.ClimateZone4Biogeography(m_gdd_zero, m_min_smoothed_tmp, c_az_thres, modelParamsP->bz_thres, c_tz_thres, c_stz_thres);
	m_tmmin = sciFn.annual_min(tmp_raw);
	m_frost_index_this_year = sciFn.FrostIndex(tmp_smoothed);

	// Calculate this year's PET for the fire model, from raw tmp and raw vpr, but using the roughness length
	// for a MAPSS climate zone based on smoothed min_tmp.  
	int mapss_climate_zone = sciFn.ClimateZone4MAPSS(m_min_smoothed_tmp, c_n_decid_bound, c_s_decid_bound, c_frost);
	for (int mo = 0; mo<12; mo++) 
	{ float vp_sat;
		assert((tmp_raw[mo] + FREEZE)>0.);
		vp_sat = sciFn.satvp(tmp_raw[mo]);
		pet[mo] = 
			calcPET(tmp_raw[mo], vpr_raw[mo], vp_sat, m_wnd[mo], 
					inputDataP->elev, c_z_all, modelParamsP->z0[mapss_climate_zone][TREE], 
					sciFn.days_per_mo[modelParamsP->southernHemisphereFlag ? (mo + 6)%12 : mo]);
	}

	pF->FireData(ppt_raw, tmp_raw, m_wnd, vpr_raw, pet, tmin_raw, tmax_raw);

	// Simulate one point for one year.
	rtnFlag = simulate1pt1yrInner(
			pF,
			pB,
			prev_vclassP,
			vtypeP,
			tree_typP,
			yrNdx,
			cen_state,
			mapss_climate_zone,
			mix_index,
			tmp_index,
			ppt_index,
			ppt_smoothed,
			tmp_smoothed,
			tmax_smoothed,
			ppt_raw,
			tmp_raw,
			tmax_raw,
			tmin_raw,
			tdmean_raw,
			vpr_raw,
			rh_raw);

	if (!rtnFlag) return(false);

	assert(!runParamsP->spaceBeforeTimeFlag);
	// Save away this year's output data.
	// The data for all the years for this cell will be written to the output file after the
	// simulation is completed for this cell.
	pS->VarDict[NPP].save(m_npp_yr, YEAR_INTERVAL);
	pS->VarDict[VTYPE].save(*vtypeP, YEAR_INTERVAL);
	pS->VarDict[BIOME].save(m_biome, YEAR_INTERVAL);
	pS->VarDict[PHYSIOGNOMIC_CLASS].save(m_physiognomic_class, YEAR_INTERVAL);
	pS->VarDict[C_VEG].save(sciFn.annual_average(m_vegc_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_FOREST].save(sciFn.annual_average(m_frstc_mo, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C_MAX_LIVE_GRASS_ABOVEGR].save(sciFn.annual_max(m_aglivc_mo), YEAR_INTERVAL);
	pS->VarDict[GFRAC].save(m_max_gfrac, YEAR_INTERVAL);
	pS->VarDict[NEP].save(m_nep_yr, YEAR_INTERVAL);
	pS->VarDict[NBP].save(m_nbp_yr, YEAR_INTERVAL);
	pS->VarDict[RSP].save(m_rsp_yr, YEAR_INTERVAL);
	pS->VarDict[BIO_CONSUME_CENTURY].save(m_bio_consume_century, YEAR_INTERVAL);
	pS->VarDict[RH].save(sciFn.annual_average(rh_raw, modelParamsP->southernHemisphereFlag), YEAR_INTERVAL);
	pS->VarDict[C3_PCT_PROD].save(m_c3pct, YEAR_INTERVAL);
	pS->VarDict[TREE_TYPE].save(*tree_typP, YEAR_INTERVAL);
	pS->VarDict[MIN_SMOOTHED_TMP].save(m_min_smoothed_tmp, YEAR_INTERVAL);
	pS->VarDict[MC_CLIMATE_ZONE].save(m_zone4biogeog, YEAR_INTERVAL);
	pS->VarDict[SNOWPACK_DAY91].save(pF->m_snow_doy[90], YEAR_INTERVAL);
	pS->VarDict[CONTINENTALITY].save(m_cont_index, YEAR_INTERVAL);
	pS->VarDict[GDD].save(m_gdd_zero, YEAR_INTERVAL);
	pS->VarDict[SOIL_TMP_MAX].save(sciFn.estimate_max_soil_tmp(tmp_smoothed), YEAR_INTERVAL);


	if (pS->m_num_yearly_vars>0) pS->saveOneYearOneCellOutputData(yrNdx);

	return(rtnFlag);
} // end of bool CENTURY_BiogeochemModel::simulate1pt1yr()




