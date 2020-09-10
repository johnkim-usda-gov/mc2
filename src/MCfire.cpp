/*
 *  MCfire.cpp
 *  mc2
 *
 *  Created by Dave Conklin on 3/14/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>	// JK
#include "assert.h"
#include "math.h"

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


MC_FireModel::MC_FireModel(Simulation * pRun, FireState * fireStateP, InputDataClass * inputP)
	// Use this constructor when instantiating a fire model object with no previous state.
{
	inputDataP = inputP;
	ProcessModelInit(pRun);
	initializeFireModelConstantsAndOutputs();

	// Initialize fire model state variables
	stateData.fireState.sum_ann_ppt = 0.;
	stateData.fireState.Pprev = 0.f;
	stateData.fireState.Kprev = 0.f;
	stateData.fireState.Ksum = 0.f;
	stateData.fireState.prev_l1hr = 0.;
	stateData.fireState.prev_dstand = 0.;
	stateData.fireState.prev_litter = 0.;
	stateData.fireState.prev_1hr = 0.;
	stateData.fireState.prev_10hr = 0.;
	stateData.fireState.prev_100hr = 0.;
	stateData.fireState.prev_1000hr = 0.;
	stateData.fireState.prev_day_mc_100hr = 35.; 
	stateData.fireState.prev_day_mc_1000hr = 35.; 
	stateData.fireState.prev_mxd = 25.;
	stateData.fireState.prev_mc_grass = 30.;
	stateData.fireState.prev_mc_tree = 80.;
	stateData.fireState.prev_depth = 0.;
	stateData.fireState.clai_in_midseason = 0.f;
	stateData.fireState.prev_snow = 0.;
	stateData.fireState.F1 = 0;
	stateData.fireState.F2 = 0;
	stateData.fireState.rand_seed = 1; 
	stateData.fireState.yrs_in_sum_ann_ppt = 0;
	stateData.fireState.ffmc_prev = 85.0f;
	stateData.fireState.dmc_prev = 6.0f;
	stateData.fireState.dc_prev = 15.0f;
	stateData.fireState.yrs_since_fire = 0; // maybe this should be a large number

	*fireStateP = stateData.fireState; // Copy fire state variables to the state variables of the caller.

} // end of constructor MC_FireModel(Simulation * pRun, FireState * fireStateP, InputDataClass * inputP)


MC_FireModel::MC_FireModel(Simulation * pRun, float cenOutvars[], InputDataClass * inputP)
	// Use this constructor when instantiating a fire model object and restoring a previous state.
{
	inputDataP = inputP;
	ProcessModelInit(pRun);
	initializeFireModelConstantsAndOutputs();

	// Restore fire model state variables.
	stateData.fireState.sum_ann_ppt = cenOutvars[FSsum_ann_ppt];
	stateData.fireState.Pprev = cenOutvars[FS_Pprev];
	stateData.fireState.Kprev = cenOutvars[FS_Kprev];
	stateData.fireState.Ksum = cenOutvars[FS_Ksum];
	stateData.fireState.prev_l1hr = cenOutvars[FSprev_l1hr];
	stateData.fireState.prev_dstand = cenOutvars[FSprev_dstand];
	stateData.fireState.prev_litter = cenOutvars[FSprev_litter];
	stateData.fireState.prev_1hr = cenOutvars[FSprev_1hr];
	stateData.fireState.prev_10hr = cenOutvars[FSprev_10hr];
	stateData.fireState.prev_100hr = cenOutvars[FSprev_100hr];
	stateData.fireState.prev_1000hr = cenOutvars[FSprev_1000hr];
	stateData.fireState.prev_day_mc_100hr = cenOutvars[FSprev_day_mc_100hr]; 
	stateData.fireState.prev_day_mc_1000hr = cenOutvars[FSprev_day_mc_1000hr]; 
	stateData.fireState.prev_mxd = cenOutvars[FSprev_mxd];
	stateData.fireState.prev_mc_grass = cenOutvars[FSprev_mc_grass];
	stateData.fireState.prev_mc_tree = cenOutvars[FSprev_mc_tree];
	stateData.fireState.prev_depth = cenOutvars[FSprev_depth];
	stateData.fireState.clai_in_midseason = cenOutvars[FSclai_in_midseason];
	stateData.fireState.prev_snow = cenOutvars[FSprev_snow];
	stateData.fireState.F1 = cenOutvars[FS_F1];
	stateData.fireState.F2 = cenOutvars[FS_F2];
	stateData.fireState.rand_seed = cenOutvars[FSrand_seed]; 
	stateData.fireState.yrs_in_sum_ann_ppt = cenOutvars[FSyrs_in_sum_ann_ppt];
	stateData.fireState.ffmc_prev = cenOutvars[FSffmc_prev];
	stateData.fireState.dmc_prev = cenOutvars[FSdmc_prev];
	stateData.fireState.dc_prev = cenOutvars[FSdc_prev];
	stateData.fireState.yrs_since_fire = cenOutvars[FSyrs_since_fire];

} // end of constructor MC_FireModel(pRun, cenOutvars[], inputP)


void MC_FireModel::initializeFireModelConstantsAndOutputs()
{
	c_stl = c_std = 0.055; 
	// mineral damping coefficent of live and dead fuels    
	float sd, sl;
	sd = sl = 0.01;
	c_etasd = 0.174 * pow(sd, -0.19);
	c_etasl = 0.174 * pow(sl, -0.19);      

	fireParams.snw0 = 3.0f;
	fireParams.snw1 = 0.0f;
	fireParams.no_melt = 1.0f;
	fireParams.melt_b = 1.5f;

	fireParams.mc_tree_min = 80.f;
	fireParams.mc_tree_max = 130.f;
	fireParams.mc_grass_min = 30.f;
	fireParams.mc_grass_max = 120.f;

	fireParams.slp = 0.f;
	fireParams.prob_thres = 45.f;

	// These are output variables.
	m_d1hr = m_d10hr = m_d100hr = m_d1000hr = NC_FILL_FLOAT;
} // end of initializeFireModelConstantsAndOutputs()


bool MC_FireModel::runModelAndCollectData(const int year_index, const int row_index, const int col_index)
{
	return(false);
} // end of MC_FireModel::runModelAndCollectData()


void MC_FireModel::FireData(float * pptP, float * tmpP, float * wndP, float * vprP, float * petP,
		float * tminP, float * tmaxP)
// FireData: Called at the beginning of each simulation year.  Estimates daily values of climate variables. 
{
	float snowpack_carryover;
	float ann_ppt;
	float satvp_tmax, satvp_tmin;

	m_fire_ros = 0.;
	m_last_yr_snow = stateData.fireState.prev_snow;

	/* Get monthly climate data */ 
	for (int i = 0; i<12; i++)
	{      
		m_tmp[i]  = *(tmpP + i);
		m_ppt[i] = *(pptP + i);
		m_wnd[i] = *(wndP + i);
		m_rh[i]   = (*(vprP + i))/sciFn.satvp(m_tmp[i]) * 100.;
		if (m_rh[i]>100.) 
		{
			m_rh[i] = 100.;
			pS->m_num_vpr_tmp_issues++;
			// printf("*** FireData(): mo, vpr[mo], tmp[mo], satvp(tmp[mo]) = %d, %f, %f, %f\n"
			//      "vpr is inconsistent with tmp\n", i, *(vprP + i), m_tmp[i], sciFn.satvp(m_tmp[i]));
		}
		m_pet[i] = *(petP + i);
		m_tmin[i] = *(tminP + i);
		m_tmax[i] = *(tmaxP + i);

		satvp_tmin = sciFn.satvp((double) m_tmin[i]); 
		m_rhmax[i] = ((*(vprP + i))/satvp_tmin) * 100.;
		if (m_rhmax[i] > 100.) m_rhmax[i] = 100.;

		satvp_tmax = sciFn.satvp((double) m_tmax[i]);
		m_rhmin[i] = ((*(vprP + i))/satvp_tmax) * 100.; 
		if (m_rhmin[i] > 100.) m_rhmin[i] = 100.;

		// estimate monthly rainfall intensity using frontal vs. convectional pet threshold    
		m_ppt_rat_in_per_hr[i] = (m_pet[i] > 100.) ? 0.25 : 0.05; // inches per hour 

	} // end of monthly loop

	// Estimate daily values
	{
		float tmin_doy[365], tmax_doy[365];

		daily_ppt();    
		daily_dat(m_tmp, m_tmp_doy);
		daily_dat(m_tmax, tmax_doy);
		daily_dat(m_tmin, tmin_doy);
		for (int doy = 0; doy<365; doy++)
		{
			m_tmp_doy_F[doy] = sciFn.CtoF(m_tmp_doy[doy]);
			m_tmin_doy_F[doy] = sciFn.CtoF(tmin_doy[doy]);
			m_tmax_doy_F[doy] = sciFn.CtoF(tmax_doy[doy]);
		}
		daily_dat(m_wnd, m_wnd_doy);
		daily_dat(m_ppt_rat_in_per_hr, m_ppt_rat_doy_in_per_hr);
		daily_dat(m_rh, m_rh_doy);
		daily_dat(m_rhmin, m_rhmin_doy);
		daily_dat(m_rhmax, m_rhmax_doy);

		/* estimate daily pet */
		/*      int doy = 0;
			for (int mo = 0; mo < 12; mo++)
			{ int month_len;
			month_len = sciFn.days_per_mo[modelParamsP->southernHemisphereFlag? (mo + 6)%12 : mo];
			for(int dom = 0; dom<month_len; dom++)  
			{
			m_pet_doy[doy] = m_pet[mo]/month_len;
			doy++;
			}
			}
			*/
	} // end of block to estimate daily values

	snow_cond(); // run snowpack and snowmelt calculations 

	snowpack_carryover = m_snow_doy[364];
	// prev_snow shouldn't be allowed to exceed one year's precip or 1000 mmH2O, whichever is greater
	if (snowpack_carryover>m_tot_snowfall && snowpack_carryover>1000.) snowpack_carryover = 
		m_tot_snowfall>1000. ? m_tot_snowfall : 1000.;
	stateData.fireState.prev_snow = snowpack_carryover;
	// One year's precip can't exceed 12 months x 3276.7 mmH2O/month = 39320.4 mmH2O.
	assert(stateData.fireState.prev_snow>=0. && stateData.fireState.prev_snow<=39320.4);

	/* calc running ave for ann ppt for kbdi func */     
	ann_ppt = 0.;      
	for (int i = 0; i<12; i++) ann_ppt += m_ppt[i];

	stateData.fireState.sum_ann_ppt += ann_ppt; stateData.fireState.yrs_in_sum_ann_ppt++;
	m_ave_ann_ppt = stateData.fireState.sum_ann_ppt/((float)stateData.fireState.yrs_in_sum_ann_ppt);

	kbdi();

} // end of MC_FireModel::FireData()


void MC_FireModel::daily_ppt()
{
	float ppt_per_event, ppt_per_event_in, rnd_hi;   
	int i, j, k, events, day[31], no_random, dup_flag;
	int juldat = 0;
	float ppt_events[12];

	for (int day_index = 0; day_index<365; day_index++) 
		m_ppt_doy[day_index] = m_ppt_doy_in[day_index] = 0.;

	for (i = 0; i<12; i++) // for each month of the year
	{ 
		int month_length = sciFn.days_per_mo[modelParamsP->southernHemisphereFlag ? ((i + 6)%12) : i];

		// estimate no. ppt events using empirical functions      
		ppt_events[i] = (m_pet[i]>100.) ? 
			0.887 + (0.052 * m_ppt[i]) : 2.319 + (0.027 * m_ppt[i]);
		if (ppt_events[i]>month_length) ppt_events[i] = month_length;
		rnd_hi = ceil(ppt_events[i]); // round up ppt_events             

		ppt_per_event = m_ppt[i]/rnd_hi; // calculate ppt per ppt event 
		assert(ppt_per_event>=0.);
		ppt_per_event_in = sciFn.mm_to_in(ppt_per_event);

		events = (int)rnd_hi;

		j = 0;
		while (j<events) // for the number of days with ppt
		{             
			while (1) // get a random julian day within current month 
			{
				portable_srand(stateData.fireState.rand_seed);
				stateData.fireState.rand_seed = no_random = portable_rand();
				day[j] = (no_random % month_length) + juldat; 
				if (j==0) break;               
				dup_flag = 0; 
				for (k = 0; k<j; k++) if(day[j]==day[k]) dup_flag = 1;

				if(dup_flag == 0) break;
			}

			m_ppt_doy[day[j]] = ppt_per_event; // assign ppt_per_event to randomly chosen day
			m_ppt_doy_in[day[j]] = ppt_per_event_in; 
			j++;
		} // end of loop on days with ppt

		juldat += month_length; // calculate julian date at begining of next month  

	} // end of loop on months of the year

} // end of MC_FireModel::daily_ppt()


int MC_FireModel::daily_values(float monthly_vals[], int mo, float t0, float no_days, float daily_vals[])
	// Interpolates values for noon on each day from the midpoint of one month thru the day before the midpoint of the next month
	// Returns last doy for which a value was filled in.
{
	int initial_doy = (int)t0;
	float t = initial_doy + 0.5; // We estimate values for noon on doy.
	int next_mo = (mo + 1)%12;
	float slope = (monthly_vals[next_mo] - monthly_vals[mo])/no_days;
	float intercept = monthly_vals[mo] - t0*slope;

	int doy = initial_doy;
	int final_doy = -1;
	while (t<(t0 + no_days))
	{
		daily_vals[doy] = slope*t + intercept;
		final_doy = doy;
		t += 1.;
		doy = (doy + 1)%365;
	}

	// Do not allow roundoff to generate values outside the range of the two monthly values.
	// This is necessary to keep relative humidity from slightly exceeding 100%.
	// That happened at row 0, col 52 of the VEMAP spinup dataset in year 183, month 6.  DRC 7/31/12
	daily_vals[initial_doy] = sciFn.clip(daily_vals[initial_doy], monthly_vals[mo], monthly_vals[next_mo]);
	daily_vals[final_doy] = sciFn.clip(daily_vals[final_doy], monthly_vals[mo], monthly_vals[next_mo]);

	return(final_doy);
} // end of MC_FireModel::daily_values()


void MC_FireModel::daily_dat(float monthly_vals[], float daily_vals[])
	// daily_dat() uses 12 monthly values to estimate 365 daily values
	// The monthly values are interpreted as the value at the midpoint of the month.
	// For months with 30 days, the midpoint is midnight between the 15th day and the
	// 16th day.  For months with 31 days, the midpoint is noon on the 16th day.
	// Values at noon on the days between the midpoints of the adjacent months are 
	// determined by linear interpolation between the values at the midpoints.
	// The value for the 12th month is also used for the month prior to the start of 
	// the 12 month series. The 12 months may start in January or July, depending on
	// whether the study area is in the northern or the southern hemisphere.
{   
	float t0;

	t0 = 0.5*sciFn.days_per_mo[modelParamsP->southernHemisphereFlag ? 6 : 0]; // will be 15.5
	for (int mo = 0; mo < 12; mo++) 
	{
		int day, last_day;
		float no_days;

		int first_month = modelParamsP->southernHemisphereFlag ? ((mo + 6)%12) : mo;
		int second_month = (mo + 1)%12;
		no_days = (sciFn.days_per_mo[first_month] + sciFn.days_per_mo[second_month])/2.0;

		last_day = daily_values(monthly_vals, mo, t0, no_days, daily_vals);
		t0 += no_days;
		day = (int)t0;
		assert( ((last_day + 1)==day) || ((mo==11) && (last_day + 366)==day) );

	}
} // end of MC_FireModel::daily_dat()


void MC_FireModel::snow_cond()
{   
	float fall_a, fall_b;
	float melt_a;
	float pot_melt, snowpack;
	float snowfall, snowmelt; 

	/* compute coefficients for snowfall and snowmelt functions */
	fall_b = -1.0 / (fireParams.snw0 - fireParams.snw1); 
	fall_a = 1.0 + fall_b * (-fireParams.snw1);
	melt_a = fireParams.melt_b *  (-fireParams.no_melt);

	/* initialize snowpack and begin daily cycle */   
	snowpack = m_last_yr_snow;
	m_tot_snowfall = 0.;

	for (int i = 0; i<365; i++)
	{
		/* estimate ppt falling as snow */      
		if (m_tmp_doy[i]<=fireParams.snw1) snowfall = m_ppt_doy[i];
		else if (m_tmp_doy[i]>=fireParams.snw0) snowfall = 0.;
		else snowfall = (fall_a + fall_b*m_tmp_doy[i])*m_ppt_doy[i];

		if (!(snowfall<=m_ppt_doy[i]))
		{
			printf("*** snow_cond(): i, snowfall, m_ppt_doy[i] = %d, %f, %f\n", i, snowfall, m_ppt_doy[i]);
			assert(0);
		}

		m_tot_snowfall += snowfall;

		/* determine potential snow melt */      
		if (m_tmp_doy[i]<=fireParams.no_melt) pot_melt = 0.;
		else pot_melt = melt_a + fireParams.melt_b*m_tmp_doy[i];

		snowpack += snowfall; /* add snowfall to snowpack */

		/* determine snowmelt and subtract from snowpack */
		if (snowpack>0.)
		{
			if (pot_melt>snowpack)
			{
				snowmelt = snowpack;
				snowpack = 0.;
			}
			else
			{
				snowmelt = pot_melt; 
				snowpack -= snowmelt;
			}
		}

		m_snow_doy[i] = snowpack;
	} // end of for (i = 0; i<365; i++)

} // end of MC_FireModel::snow_cond()


void MC_FireModel::kbdi()
{
	float Tcurr, Pann, Pcurr, Kprev;
	float Snow;
	float num, dem, factor;
	int F1; 
	float sum_prev, sum_curr;
	float Padj;
	float Kcurr;
	float Kmax;

	Pann = m_ave_ann_ppt*0.0393700787f; // convert from mm to inches 
	Kmax = 0.;

	/* Initialization for start of year */
	/* Initialization for start of run is done by ProcessModelInit(). */    
	Kprev  = stateData.fireState.Kprev;
	F1     = stateData.fireState.F1;
	sum_prev = stateData.fireState.Ksum;

	for (int day = 0; day < 365; day++)
	{ 
		Pcurr = m_ppt_doy[day]*0.0393700787f; // convert from mm to inches
		Tcurr = (m_tmp_doy[day]*9./5.) + 32.; // convert from C to F 
		Snow = m_snow_doy[day]*0.0393700787f; // convert from mm to inches

		if (Tcurr < 50.) Tcurr = 50.;

		/* Calculate drought factor */

		num = .0486 * Tcurr;
		num = exp(num);
		num = .968 * num;
		num = num - 8.30;
		num = (800. - Kprev) * num;

		dem = -.0441 * Pann;
		dem = exp(dem);
		dem = 10.88 * dem;
		dem = 1 + dem;

		factor = (num / dem) * .001;


		/* Adjust Pcurr */         
		if (Snow > 0.0) Padj = Pcurr;            
		else if (Pcurr > 0.0)
		{                    
			sum_curr = sum_prev + Pcurr;

			if (sum_curr<=0.2f) Padj = 0.;                     
			else if (F1==0)
			{
				Padj = sum_curr - .2f;
				F1 = 1;
			}                  
			else Padj = Pcurr;
		}         
		else
		{
			Padj = 0.;
			sum_prev = 0.; 
			F1 = 0;
		}

		sum_prev += Pcurr;


		/* Current day's KDBI */         
		Kcurr = (Kprev - (Padj*100.f)) + factor;      
		if (Kcurr < 0.) Kcurr = 0.;

		m_kbdi_doy[day] = Kcurr;

		if (Kcurr > Kmax) Kmax = Kcurr;

		Kprev = Kcurr; // prepare for next pass thru the loop

	} // end of loop thru the days of the year

	m_kbdi_max = Kmax;
	stateData.fireState.Kprev = Kcurr;
	stateData.fireState.F1 = F1; 
	stateData.fireState.Ksum = sum_prev; 

} // end of MC_FireModel::kbdi()


void MC_FireModel::DFuelMC()
	// For each day of the year, estimates percent moisture content
	// of dead fuels in the 1,10,100, and 1000-hr
	// time lag classes and flammability
{

	// calculate Canadian FWI indices    
	fwi();

	// calculate dead fuel moisture    
	fuel_mc();

	// estimate fine fuel flammability    
	// flammability(); needed for VINCERA fire model option

	// find month with minimum 1000-hr mc    
	min_mc_month();

} // end of MC_FireModel::DFuelMC()


void MC_FireModel::fwi()
	// Calculate Canadian FWI indices.
	// Uses metric units.
{

	int mon, day, n_days, julday;
	float el[13],fl[13],ffmc_old,dmc_old,dc_old;
	float dc = -1.;
	float dmc = -1.;
	float ffmc = -1.;
	float isi,bui,fwi,dsr,t,ra,wmo,ed,ew,z,x,wm,rk,pr,rw,wmi,b,wmr,pe;
	float smi,dr,fm,sf,p,cc,bb; // ,wnd;
	int days_per_month[] = DAYS_PER_MONTH;

	ffmc_old = stateData.fireState.ffmc_prev; 
	dmc_old  =  stateData.fireState.dmc_prev;
	dc_old   = stateData.fireState.dc_prev;

	el[0]=6.5;el[1]=7.5;el[2]=9.0;el[3]=12.8;el[4]=13.9;el[5]=13.9;
	el[6]=12.4;el[7]=10.9;el[8]=9.4;el[9]=8.0;el[10]=7.0;el[11]=6.0;

	fl[0]=-1.6;fl[1]=-1.6;fl[2]=-1.6;fl[3]=0.9;fl[4]=3.8;fl[5]=5.8;
	fl[6]=6.4;fl[7]=5.0;fl[8]=2.4;fl[9]=0.4;fl[10]=-1.6;fl[11]=-1.6;

	julday = 0;

	for (mon = 0; mon < 12; mon++)
	{
		n_days = days_per_month[modelParamsP->southernHemisphereFlag ? (mon + 6)%12 : mon];

		for (day = 0; day < n_days; day++)
		{         
			// wnd = m_wnd_doy[julday] * 3.6;

			// CALC FFMC          
			wmo = 147.2f * (101.f - ffmc_old) / (59.5f + ffmc_old);         
			if (m_ppt_doy[julday]>0.5)
			{ 
				ra = m_ppt_doy[julday] - 0.5;               
				if (wmo > 150.f) wmo += 
					42.5*ra*exp(-100.0/(251-wmo))*(1.0-exp(-6.93/ra)) + 0.0015*(wmo-150)*(wmo-150)*sqrt(ra);
				else wmo += 42.5*ra*exp(-100.0/(251-wmo))*(1.0-exp(-6.93/ra));
			}         
			if (wmo>250.f) wmo = 250.f;

			ed = 0.942*pow(m_rh_doy[julday],0.679) + (11.0*exp((m_rh_doy[julday] - 100.0)/10.0)) +
				0.18*(21.1 - m_tmp_doy[julday])*(1.0-1.0/exp(m_rh_doy[julday]*0.115));
			ew = 0.618*pow(m_rh_doy[julday],0.753) + (10.0*exp((m_rh_doy[julday] - 100.0)/10.0)) +
				0.18*(21.1 - m_tmp_doy[julday])*(1.0 - 1.0/exp(m_rh_doy[julday]*0.115));

			if (wmo<ed && wmo<ew)
			{ 
				z = 0.424*(1.0 - pow(((100.0 - m_rh_doy[julday])/100.0),1.7)) +
					0.0694*sqrt(m_wnd_doy[julday])*(1.0 - pow((100.0 - m_rh_doy[julday])/100.0, 8.0));                  
				x = z*0.581*exp(0.0365*m_tmp_doy[julday]);                  
				wm = ew - (ew - wmo)/pow(10.0, x);
			}         
			else if(wmo > ed)
			{
				z  = 0.424*(1.0 - pow((m_rh_doy[julday]/100.), 1.7)) + 0.0694*sqrt(m_wnd_doy[julday])*(1.0 - pow(m_rh_doy[julday]/100, 8.0));
				x  = z*0.581*exp(0.0365*m_tmp_doy[julday]);
				wm = ed + (wmo - ed)/pow(10.0, x);
			}
			else wm = wmo;

			ffmc = 59.5*(250.0 - wm)/(147.2 + wm);
			if (ffmc>101.0f) ffmc = 101.0f;
			if(ffmc<0.0) ffmc = 0.f;

			assert(!isnan(ffmc));
			m_ffmc_doy[julday] = ffmc;

			// CALC DMC 
			t = m_tmp_doy[julday];
			if (t<-1.1f) t = -1.1;

			rk = 1.894*(t + 1.1)*(100.0 -m_rh_doy[julday])*el[mon]*0.0001;

			if (m_ppt_doy[julday]<=1.5) pr = dmc_old;
			else
			{
				ra = m_ppt_doy[julday];
				rw = 0.92*ra - 1.27;
				wmi = 20.0 + 280.0/exp(0.023*dmc_old);

				if (dmc_old<=33) b = 100.0/(0.5 + 0.3*dmc_old);
				else if (dmc_old <= 65) b = 14.0 - 1.3*log(dmc_old);
				else b = 6.2*log(dmc_old) - 17.2;

				wmr = wmi + 1000.0*rw/(48.77 + b*rw);

				pr = 43.43*(5.6348 - log(wmr - 20));
			}

			if (pr<0) pr = 0;

			dmc = pr + rk;
			if (dmc<0) dmc = 0;
			if (dmc > 250.) dmc = 250.;

			// CALC DC 
			if(m_tmp_doy[julday]<-2.8) t = -2.8;

			pe = (.36*(t + 2.8) + fl[mon])/2.0;
			if (pe<0.0) pe = 0.0;

			if(m_ppt_doy[julday]<=2.8) dr = dc_old;
			else
			{
				ra = m_ppt_doy[julday];
				rw = 0.83*ra - 1.27;
				smi = 800*exp(-dc_old/400);
				dr = dc_old - 400.0*log(1.0 + 3.937*rw/smi);
				if(dr < 0) dr = 0;
			}

			dc = dr + pe;
			if (dc<0) dc = 0;
			if (dc>2000.) dc = 2000.;

			m_dc_doy[julday] = dc;

			// CALC ISI   
			fm  = 147.2*(101.0 - ffmc)/(59.5 + ffmc);
			sf  = 19.115*exp(-0.1386*fm)*(1.0 + pow(fm,5.31)/4.93e07);
			isi = sf*exp(0.05039*m_wnd_doy[julday]);

			// CALC BUI 
			if (dmc==0 && dc==0) bui = 0;
			else bui = 0.8*dc*dmc/(dmc + 0.4*dc);
			if (bui< dmc)
			{ 
				p = (dmc - bui)/dmc;
				cc  = 0.92+pow((0.0114*dmc),1.7);
				bui = dmc - cc*p;
				if (bui<0) bui = 0;
			}

			m_bui_doy[julday] = bui;

			// CALC FWI AND DSR 
			if (bui>80) bb= 0.1*isi*(1000.0/(25.0 + 108.64/exp(0.023*bui)));
			else bb= 0.1*isi*(0.626*pow(bui,0.809) + 2.0);
			if (bb<=1) fwi = bb;
			else fwi = exp(2.72*pow(0.434*log(bb),0.647));

			dsr = 0.0272*pow(fwi, 1.77);

			// SET VALUES FOR NEXT DAY 
			ffmc_old = ffmc; 
			dmc_old  = dmc;  
			dc_old   = dc;
			julday += 1;

		} // END DAILY LOOP 
	} // END MONTHLY LOOP 


	// assign startup vals for next year    
	stateData.fireState.ffmc_prev = ffmc; 
	stateData.fireState.dmc_prev =  dmc;
	stateData.fireState.dc_prev = dc;

	// find annual max ffmc and bui vals for output 

	m_ffmc_ann_max = m_bui_ann_max = 0.;   
	for (day = 0; day < 365; day++)
	{
		if (m_ffmc_doy[day] > m_ffmc_ann_max) m_ffmc_ann_max = m_ffmc_doy[day];         
		if (m_bui_doy[day] > m_bui_ann_max) m_bui_ann_max = m_bui_doy[day];
	}

} // MC_FireModel::end of fwi()


void MC_FireModel::fuel_mc()
	// Calculate daily moisture values of 1-hr, 10-hr, 100-hr, and 1000-hr fuels.   
{

	int i, j;
	float arg;
	float ppt_dur;
	float jdate, phi, decl, daylit;
	float emc_bar; 
	float bnd_100, ymc_100;
	float ymc_1000;
	float bnd_sum, bnd_bar;
	float bnd_1000[365];

	ymc_100 = stateData.fireState.prev_day_mc_100hr;
	ymc_1000 = stateData.fireState.prev_day_mc_1000hr;

	for (i = 0; i<365; i++)
	{ float rh_corr, mc_1hr, mc_10hr, mc_100hr, mc_1000hr, emc_min, emc_max, emc_min_uc;

		m_temp_corr_doy[i] = m_tmax_doy_F[i] + 15.f; // correct daytime temperature       
		rh_corr = m_rhmin_doy[i]*0.87f; // correct daytime relative humidity 

		/* estimate precipitation duration, constrain to 8 hr */
		ppt_dur = floor((m_ppt_doy_in[i]/m_ppt_rat_doy_in_per_hr[i]) + 0.49f);       
		if (ppt_dur>8.0f) ppt_dur = 8.0f;

		if (m_tmin_doy_F[i]<32.f)
		{ // set ppt_dur and rh_corr for frozen fuels 
			ppt_dur = 0.f;
			rh_corr = 100.f;
		}

		/* calculate minimum equilibrium moisture content using corrected tmp and rh */      
		if (rh_corr<10.f) emc_min = 0.03229 + 0.281073 * rh_corr - 0.000578 * m_temp_corr_doy[i] * rh_corr;
		else if (rh_corr<50.f) emc_min = 2.22749 + 0.160107 * rh_corr - 0.014784 * m_temp_corr_doy[i];
		else emc_min = 21.06060 + 0.005565*pow(rh_corr, 2.) - 0.00035*rh_corr*m_temp_corr_doy[i] - 0.483199*rh_corr;

		/* calculate minimum equilibrium moisture content using uncorrected tmp and rh */
		if (m_rhmin_doy[i]<10.f) emc_min_uc = 0.03229 + 0.281073*m_rhmin_doy[i] - 0.000578*m_tmax_doy_F[i]*m_rhmin_doy[i];
		else if (m_rhmin_doy[i]<50.f) emc_min_uc = 2.22749 + 0.160107*m_rhmin_doy[i] - 0.014784*m_tmax_doy_F[i];
		else emc_min_uc = 21.06060 + 0.005565*pow(m_rhmin_doy[i], 2.) - 0.00035*m_rhmin_doy[i]*m_tmax_doy_F[i] - 0.483199*m_rhmin_doy[i];

		/* calculate maximum equilibrium moisture content */
		if (m_rhmax_doy[i]<10.f) emc_max = 0.03229 + 0.281073*m_rhmax_doy[i] - 0.000578*m_tmin_doy_F[i]*m_rhmax_doy[i];
		else if (m_rhmax_doy[i]<50.f) emc_max = 2.22749 + 0.160107*m_rhmax_doy[i] - 0.014784*m_tmin_doy_F[i];
		else emc_max = 21.06060 + 0.005565*pow(m_rhmax_doy[i], 2.) - 0.00035*m_rhmax_doy[i]*m_tmin_doy_F[i] - 0.483199*m_rhmax_doy[i];


		/* 1-HOUR AND 10-HOUR FUEL MOISTURE CONTENT */                                 
		if ((m_ppt_doy_in[i]>0.f) || (m_tmin_doy_F[i]<32.f)) mc_1hr = mc_10hr = 35.; // wet or frozen fuels 
		else
		{ 
			mc_1hr = 1.0329 * emc_min;
			mc_10hr = 1.2815 * emc_min;
		}      
		fuel_data.mc_1hr[i]  = mc_1hr;
		fuel_data.mc_10hr[i] = mc_10hr;      

		/* calculate daylength */      
		jdate = ((float) i) + 1.;
		phi = inputDataP->lat * 0.01745;
		decl = 0.41008 * sin((jdate - 82.) * 0.01745);
		arg = tan(phi) * tan(decl);      
		if (arg < -1.) arg = -1.;
		else if (arg > 1.) arg = 1.;
		daylit = 24. * (1. - acos(arg) / 3.1416);

		/* calculate weighted 24-hour average emc */
		emc_bar = (daylit * emc_min_uc + (24. - daylit) * emc_max) / 24.;

		/* 100-HR MC */
		/* calculate 100-hr boundary condition for current day */
		bnd_100 = ((24. - ppt_dur) * emc_bar + ppt_dur * (0.5 * ppt_dur + 41.)) / 24.;

		mc_100hr = ymc_100 + (bnd_100 - ymc_100) * (1. - 0.87 * exp(-0.24));
		if (mc_100hr > 35.) mc_100hr = 35.;
		fuel_data.mc_100hr[i] = mc_100hr;
		ymc_100 = mc_100hr;

		/* 1000-HR MC */         
		/* calculate 1000-hr boundary condition for current day */         
		bnd_1000[i] = ((24. - ppt_dur) * emc_bar + ppt_dur * (2.7 * ppt_dur + 76.)) / 24.;

		if (i < 7) mc_1000hr = ymc_1000; // set mc for first week of year           
		else
		{ // calculate fuel moisture content 
			/* calculate prior seven-day average boundary condition */            
			bnd_sum = 0.;
			for (j = i; j > i - 7; j--) bnd_sum += bnd_1000[j];
			bnd_bar = bnd_sum / 7.;

			/* calculate fuel moisture content */            
			ymc_1000 = fuel_data.mc_1000hr[i - 7];
			mc_1000hr = ymc_1000 + (bnd_bar - ymc_1000)*(1. - 0.82 * exp(-0.168));            
			if (mc_1000hr > 35.) mc_1000hr = 35.;
		}
		fuel_data.mc_1000hr[i] = mc_1000hr;

	} // end of loop thru 365 days of year

	stateData.fireState.prev_day_mc_100hr = fuel_data.mc_100hr[364];
	stateData.fireState.prev_day_mc_1000hr = fuel_data.mc_1000hr[364];

} // end of MC_FireModel::fuel_mc()


/* Needed for VINCERA fire model option.
   void MC_FireModel::flammability()
   {   
   float fuel_tmp, mc;   
   float qign, chi;
   float pnorm1 = 0.00232;
   float pnorm2 = 0.99767;
   float pnorm3 = 0.0000185;

   for (int doy = 0; doy<365; doy++)
   {      
// calculate heat required to produce ignition 
fuel_tmp = m_temp_corr_doy[doy];; 
mc =  fuel_data.mc_1hr[doy];

qign = 144.5 - 
(0.266 * fuel_tmp) -
(0.00058 * pow(fuel_tmp, 2.0)) -
(0.01 * fuel_tmp * mc) +
(18.54 * (1.0 - exp(-0.151 * mc))) +
(6.4 * mc);

if (qign > 344.) qign = 344.;  


chi = (344.0 - qign) / 10.; // intermediate calculation 

// calculate probability of firestart 
if ((pow(chi, 3.6) * pnorm3)<=pnorm1) m_p_flamm = 0.;
else m_p_flamm = (pow(chi, 3.6) * pnorm3 - pnorm1) * 100. / pnorm2;
if (m_p_flamm < 0.) m_p_flamm = 0.;
else if (m_p_flamm > 100.) m_p_flamm = 100.;

} // end of daily loop

} // end of MC_FireModel::flammability()
*/


void MC_FireModel::min_mc_month()
{
	int day, month, min_day, min_mo, month_length; 
	float min_mc;

	min_mc = 35.;
	min_day = 365;

	// find day for min mc_1000hr 

	for (day = 0; day < 365; day++)
	{      
		if (fuel_data.mc_1000hr[day]<min_mc)
		{
			min_mc = fuel_data.mc_1000hr[day];
			min_day = day;
		} 

	}

	// find month for min mc_1000hr 
	min_mo = -1;
	day = month = 0;
	while (min_mo<0 && month<12)
	{
		month_length = sciFn.days_per_mo[modelParamsP->southernHemisphereFlag ? (month + 6)%12 : month];
		day += month_length;
		if (min_day<day) min_mo = month;
		month++;
	}
	assert(min_mo>=0 && min_mo<12);

	m_min_mc_yr = min_mc;
	m_month_of_min_mc = min_mo;

} // end of MC_FireModel::min_mc_month()


float MC_FireModel::liveFuelMoistureContent(float stress, float mc_min, float mc_max)
	// Estimates percent moisture content of live herbaceous and tree fuel classes
{
	float x, y, mc;

	x = stress * 100.;
	y = 16.99365 + 84.59560 / (1.+ exp(-(x - 39.88160) / -5.19634));
	if (y>100.) y = 100.;
	else if (y<0.) y = 0.;

	mc = ((mc_max - mc_min) * (y / 100.)) + mc_min;
	assert(mc>=mc_min && mc<=mc_max);

	return(mc);
} // end of MC_FireModel::liveFuelMoistureContent()


void MC_FireModel::FuelLoad(int vtype, float clai, int month_length)
	/* FuelLoad: estimates loading in all live and dead
	   fuel classes from CENTURY estimates of per unit
	   weights of carbon pools and various allometric 
	   functions
	   */
{
	float littr, dwood, ltree, dbh;
	float frac_1hr, frac_10hr, frac_100hr, frac_1000hr, depth_ratio, thick_ratio, cl_ratio;
	float * vveg2loadP;
	const float vveg2loadGlobal[][7] = {
		{.00, .00, .00, .00, .0, .000, .0}, // 0 unused
		{.00, .00, .00, .00, .0, .000, .0}, // 1 ice aka barren
		{.25, .25, .25, .25, .4, .022, .8}, // 2 tundra aka alpine
		{.27, .20, .24, .29, .042, .022, .8}, // 3 taiga-tundra
		{.27, .20, .24, .29, .042, .022, .8}, // 4 boreal needleleaf forest
		{.37, .26, .27, .10, .042, .043, .5}, // 5 boreal mixed woodland
		{.27, .20, .24, .29, .042, .022, .8}, // 6 subalpine forest
		{.20, .17, .20, .43, .042, .062, .7}, // 7 maritime needleleaf forest
		{.39, .28, .14, .19, .042, .043, .5}, // 8 temperate needleleaf forest
		{.62, .20, .18, .00, .042, .033, .4}, // 9 temperate deciduous broadleaf forest
		{.37, .26, .27, .10, .042, .043, .5}, // 10 temperate cool mixed forest
		{.45, .37, .16, .02, .042, .043, .5}, // 11 temperate warm mixed forest
		{.61, .31, .06, .02, .042, .062, .4}, // 12 temperate needleleaf woodland
		{.83, .07, .05, .05, .4, .033, .4}, // 13 temperate deciduous broadleaf woodland
		{.83, .07, .05, .05, .4, .033, .4}, // 14 temperate cool mixed woodland
		{.70, .18, .12, .00, .4, .043, .5}, // 15 temperate warm mixed woodland
		{.72, .24, .02, .02, .4, .043, .7}, // 16 C3 shrubland and temperate shrubland
		{.93, .05, .01, .01, 1.0, .022, .8}, // 17 C3 grassland
		{.75, .24, .01, .00, .042, .043, .7}, // 18 temperate desert
		{.39, .28, .14, .19, .042, .043, .5}, // 19 subtropical needleleaf forest
		{.62, .20, .18, .00, .042, .033, .4}, // 20 subtropical deciduous broadleaf forest 
		{.45, .37, .16, .02, .042, .043, .5}, // 21 subtropical evergreen broadleaf forest
		{.45, .37, .16, .02, .042, .043, .5}, // 22 subtropical mixed forest
		{.78, .19, .01, .02, .4, .062, .4}, // 23 subtropical needleleaf woodland
		{.83, .07, .05, .05, .4, .033, .4}, // 24 subtropical deciduous broadleaf woodland
		{.70, .18, .12, .00, .4, .043, .5}, // 25 subtropical evergreen broadleaf woodland
		{.70, .18, .12, .00, .4, .043, .5}, // 26 subtropical mixed woodland
		{.72, .24, .02, .02, .4, .043, .7}, // 27 C4 shrubland and subtropical shrubland
		{.92, .07, .01, .00, 1.0, .022, .8}, // 28 C4 grassland
		{.75, .24, .01, .00, .042, .043, .7}, // 29 subtropical desert
		{.62, .20, .18, .00, .042, .033, .4}, // 30 tropical evergreen broadleaf forest 
		{.83, .07, .05, .05, .4, .033, .4}, // 31 tropical deciduous woodland
		{.83, .07, .05, .05, .4, .033, .4}, // 32 tropical savanna
		{.72, .24, .02, .02, .4, .043, .7}, // 33 tropical shrubland
		{.92, .07, .01, .00, 1.0, .022, .8}, // 34 tropical grassland
		{.75, .24, .01, .00, .042, .043, .7}, // 35 tropical desert
		{.39, .28, .14, .19, .042, .043, .5}, // 36 cool moist needleleaf forest
		{.72, .24, .02, .02, .4, .043, .7}, // 37 unused or Lynx_AgricultureGrazing 
		{.93, .05, .01, .01, 1.0, .022, .8}, // 38 subalpine meadow
		{.00, .00, .00, .00, .0, .000, .0}, // 39 water and wetlands
		{.00, .00, .00, .00, .0, .000, .0}, // 40 natural barren 
		{.00, .00, .00, .00, .0, .000, .0}, // 41 developed
		{.27, .20, .24, .29, .042, .022, .8}, // 42 larch forest 
		{.20, .17, .20, .43, .042, .062, .7}, // 43 Sitka spruce zone, from 7 maritime needleleaf forest
		{.20, .17, .20, .43, .042, .062, .7}, // 44 western hemlock zone, from 7 maritime needleleaf forest
		{.20, .17, .20, .43, .042, .062, .7}, // 45 Pacific silver fir zone, 7 maritime needleleaf forest
		{.27, .20, .24, .29, .042, .022, .8}, // 46 mountain hemlock zone, from 6 subalpine forest
		{.27, .20, .24, .29, .042, .022, .8}, // 47 subalpine fir zone, from 6 subalpine forest
		{.37, .26, .27, .10, .042, .043, .5}, // 48 subalpine parkland zone, from 5 boreal mixed woodland
		{.39, .28, .14, .19, .042, .043, .5}, // 49 cool dry needleleaf forest
		// {.61, .31, .06, .02, .042, .062, .4}, // 50 XERIC_NEEDLELEAF_WOODLANDveg
	}; // end of vveg2loadGlobal[][7]

	littr = m_mlittr + m_slittr;
	dwood = littr + m_dstnd + m_dwod1 + m_dwod100;
	ltree = m_lleaf + m_lwod1 + m_lwod100;

	/* get veg type dependent parameters */
	switch (runParamsP->baseCalibration)
	{
		case mc2WCR:
		case mc2ConUS:
		case mc2GLOBAL:
		case mc2ConUS_LC:
		case mc2California:
			vveg2loadP = (float *)&(vveg2loadGlobal[0][0]);
			break;
			/*      
				case mc2CA08: 
				case mc2YOSE:
				assert(0); // vveg2loadP = vveg2loadCA08; 
				break;

				case mc2VINCERA:
				assert(0); // vveg2loadP = vveg2loadVINCERA;
				break;
				*/
		default: assert(0); break;
	}
	frac_1hr = *(vveg2loadP + vtype*7 + 0);
	frac_10hr = *(vveg2loadP + vtype*7 + 1);
	frac_100hr = *(vveg2loadP + vtype*7 + 2);
	frac_1000hr = *(vveg2loadP + vtype*7 + 3);
	assert(sciFn.close_enough(frac_1hr + frac_10hr + frac_100hr + frac_1000hr, 1.0, .00001) || frac_1hr==0.);

	depth_ratio = *(vveg2loadP + vtype*7 + 4);
	thick_ratio = *(vveg2loadP + vtype*7 + 5);
	cl_ratio = *(vveg2loadP + vtype*7 + 6);   

	/* Are there trees? */
	if (ltree > 0.0)
	{ // There are trees. 
		float lwood, lstem;

		/* estimate dbh and height */
		if (modelParamsP->code_flags[ALT_TREE_ALLOMETRY_CODE_FLAG]) 
			sciFn.allometry_DavidKing(m_lwod1 + m_lwod100, stateData.deciduous, &m_tree_ht_m, &dbh);
		else sciFn.tree_dim(clai, ltree, stateData.deciduous, &m_tree_ht_m, &dbh);

		m_bark_thick = dbh * thick_ratio; // estimate bark thickness using vveg thickness to dbh ratio
		m_crown_len_m = m_tree_ht_m * cl_ratio; // estimate crown length using vveg length to height ratio 
		lwood = m_lwod1 + m_lwod100;
		sciFn.live_wood(dbh, stateData.deciduous, lwood, &m_lbranch, &lstem); // estimate loads of live branch and stem wood classes 
	}
	else dbh = m_tree_ht_m = m_crown_len_m = m_bark_thick = m_lbranch = 0.0; // No trees.

	/* dead fuels */   
	if (dwood>0.0) 
	{
		// dead_wood(data_point, mo, yr);
		float tot_dfuel;

		tot_dfuel = littr + m_dstnd + m_dwod1 + m_dwod100;      
		if (tot_dfuel > 10000.) tot_dfuel = 10000.;

		if (modelParamsP->code_flags[ALT_FUEL_LOAD_CODE_FLAG])
		{
			m_d1hr = littr + m_dstnd;
			m_d10hr = 0.5*m_dwod1;
			m_d100hr = 0.5*m_dwod1;
			m_d1000hr = m_dwod100; 
		}
		else 
		{
			m_d1hr    = frac_1hr    * tot_dfuel;
			m_d10hr   = frac_10hr   * tot_dfuel;
			m_d100hr  = frac_100hr  * tot_dfuel;
			m_d1000hr = frac_1000hr * tot_dfuel;
		}
		m_fine_fuel_frac = (m_lgras + m_dstnd) / (m_lgras + m_dstnd + m_dwod1);
		// data_point->vemap2_mo.fine_frac[mo] = grass_frac;
	}
	else m_fine_fuel_frac = m_d1hr = m_d10hr = m_d100hr = m_d1000hr = 0.0;

	/* estimate fuel bed depth */
	m_fuel_depth = sciFn.bed_depth(m_lgras + m_dstnd + m_dwod1 + m_dwod100, depth_ratio);

	// Get daily values of dead fuel classes, converting from g m-2 to lbs ft-2
	get_daily(stateData.fireState.prev_1hr*0.000204849, m_d1hr*0.000204849, month_length, m_w1p_dom);
	get_daily(stateData.fireState.prev_10hr*0.000204849, m_d10hr*0.000204849, month_length, m_w10_dom);
	get_daily(stateData.fireState.prev_100hr*0.000204849, m_d100hr*0.000204849, month_length, m_w100_dom);
	get_daily(stateData.fireState.prev_1000hr*0.000204849, m_d1000hr*0.000204849, month_length, m_w1000_dom);
	get_daily(stateData.fireState.prev_l1hr*0.000204849, m_lgras*0.000204849, month_length, m_wherbp_dom);

	/* create pseudo_daily values for fuel moisture and loading */
	get_daily(stateData.fireState.prev_depth, m_fuel_depth, month_length, m_depth_dom);
	get_daily(stateData.fireState.prev_mc_grass, m_mc_grass, month_length, m_mc_grass_dom);
	get_daily(stateData.fireState.prev_mc_tree, m_mc_tree, month_length, m_mc_tree_dom);

} // end of MC_FireModel::FuelLoad()


/*
   void MC_FireModel::erc_G(int month_length)
   { 
   float sum_erc;   

// Calculate ERC using fuel model G   
sum_erc = 0.;

// surface to volume ratios  
m_sg1    = 2000.;
m_sg10   = 109.;
m_sg100  = 30.;
m_sg1000 = 8.;
m_sgherb = 1500.;
m_sgwood = 2000.;
m_mxd_frac = 0.25; // dead fuel moisture of extinction 
m_hl = m_hd = 8000.; // dead and live fuel heat of combustion 
for (int dom = 0; dom<month_length; dom++)
{      
//assign(data_point, mo, dom);

// FUEL MODEL G  // from Jim Lenihan, 7/12/2011
// fuel loadings in tons per acre 

m_w1p_dom[dom]    = 2.5 // tons acre-1  * .0459137 // (lbs ft-2)/(tons acre-1);
m_w10_dom[dom]    = 2.0 // tons acre-1  * .0459137 // (lbs ft-2)/(tons acre-1);
m_w100_dom[dom]   = 5.0 // tons acre-1  * .0459137 // (lbs ft-2)/(tons acre-1);
m_w1000_dom[dom]  = 12.0 // tons acre-1 * .0459137 // (lbs ft-2)/(tons acre-1);
m_wherbp_dom[dom] = 0.5 // tons acre-1 * .0459137 // (lbs ft-2)/(tons acre-1);
m_wwood_dom[dom]  = 0.5 // tons acre-1 * .0459137 // (lbs ft-2)/(tons acre-1);
m_depth_dom[dom] = 1.0; // fuel bed depth in ft 

// prelim(dom, mo, yr); // make some preliminary calculations of spread model input vars 

// FIRE BEHAVIOR ROUTINES //
// spread(data_point, dom, mo, yr); // rate of spread 
// release(dom); // energy release 
// intensity(data_point, dom, mo, yr); // fireline intensity and related measures 
// sum_erc += .fire_behav.erc_daily[31];
} // end of daily loop using fuel model G
m_erc_g = NC_FILL_FLOAT; // sum_erc/month_length;
// end of logic to calculate reference value of ERC using fuel model G

} // end of MC_FireModel::erc_G()
*/

void MC_FireModel::get_daily(float prev_mo_val, float curr_mo_val, int mo_len, float daily_vals[])
	// Given the values of a variable on the last day of the previous month and 
	// the last day of the current month, estimate the value on every day of the
	// current }month by linear interpolation.
{
	float diff, delt;

	assert(28<=mo_len && mo_len<=31);

	diff = curr_mo_val - prev_mo_val;
	delt = diff/mo_len;

	for (int day = 0; day<mo_len; day++) daily_vals[day] = prev_mo_val + (day + 1)*delt;

} // end of MC_FireModel::get_daily()


void MC_FireModel::FireOccur(int mo, int yrs_since_fire, int vtype, int * fire_domP, int * fire_doyP)
	// Returns dom and doy of fire.  If no fire, returns NO_FIRE for both.
{
	int mo_len, doy;
	int fire_dom, fire_doy;
	bool fire_set_flag;

	fire_set_flag = runParamsP->fire_set_interval>0 && runParamsP->fire_set_interval<=yrs_since_fire;
	fuel_characteristics(vtype);

	mo_len = sciFn.days_per_mo[modelParamsP->southernHemisphereFlag ? ((mo + 6)%12) : mo];
	doy = sciFn.doy_of_dom0(mo, modelParamsP->southernHemisphereFlag) - 1;   
	int dom = -1;
	fire_dom = fire_doy = NO_FIRE;
	while ((fire_dom==NO_FIRE) && (dom<(mo_len - 1)))
	{
		float ros, sgbrt;
		dom++;
		doy++;

		ros_dom(dom, doy, &ros, &sgbrt);
		if (ros<=0.) continue;

		if ((fire_set_flag && doy==(runParamsP->fire_set_jday - 1)) // doy Jan 1 is 0; fire_set_jday Jan 1 is 1
				|| (m_ffmc_doy[doy]>=m_ffmc_threshold && m_bui_doy[doy]>=m_bui_threshold))
		{ // Simulate a fire today.
			fire_dom = dom;
			fire_doy = doy;
			m_fire_ros = ros;
			m_fire_sgbrt = sgbrt;
			erc_dom(dom, doy, &m_fire_erc, &m_fire_tau, &m_fire_ire); 
			float fd = m_fire_ros*m_fire_tau; // depth of flaming zone in ft
			m_fire_fli = m_fire_ire*(fd/60.);
		}         
	}

	*fire_domP = fire_dom;
	*fire_doyP = fire_doy;

} // end of MC_FireModel::FireOccur()


void MC_FireModel::ros_dom(int dom, int doy, float * rosP, float * sgbrtP)
{
	float ros = 0.;
	float w1p, w10, w100, w1000, wherbp, wwood; // lbs per sq. ft

	// Make local copies because they may be adjusted below.
	w1p = m_w1p_dom[dom];
	w10 = m_w10_dom[dom];;
	w100 = m_w100_dom[dom];
	w1000 = m_w1000_dom[dom];
	wherbp = m_wherbp_dom[dom];
	wwood = 0.0; // m_wwood_dom[dom];

	// total dead, total live, and total fuel load 
	float wtotd, wtotl, wtot; // lbs ft-2
	wtotd = w1p + w10 + w100 + w1000;
	wtotl = wherbp + wwood;
	wtot = wtotd + wtotl;
	if (wtot<0.1) 
	{
		*rosP = 0.;
		*sgbrtP = NC_FILL_FLOAT; // token for no data
		return;
	}

	// 1988 revision: duff drying with drought    
	if (m_kbdi_doy[doy] > 100.)      
	{ float w1f, w10f, w100f, w1000f, unit_incr, tot_incr;
		w1f    =w1p / wtotd;
		w10f   = w10 / wtotd;
		w100f  = w100 / wtotd;
		w1000f = w1000 / wtotd;

		unit_incr = w1p / 700.;
		tot_incr = (m_kbdi_doy[doy] - 100.) * unit_incr;

		w1p   += (w1f    * tot_incr);
		w10   += (w10f   * tot_incr);
		w100  += (w100f  * tot_incr);
		w1000 += (w1000f * tot_incr);

		wtotd = w1p + w10 + w100 + w1000;
		wtot = wtotd + wtotl;
	}

	// Calculate surface area of each class.
	float sa1, sa10, sa100, saherb, sawood, sadead, salive, rhol, rhod;
	rhol = rhod = 32.;
	sa1 = (w1p/rhod) * m_sg1;
	sa10 = (w10/rhod) * m_sg10;
	sa100 = (w100/rhod) * m_sg100;
	saherb = (wherbp/rhol) * m_sgherb;
	sawood = (wwood/rhol) * m_sgwood;
	sadead = sa1 + sa10 + sa100;
	salive = saherb + sawood;

	// Calculate weighting factors of each fuel class.
	float f1, f10, f100, fherb, fwood, fdead, flive;
	if (sadead<=0.) f1 = f10 = f100 = 0.;
	else
	{
		f1 = sa1 / sadead;
		f10 = sa10 / sadead;
		f100 = sa100 / sadead;
	}

	if (salive<=0.) fherb = fwood = 0.;
	else 
	{ 
		fherb = saherb / salive; 
		fwood = sawood / salive;
	}

	if ((sadead + salive) <= 0.0) fdead = flive = 0.;
	else
	{
		fdead = sadead / (sadead + salive);
		flive = salive / (sadead + salive);
	}


	// net fuel loading of each fuel class 
	float std, stl, w1n, w10n, w100n, wherbn, wwoodn; // lbs per sq ft
	std = stl = 0.0555;
	w1n = w1p * (1. - std);
	w10n = w10 * (1. - std);
	w100n = w10 * (1. - std);
	wherbn = wherbp * (1. - stl);
	wwoodn = wwood * (1. - stl);

	// weighted net loadings of dead and live fuels 
	float wdeadn, wliven;   
	wdeadn = (f1 * w1n) + (f10 * w10n) + (f100 * w100n);
	wliven = (fwood * wwoodn) + (fherb * wherbn);

	// dead and live fuel characteristic surface area-to-volume ratios 
	float sgbrd, sgbrl;
	sgbrd = (f1 * m_sg1) + (f10 * m_sg10) + (f100 * m_sg100);
	sgbrl = (fherb * m_sgherb) + (fwood * m_sgwood);

	float sgbrt = (fdead * sgbrd) + (flive * sgbrl); // characteristic surface area-to-volume ratio 

	float betop = 3.348 * pow(sgbrt, -0.8189); // optimum packing ratio 

	/* packing ratio */
	float depth_in_feet = m_depth_dom[dom]*3.28;
	float rhobed, rhobar, betbar;
	rhol = rhod = 32.;
	rhobed = (wtot - w1000) / depth_in_feet;
	rhobar = ((wtotl * rhol) + (wtotd * rhod)) / wtot;
	betbar = rhobed / rhobar;

	// float propr_rat = betbar / betop; // calculate packing ratio : optimum packing ratio 
	float gmamx = (pow(sgbrt, 1.5))/(495. + 0.0594 * pow(sgbrt, 1.5)); // maximum reaction velocity

	// optimum reaction velocity
	float ad = 133.*pow(sgbrt, -0.7913);  
	float gmaop = gmamx * pow((betbar/betop), ad) * exp(ad * (1.0 - betbar/betop));

	// no wind propagating flux ratio 
	float zeta = exp((0.792 + 0.681 * pow(sgbrt, 0.5)) * (betbar + 0.1))/(192. + 0.2595 * sgbrt);

	/* heating numbers of each fuel class */
	float hn1, hn10, hn100, hnherb, hnwood;          
	hn1 = w1n * exp(-138. / m_sg1); 
	hn10 = w10n * exp(-138. / m_sg10);
	hn100 = w100n * exp(-138. / m_sg100);
	hnherb = wherbn * exp(-500. / m_sgherb);
	hnwood = wwoodn * exp(-500. / m_sgwood);

	/* ratio of dead-to-live fuel heating numbers */             
	float wrat = ((hnherb + hnwood) <= 0.) ? -1.: (hn1 + hn10 + hn100) / (hnherb + hnwood);

	// convert moisture contents from % to fractions
	float mc_1hr_frac, mc_10hr_frac, mc_100hr_frac, mc_grass_frac, mc_tree_frac;
	mc_1hr_frac = fuel_data.mc_1hr[doy] * .01;
	mc_10hr_frac = fuel_data.mc_10hr[doy] * .01;
	mc_100hr_frac = fuel_data.mc_100hr[doy]*0.01;
	mc_grass_frac = m_mc_grass*0.01;
	mc_tree_frac = m_mc_tree*0.01;

	// weighted dead_fuel moisture content for live-fuel extinction moisture 
	float mclfe = ((mc_1hr_frac * hn1) + (mc_10hr_frac * hn10) + (mc_100hr_frac * hn100)) 
		/ (hn1 + hn10 + hn100);

	// moisture extinction of live fuels 
	if (wrat <= 0.) m_mxl_frac = -1.;
	else
	{
		m_mxl_frac = 2.9 * wrat * (1. - mclfe / m_mxd_frac) - 0.226;
		if (m_mxl_frac<m_mxd_frac) m_mxl_frac = m_mxd_frac;
	}      

	// weighted moisture content of dead and live fuels
	float wtmcd, wtmcl;
	wtmcd = (f1 * mc_1hr_frac) + (f10 * mc_10hr_frac) + (f100 * mc_100hr_frac);
	wtmcl = (fherb * mc_grass_frac) + (fwood * mc_tree_frac);

	// mineral damping coefficent of live and dead fuels    
	float sd, sl, etasd, etasl;
	sd = sl = 0.01;
	etasd = 0.174 * pow(sd, -0.19);
	etasl = 0.174 * pow(sl, -0.19);

	// moisture damping coefficents of dead and live fuels 
	float dedrt, livrt, etamd, etaml;   
	dedrt = wtmcd/m_mxd_frac;
	livrt = wtmcl/m_mxl_frac;
	etamd = 1. - 2.59 * dedrt + 5.11 * pow(dedrt, 2.) - 3.52 * pow(dedrt, 3.);
	if (etamd < 0.) etamd = 0.;
	else if (etamd > 1.) etamd = 1.;
	etaml = 1. - 2.59 * livrt + 5.11 * pow(livrt, 2.) - 3.52 * pow(livrt, 3.);
	if (etaml < 0.) etaml = 0.;
	else if (etaml > 1.) etaml = 1.;

	// wind effect multiplier  
	float wnd_mph, b_eff, c_var, e_eff, ufact, phiwnd;  
	wnd_mph = sciFn.m_per_sec_to_mph(m_wnd_doy[doy]); // m_wnd_doy[doy]*2.232; // Convert from m/s to mph.
	b_eff = 0.02526 * pow(sgbrt, 0.54);
	c_var = 7.47 * exp(-0.133 * pow(sgbrt, 0.55));
	e_eff = 0.715 * exp(-3.59 * pow(10.,-4.) * sgbrt);
	ufact = c_var * pow((betbar / betop), (e_eff * -1.));
	phiwnd = ufact * pow((wnd_mph * 88. * m_wndfac), b_eff);

	// slope effect multiplier coefficient 
	float slpfct, phislp;
	if (fireParams.slp <= 25.) slpfct = 0.267;
	else if (fireParams.slp <= 40.) slpfct = 0.533;
	else if (fireParams.slp <= 55.) slpfct = 1.068;
	else if (fireParams.slp <= 75.) slpfct = 2.134;
	else slpfct = 4.273;
	phislp = slpfct * pow(betbar, -0.3);

	// reaction intensity 
	float ir = gmaop * ((wdeadn * m_hd * etasd * etamd) + (wliven * m_hl * etasl * etaml));

	// correction to phiwnd for high winds 
	if ((wnd_mph * 88. * m_wndfac) > (0.9 * ir)) phiwnd = ufact * pow((0.9 * ir), b_eff);

	// heat sink 
	float htsink = rhobed *  
		(fdead * 
		 ((f1    * exp(-138. / m_sg1)    * (250. + 1116 * mc_1hr_frac)) +
		  (f10   * exp(-138. / m_sg10)   * (250. + 1116 * mc_10hr_frac)) +
		  (f100  * exp(-138. / m_sg100)  * (250. + 1116 * mc_100hr_frac)))) +
		(flive *
		 ((fherb * exp(-138. / m_sgherb) * (250. + 1116 * mc_grass_frac)) +
		  (fwood * exp(-138. / m_sgwood) * (250. + 1116 * mc_tree_frac))));

	// rate of spread in ft/min       
	ros = ir * zeta * (1. + phislp + phiwnd) / htsink;

	*rosP = ros;
	*sgbrtP = sgbrt;

} // end of MC_FireModel:: ros_dom()


void MC_FireModel::erc_dom(int dom, int doy, float * ercP, float * tauP, float *ireP)
{
	/* weighting factors of each fuel class */
	float wtotd = m_w1p_dom[dom] + m_w10_dom[dom] + m_w100_dom[dom] + m_w1000_dom[dom];
	float wtotl = m_wherbp_dom[dom]; // + m_wwood_dom[dom];
	float f1e, f10e, f100e, f1000e, fherbe, fwoode;
	if (wtotd <= 0.) f1e = f10e = f100e = f1000e = 0.;   
	else
	{ 
		f1e = m_w1p_dom[dom] / wtotd;
		f10e = m_w10_dom[dom] / wtotd;
		f100e = m_w100_dom[dom] / wtotd;
		f1000e = m_w1000_dom[dom] / wtotd;
	}   
	if (wtotl <= 0.) fherbe = fwoode = 0.;
	else
	{
		fherbe = m_wherbp_dom[dom] / wtotl;
		fwoode = 0; // m_wwood_dom[dom] / wtotl;
	}

	/* weighting factors of dead and live fuels */ 
	float wtot = wtotl + wtotd;  
	float fdeade, flivee;
	if (wtot <= 0.0) fdeade = flivee = 0.;   
	else
	{
		fdeade = wtotd / wtot;
		flivee = wtotl / wtot;
	}

	/* net loadings of dead and live fuels */
	float wdedne, wlivne;   
	wdedne = wtotd * (1.0 - c_std);
	wlivne = wtotl * (1.0 - c_stl);

	/* dead and live fuel characteristic surface area-to-volume ratios */
	float sgbrde, sgbrle;
	sgbrde = (f1e * m_sg1) + (f10e * m_sg10) + (f100e * m_sg100) + (f1000e * m_sg1000);
	sgbrle = (fherbe * m_sgherb) + (fwoode * m_sgwood);

	float sgbrte = (fdeade * sgbrde) + (flivee * sgbrle); /* characteristic surface area-to-volume ratio */
	float betope = 3.348 * pow(sgbrte, -0.8189); /* optimum packing ratio */

	/* packing ratio */
	float depth_in_feet = m_depth_dom[dom]*3.28;
	float rhol, rhod, rhobed, rhobar, betbar;
	rhol = rhod = 32.;
	rhobed = (wtot - m_w1000_dom[dom]) / depth_in_feet;
	rhobar = ((wtotl * rhol) + (wtotd * rhod)) / wtot;
	betbar = rhobed / rhobar;

	// float propr_rate = betbar / betope; /* calculate packing ratio : optimum packing ratio */
	float gmamxe = pow(sgbrte, 1.5) / (495. + 0.0594 * pow(sgbrte, 1.5)); /* maximum reaction velocity */
	float ade = 133. * pow(sgbrte, -0.7913); /* optimum reaction velocity */
	float gmaope = gmamxe * pow((betbar / betope), ade) * exp(ade * (1.0 - betbar / betope));

	// convert moisture contents from % to fractions
	float mc_1hr_frac, mc_10hr_frac, mc_100hr_frac, mc_1000hr_frac, mc_grass_frac, mc_tree_frac;
	mc_1hr_frac = fuel_data.mc_1hr[doy] * .01;
	mc_10hr_frac = fuel_data.mc_10hr[doy] * .01;
	mc_100hr_frac = fuel_data.mc_100hr[doy]*0.01;
	mc_1000hr_frac = fuel_data.mc_1000hr[doy]*0.01;
	mc_grass_frac = m_mc_grass*0.01;
	mc_tree_frac = m_mc_tree*0.01;

	/* weighted moisture content of dead and live fuels */
	float wtmcde = (f1e * mc_1hr_frac) + (f10e * mc_10hr_frac)
		+ (f100e * mc_100hr_frac) + (f1000e * mc_1000hr_frac);
	float wtmcle = (fherbe * mc_grass_frac) + (fwoode * mc_tree_frac);

	/* moisture damping coefficents of dead and live fuels */   
	float dedrte = (wtmcde / m_mxd_frac);
	float livrte = (wtmcle / m_mxl_frac);

	float etamde = 1. - 2.0 * dedrte + 1.5 * pow(dedrte, 2.) - 0.5 * pow(dedrte, 3.);   
	if (etamde < 0.) etamde = 0.;  
	else if (etamde > 1.) etamde = 1.;  

	float etamle = 1. - 2.0 * livrte + 1.5 * pow(livrte, 2.) - 0.5 * pow(livrte, 3.);
	if (etamle < 0.) etamle = 0.;
	else if (etamle > 1.) etamle = 1.;

	float ire = gmaope*((fdeade*wdedne*m_hd*c_etasd*etamde) + (flivee*wlivne*m_hl*c_etasl*etamle)); // reaction intensity
	float tau = 384./m_fire_sgbrt; // residence time of flaming front in minutes
	float erc = 0.04*ire*tau; // energy release component

	*ercP = erc;
	*tauP = tau;
	*ireP = ire;

} // end of MC_FireModel::erc_dom()


void MC_FireModel::fuel_characteristics(int vtype)
{
	float * vveg2fuelP;
	const float vveg2fuelGlobal[][10] = {
		{.00, .00, .00, .00, .0, .000, .0, .0, .0, .0}, // 0 unused
		{.00, .00, .00, .00, .0, .000, .0, .0, .0, .0}, // 1 ice aka barren
		{1959., 109., 30., 8., 1462., 1950., 8001., 8001., 30., 0.5}, // 2 tundra aka alpine
		{1852., 109., 30., 8., 1470., 1984., 8039., 8039., 30., 0.4}, // 3 taiga-tundra
		{1852., 109., 30., 8., 1470., 1984., 8039., 8039., 30., 0.4}, // 4 boreal needleleaf forest
		{1634., 109., 30., 8., 1478., 1967., 8026., 8026., 30., 0.4}, // 5 boreal mixed woodland
		{1852., 109., 30., 8., 1470., 1984., 8039., 8039., 30., 0.4}, // 6 subalpine forest
		{1960., 109., 30., 8., 1488., 2045., 8068., 8068., 23., 0.4}, // 7 maritime needleleaf forest
		{1937., 109., 30., 8., 1489., 2120., 8053., 8053., 21., 0.4}, // 8 temperate needleleaf forest
		{1730., 109., 30., 8., 1499., 2003., 8006., 8006., 30., 0.5}, // 9 temperate deciduous broadleaf forest
		{1634., 109., 30., 8., 1478., 1967., 8026., 8026., 30., 0.4}, // 10 temperate cool mixed forest
		{1670., 109., 30., 8., 1499., 2012., 8210., 8210., 30., 0.4}, // 11 temperate warm mixed forest
		{2072., 109., 30., 8., 1451., 2059., 8293., 8293., 16., 0.6}, // 12 temperate needleleaf woodland
		{1906., 109., 30., 8., 1442., 1909., 8016., 8016., 17., 0.6}, // 13 temperate deciduous broadleaf woodland
		{1906., 109., 30., 8., 1442., 1909., 8016., 8016., 17., 0.6}, // 14 temperate cool mixed woodland
		{1433., 109., 30., 8., 1386., 2000., 8680., 8680., 15., 0.6}, // 15 temperate warm mixed woodland
		{2326., 109., 30., 8., 1497., 1978., 8020., 8020., 16., 0.6}, // 16 C3 shrubland
		{2020., 109., 30., 8., 1498., 2021., 8014., 8014., 16., 0.6}, // 17 C3 grassland
		{2425., 109., 30., 8., 1488., 1750., 8072., 8072., 15., 0.6}, // 18 temperate desert
		{1937., 109., 30., 8., 1489., 2120., 8053., 8053., 21., 0.5}, // 19 subtropical needleleaf forest
		{1730., 109., 30., 8., 1499., 2003., 8006., 8006., 30., 0.5}, // 20 subtropical deciduous broadleaf forest 
		{1670., 109., 30., 8., 1499., 2012., 8210., 8210., 30., 0.4}, // 21 subtropical evergreen broadleaf forest
		{1670., 109., 30., 8., 1499., 2012., 8210., 8210., 30., 0.4}, // 22 subtropical mixed forest
		{2232., 109., 30., 8., 1500., 2023., 8001., 8001., 16., 0.6}, // 23 subtropical needleleaf woodland
		{1906., 109., 30., 8., 1442., 1909., 8016., 8016., 17., 0.6}, // 24 subtropical deciduous broadleaf woodland
		{1433., 109., 30., 8., 1386., 2000., 8680., 8680., 15., 0.6}, // 25 subtropical evergreen broadleaf woodland
		{1433., 109., 30., 8., 1386., 2000., 8680., 8680., 15., 0.6}, // 26 subtropical mixed woodland
		{2326., 109., 30., 8., 1497., 1978., 8020., 8020., 16., 0.6}, // 27 C4 shrubland
		{2040., 109., 30., 8., 1495., 2003., 8028., 8028., 15., 0.6}, // 28 C4 grassland
		{2425., 109., 30., 8., 1488., 1750., 8072., 8072., 15., 0.6}, // 29 subtropical desert
		{1730., 109., 30., 8., 1499., 2003., 8006., 8006., 30., 0.5}, // 30 tropical evergreen broadleaf forest 
		{1906., 109., 30., 8., 1442., 1909., 8016., 8016., 17., 0.6}, // 31 tropical deciduous woodland
		{1906., 109., 30., 8., 1442., 1909., 8016., 8016., 17., 0.6}, // 32 tropical savanna
		{2326., 109., 30., 8., 1497., 1978., 8020., 8020., 16., 0.6}, // 33 tropical shrubland
		{2040., 109., 30., 8., 1495., 2003., 8028., 8028., 15., 0.6}, // 34 tropical grassland
		{2425., 109., 30., 8., 1488., 1750., 8072., 8072., 15., 0.6}, // 35 tropical desert
		{1937., 109., 30., 8., 1489., 2120., 8053., 8053., 21., 0.4}, // 36 cool moist needleleaf forest
		{2326., 109., 30., 8., 1497., 1978., 8020., 8020., 16., 0.6}, // 37 unused or Lynx_AgricultureGrazing 
		{2500., 109., 30., 8., 1500., 2500., 8000., 8000., 15., 0.4}, // 38 subalpine meadow
		{.00, .00, .00, .00, .0, .000, .0, .0, .0, .0}, // 39 water and wetlands
		{.00, .00, .00, .00, .0, .000, .0, .0, .0, .0}, // 40 natural barren 
		{.00, .00, .00, .00, .0, .000, .0, .0, .0, .0}, // 41 developed
		{1852., 109., 30., 8., 1470., 1984., 8039., 8039., 30., 0.4}, // 42 larch forest 
		{1960., 109., 30., 8., 1488., 2045., 8068., 8068., 23., 0.4}, // 43 Sitka spruce zone, from 7 maritime needleleaf forest
		{1960., 109., 30., 8., 1488., 2045., 8068., 8068., 23., 0.4}, // 44 western hemlock zone, from 7 maritime needleleaf forest
		{1960., 109., 30., 8., 1488., 2045., 8068., 8068., 23., 0.4}, // 45 Pacific silver fir zone, from 7 maritime needleleaf forest
		{1852., 109., 30., 8., 1470., 1984., 8039., 8039., 30., 0.4}, // 46 mountain hemlock zone, from 6 subalpine forest
		{1852., 109., 30., 8., 1470., 1984., 8039., 8039., 30., 0.4}, // 47 subalpine fir zone, from 6 subalpine forest
		{1634., 109., 30., 8., 1478., 1967., 8026., 8026., 30., 0.4}, // 48 subalpine parkland zone, from 5 boreal mixed woodland
		{1937., 109., 30., 8., 1489., 2120., 8053., 8053., 21., 0.4}, // 49 cool dry needleleaf forest
		// {2072., 109., 30., 8., 1451., 2059., 8293., 8293., 16., 0.6}, // 50 XERIC_NEEDLELEAF_WOODLANDveg
	}; // end of vveg2fuelGlobal[][10]

	/* get veg type dependent parameters */
	switch (runParamsP->baseCalibration)
	{
		// case mc2WWETAC:
		// case mc2US50km:
		case mc2WCR:
		case mc2GLOBAL:
		case mc2ConUS:
		case mc2ConUS_LC:
		case mc2California:
			vveg2fuelP = (float *)&(vveg2fuelGlobal[0][0]);
			break;
			/*      
				case mc2CA08: 
				case mc2YOSE:
				assert(0); // vveg2fuelP = vveg2loadCA08; 
				break;

				case mc2VINCERA:
				assert(0); // vveg2loadP = vveg2loadVINCERA;
				break;
				*/      
		default: assert(0); break;
	}
	m_sg1 = *(vveg2fuelP + vtype*10 + 0);
	m_sg10 = *(vveg2fuelP + vtype*10 + 1);
	m_sg100 = *(vveg2fuelP + vtype*10 + 2);
	m_sg1000 = *(vveg2fuelP + vtype*10 + 3);
	m_sgwood = *(vveg2fuelP + vtype*10 + 4);
	m_sgherb = *(vveg2fuelP + vtype*10 + 5);
	m_hl = *(vveg2fuelP + vtype*10 + 6);
	m_hd = *(vveg2fuelP + vtype*10 + 7);
	m_mxd_frac = *(vveg2fuelP + vtype*10 + 8)*0.01;
	m_wndfac = *(vveg2fuelP + vtype*10 + 9);

} // end of MC_FireModel::fuel_characteristics(int vtype)


float MC_FireModel::part_burn(int doy, int dom, int burn_count, bool try2suppressFlag, int vtype)
	// Returns fraction of cell which burned in the current year
{
	if (vtype>MAX_VTYPE) 
	{ 
		assert(vtype==50); 
		vtype = TEMPERATE_EN_WOODLANDveg; 
	}

	float part;
	float max_fri, min_fri, curr_fri;
	float dscalar;

	//   if (biogeog_option==VINCERA_BIOGEOG) { return(part_burnVINCERA(data_point, yr, mo, day)); };

	// get min and max fire return intervals for current veg class 
	fire_return_interval(vtype, &min_fri, &max_fri);

	// calc curr mfri   
	if (m_fine_fuel_frac < .70) dscalar = (m_bui_doy[doy] - m_bui_threshold) / m_bui_threshold;
	else dscalar = (m_ffmc_doy[doy] - m_ffmc_threshold) / m_ffmc_threshold;

	curr_fri = max_fri - ((max_fri - min_fri) * dscalar);
	if (curr_fri < min_fri) curr_fri = min_fri;
	else if (curr_fri > max_fri) curr_fri = max_fri;
	// assert(min_fri<=curr_fri && curr_fri<=max_fri && curr_fri>0.);
	if (!(min_fri<=curr_fri && curr_fri<=max_fri && curr_fri>0.))
	{
		printf("*** part_burn(): min_fri, curr_fri, max_fri, vtype = %f, %f, %f, %d\n",
				min_fri, curr_fri, max_fri, vtype);
		printf("m_bui_threshold, m_ffmc_threshold = %f, %f\n", m_bui_threshold, m_ffmc_threshold);
		assert(0);
	}

	// calc portion of cell burned          
	part = runParamsP->fire_set_interval>=0 ? 1. : (burn_count + 1.) / curr_fri; 
	if (part < 0.) part = 0.;
	else if (part > 1.0) part = 1.0;

	// When appropriate, apply the fire suppression rule.
	if (try2suppressFlag && part>runParamsP->suppressed_fire_cell_fraction)
	{ // Apply fire suppression rule.
		if ((m_fire_fli<runParamsP->fire_suppression_fli_threshold 
					&& m_fire_ros<runParamsP->fire_suppression_ros_threshold) 
				|| m_fire_erc<runParamsP->fire_suppression_erc_threshold) 
		{ // Simulate this fire as a suppressed fire.
			part = runParamsP->suppressed_fire_cell_fraction;
		}
	} // end of fire suppression attempt

	return(part);
} // end of part_burn()


void MC_FireModel::fire_return_interval(int vtype, float * min_friP, float * max_friP)
{
	float * vveg2mfriP;

	/* get veg type dependent parameters */
	switch (runParamsP->baseCalibration)
	{
		case mc2WCR:
		case mc2GLOBAL:
		case mc2ConUS:
		case mc2ConUS_LC:
		case mc2California:
			vveg2mfriP = (float *)&(modelParamsP->vveg2mfri[0][0]);
			break;
			/*      
				case mc2CA08: 
				case mc2YOSE:
				assert(0); // vveg2fuelP = vveg2mfriCA08; 
				break;

				case mc2VINCERA:
				assert(0); // vveg2loadP = vveg2mfriVINCERA;
				break;
				*/      
		default: assert(0); break;
	}

	*min_friP = *(vveg2mfriP + vtype*3 + 1);
	*max_friP = *(vveg2mfriP + vtype*3 + 2);
	if ((*min_friP)>(*max_friP) || (*max_friP)<=0.f || (*min_friP)<=0.f )
	{
		printf("*** fire_return_interval(): vtype, *min_friP, *max_friP = %d, %f, %f\n",
				vtype, *min_friP, *max_friP);
		assert(0);
	}

} // end of MC_FireModel::fire_return_interval()


void MC_FireModel::FireEffect(int fire_dom, int fire_doy, int vtype)
{
	/* KEY TO CONSUME AND KILLED VARS 

	   consumed[0] = % live grass consumption
	   consumed[1] = % live tree leaf consumption
	   consumed[2] = % live tree branch consumption
	   consumed[3] = % dead 1-hr consumption
	   consumed[4] = % dead 10-hr consumption
	   consumed[5] = % dead 100-hr consumption
	   consumed[6] = % dead 1000-hr consumption

	   killed[0]  = % live tree leaf killed
	   killed[1]  = % live tree branch killed
	   killed[2]  = % live tree stem killed
	   killed[3]  = % live root killed
	   */

	m_consumed[0] = 90.; /* 90% live grass consumed every fire event */

	// Is this a crown fire?
	bool crown_fire_event;
	{
		float mc_tree_frac = m_mc_tree_dom[fire_dom]*0.01f;

		float h = 460. + (26.*mc_tree_frac); // estimate heat of canopy ignition (kJ/kg) 

		// estimate critical fireline intensity for crown ignition (kW/m)
		float z = m_tree_ht_m - m_crown_len_m;  
		float fli_crit_British = pow((0.010 * z * h), 1.5); // BTU s-1 ft-1
		float fli_crit_metric = sciFn.btu_per_sec_per_ft_to_kW_per_m(fli_crit_British);

		crown_fire_event = m_fire_fli>fli_crit_metric;
	} // end of block for crown fire

	if (crown_fire_event)
	{
		/* if crown fire, all live compartments die, live 
		   leaves and branches are consumed, live large
		   wood and live roots transfered to dead compartment */ 
		m_consumed[1]  = modelParamsP->crown_fire_mortality_pct;
		m_consumed[2]  = modelParamsP->crown_fire_mortality_pct;

		m_killed[0] = 0.;
		m_killed[1] = 0.;
		m_killed[2] = modelParamsP->crown_fire_mortality_pct;
		m_killed[3] = modelParamsP->crown_fire_mortality_pct;
	} // end of crown fire block
	else
	{ // Not a crown fire. Is there complete mortality?
		bool mort;
		float wnd_mph = sciFn.m_per_sec_to_mph(m_wnd_doy[fire_doy]);
		float scorch_ht_ft = (63. / (158. - m_tmp_doy_F[fire_doy])) *
			(pow(m_fire_fli, (7. / 6.)) /
			 pow((m_fire_fli + pow(wnd_mph, 3.)), .5));
		float scorch_ht_m = sciFn.ft_to_m(scorch_ht_ft);
		float crown_kill = sciFn.crown_kill(m_tree_ht_m, m_crown_len_m, scorch_ht_m);
		float ck_pct = crown_kill*100.; // convert crown kill from fraction to percentage 
		float prob_mort, term1, term2, term3, term4;
		if (ck_pct<=0.0) prob_mort = 0.0;
		else
		{
			term1 = 1.910 * m_bark_thick;
			term2 = -0.1775 * pow(m_bark_thick, 2.);
			term3 = -0.000541 * pow(ck_pct, 2.);
			term4 = -1.466 + term1 + term2 + term3;
			prob_mort = (1. / (1. + exp(term4))) * 100.;
		}
		prob_mort = sciFn.clip(prob_mort, 0., 100.);
		mort = prob_mort > fireParams.prob_thres;

		if (mort)
		{
			/* Not crown fire but complete mortality, so all live
			   compartments are transfered to dead compartments */ 
			m_consumed[1]  = 0.;  
			m_consumed[2]  = 0.;  

			m_killed[0] = 100.;
			m_killed[1] = 100.;
			m_killed[2] = 100.;
			m_killed[3] = 100.;
		} // end of block for complete mortality without crown fire
		else
		{ /* live leaves, fine branches, and roots die
		     proportion to crown kill */ 
			m_consumed[1]  = 0.;   
			m_consumed[2]  = 0.;   

			float pct_killed = (crown_kill>modelParamsP->fire_min ? crown_kill : modelParamsP->fire_min)*100.;
			for (int i = 0; i < 4; i++) m_killed[i] = pct_killed; 
		} // end of block for fire without complete mortality
	} // end of blocks for calculating fire effects on live vegetation

	float c_lgrass   = sciFn.g_per_m2_to_tons_per_acre(m_lgras)  * (m_consumed[0] / 100.);
	float c_lleaf   = sciFn.g_per_m2_to_tons_per_acre(m_lleaf)   * (m_consumed[1] / 100.);
	float c_lbranch = sciFn.g_per_m2_to_tons_per_acre(m_lbranch) * (m_consumed[2] / 100.);

	// Calculate consumption of dead fuels, using CONSUME v.2.1 equations (Ottmar et al.).
	// These variables are calculated in the consumption block, used in the emissions block.
	float dia_reduc, c_1hr, c_10hr, c_100hr, c_1000hr, d100, d1000; 
	{
		/* 1-HR FUEL CONSUMPTION */
		/* Assume 90% of 1-hr dead fuels are consumed */ 
		float d1 = sciFn.g_per_m2_to_tons_per_acre(m_d1hr); 
		c_1hr = d1 * 0.9;
		m_consumed[3] = 90.0;    


		/* 10-HR FUEL CONSUMPTION */
		float d10 = sciFn.g_per_m2_to_tons_per_acre(m_d10hr); 
		c_10hr = -0.048132 + (0.917393 * d10); 
		if (c_10hr < 0.0) c_10hr = 0.0;
		m_consumed[4] = d10>0. ? (c_10hr/d10)*100. : 0.;


		/* 100-HR FUEL CONSUMPTION */ 
		d100 = sciFn.g_per_m2_to_tons_per_acre(m_d100hr); 
		c_100hr = -0.124649 + (0.869309 * d100) - (0.004804 * m_mc_duff);
		if (c_100hr<0.0) c_100hr = 0.0;
		m_consumed[5] = d100>0. ? (c_100hr/d100)*100. : 0.;


		/* 1000-HR FUEL CONSUMPTION */
		float preburn_dia, QMD, term5;
		d1000 = sciFn.g_per_m2_to_tons_per_acre(m_d1000hr); 
		preburn_dia = QMD = 6.6;
		if (m_mc_duff <= 70.) dia_reduc = -1.465442 + (0.466083 * preburn_dia) - (0.014756 * m_mc_duff);
		else
		{
			term5 = 0.03 * (m_mc_duff - 70.);
			dia_reduc = 0.5779 * exp(-term5);
		}  
		float vol_reduc = (1. - pow(((QMD - dia_reduc) / QMD), 2.));
		c_1000hr = vol_reduc * d1000;
		if (c_1000hr < 0.0) c_1000hr = 0.0;
		m_consumed[6] = d1000>0. ? (c_1000hr/d1000)*100. : 0.;


		/* LARGE DEAD WOOD CONSUMPTION */
		// cf_ldw = (c_10hr + c_100hr + c_1000hr) / (d10 + d100 + d1000);
	} // end of block to calculate consumption of dead fuels

	// emissions();
	{
		float fc_1hr, fc_10hr, fc_100hr, fc_1000hr;
		float fc_lgrass, fc_lleaf, fc_lbranch;
		float term1, term2, term3, term4;
		float f_portion, f_dia_reduc;
		float fcf_1000hr;
		float tot_c, tot_fc, tot_sc;
		float fem_co2, sem_co2, fem_co, sem_co, fem_ch4, sem_ch4, fem_nmhc, sem_nmhc;
		float fem_pm, sem_pm, fem_pm2p5, sem_pm2p5, fem_pm10, sem_pm10;

		/* assume all dead and live fine fuel consumption due to flaming combustion */
		fc_1hr = c_1hr;
		fc_10hr = c_10hr;
		fc_lgrass = c_lgrass;
		fc_lleaf = c_lleaf;
		fc_lbranch = c_lbranch;

		/* flaming large fuel diameter reduction */
		term1 = ((20.0 - c_100hr) / 20.0) - 1.0;
		term2 = term1 / 0.2313;
		if (term2<=0.0) f_portion = 0.0; /* ??? */
		else
		{
			term3 = pow(term2, -2.260);
			term4 = -1.0 * term3;
			f_portion = 1.0 - exp(term4);
		}
		f_dia_reduc = dia_reduc * f_portion;


		/* flaming 100-hr fuel combustion */
		if (f_dia_reduc > 1.68) fc_100hr = c_100hr;
		else
		{
			term1 = pow((1.68 - f_dia_reduc), 2.0);
			term2 = term1 / 2.8224;
			term3 = 1.0 - term2;
			fc_100hr = d100 * term3;
		}

		/* flaming 1000-hr combustion */
		term1 = pow((6.6 - f_dia_reduc), 2.0);
		term2 = term1 / 43.56;
		fcf_1000hr = 1.0 - term2;
		fc_1000hr = d1000 * fcf_1000hr;

		/* smoldering  = total - flaming combustion */
		tot_c  = c_1hr + c_10hr + c_100hr + c_1000hr;
		tot_c  = tot_c + c_lgrass + c_lleaf + c_lbranch;

		tot_fc = fc_1hr + fc_10hr + fc_100hr + fc_1000hr;
		tot_fc = tot_fc + fc_lgrass + fc_lleaf + fc_lbranch;

		tot_sc = tot_c - tot_fc;

		/* calculate emissions (lbs/acre) */
		fem_co2 = tot_fc * emfac(vtype, CO2em, FLAMING);
		sem_co2 = tot_sc * emfac(vtype, CO2em, SMOLDERING);
		float em_co2  = fem_co2 + sem_co2;

		fem_co = tot_fc * emfac(vtype, COem, FLAMING);
		sem_co = tot_sc * emfac(vtype, COem, SMOLDERING);
		float em_co  = fem_co + sem_co;

		fem_ch4 = tot_fc * emfac(vtype, CH4em, FLAMING);
		sem_ch4 = tot_sc * emfac(vtype, CH4em, SMOLDERING);
		float em_ch4  = fem_ch4 + sem_ch4;

		fem_nmhc = tot_fc * emfac(vtype, NMHCem, FLAMING);
		sem_nmhc = tot_sc * emfac(vtype, NMHCem, SMOLDERING);
		float em_nmhc  = fem_nmhc + sem_nmhc;

		fem_pm = tot_fc * emfac(vtype, PMem, FLAMING);
		sem_pm = tot_sc * emfac(vtype, PMem, SMOLDERING);
		float em_pm  = fem_pm + sem_pm;

		fem_pm2p5 = tot_fc * emfac(vtype, PM2P5em, FLAMING);
		sem_pm2p5 = tot_sc * emfac(vtype, PM2P5em, SMOLDERING);
		float em_pm2p5  = fem_pm2p5 + sem_pm2p5;

		fem_pm10 = tot_fc * emfac(vtype, PM10em, FLAMING);
		sem_pm10 = tot_sc * emfac(vtype, PM10em, SMOLDERING);
		float em_pm10  = fem_pm10 + sem_pm10;

		/* convert lbs/acre to gm/m2 */
		m_em_co2 = sciFn.lbs_per_ac_to_g_per_m2(em_co2);
		m_em_co = sciFn.lbs_per_ac_to_g_per_m2(em_co);
		m_em_ch4 = sciFn.lbs_per_ac_to_g_per_m2(em_ch4);
		m_em_nmhc = sciFn.lbs_per_ac_to_g_per_m2(em_nmhc);
		m_em_pm = sciFn.lbs_per_ac_to_g_per_m2(em_pm);
		m_em_pm2p5 = sciFn.lbs_per_ac_to_g_per_m2(em_pm2p5);
		m_em_pm10 = sciFn.lbs_per_ac_to_g_per_m2(em_pm10);

	} // end of emissions block

	// black_carbon
	{
		/* KEY TO BLACK CARBON VARS

		   blkc[0] = % live grass pool converted to black carbon
		   blkc[1] = % live tree leaf converted to black carbon
		   blkc[2] = % live tree branch converted to black carbon
		   blkc[3] = % dead 1-hr converted to black carbon
		   blkc[4] = % dead 10-hr converted to black carbon
		   blkc[5] = % dead 100-hr converted to black carbon
		   blkc[6] = % dead 1000-hr converted to black carbon
		   */

		for (int i = 0; i < 7; i++) if (m_consumed[i] == 0.) m_blkc[i] = 0.;
		else
		{
			float cterm;
			float term1, term2,term3,term4;

			/* Kuhlbusch et al 1995 GBC equation 3 */
			if (m_consumed[i] < 70.) cterm = 70.;  
			else cterm = m_consumed[i];

			term1 = 28.5;
			term2 = 88.2 - cterm;
			term3 = pow(1.3, term2) + 1.;
			term4 = (100. - cterm) / 100.;
			m_blkc[i] = ((term1 / term3) * term4)/100.f;
		}

	} // end of black carbon block

} // end of MC_FireModel::FireEffect()


float MC_FireModel::emfac(int vtype, EmissionSpeciesEnum em_species, EmissionSourceEnum em_source)
{
	EmissionFuelTypeEnum em_fuel_type;

	float emfacGLOBAL[NUM_EM_FUEL_TYPES][NUM_EM_SOURCES][NUM_EM_SPECIES] = {
		//typedef enum {PMem=0, PM10em, PM2P5em, COem, CO2em, CH4em, NMHCem} EmissionSpeciesEnum; 

		// em fuel type DOUGLAS_FIR_SLASH 
		{{24.7, 16.6, 14.9, 143., 3385., 4.6, 4.2}, // flaming
			{35.0, 27.6, 26.1, 463., 2804., 15.2, 8.4}}, // smoldering

		// em fuel type HARDWOODS_SLASH
		{{23.0, 14.0, 12.2, 92., 3389., 4.4, 5.2}, // flaming
			{38.0, 25.9, 23.4, 466., 2851., 19.6, 14.0}}, // smoldering

		// em fuel type PONDEROSA_LODGEPOLE_PINE_SLASH
		{{18.8, 11.5, 10.0, 89., 3401., 3.0, 3.6}, // flaming
			{48.6, 36.7, 34.2, 285., 2971., 14.6, 9.6}}, // smoldering

		// em fuel type MIXED_CONIFER_SLASH
		{{22.0, 11.7, 9.6, 53., 3458., 3.0, 3.2}, // flaming
			{33.6, 25.3, 23.6, 273., 3023., 17.6, 13.2}}, // smoldering

		// em fuel type JUNIPER_SLASH
		{{21.9, 15.3, 13.9, 82., 3401., 3.9, 5.5}, // flaming
			{35.1, 25.8, 23.8, 250., 3050., 20.5, 15.5}}, // smoldering

		// em fuel type SAGEBRUSH
		{{45.0, 31.8, 29.1, 155., 3197., 7.4, 6.8}, // flaming
			{45.3, 29.6, 26.4, 212., 3118., 12.4, 14.5}}, // smoldering

		// em fuel type CHAPARRAL
		{{31.6, 16.5, 13.5, 119., 3326., 3.4, 17.2}, // flaming
			{40.0, 24.7, 21.6, 197., 3144., 9.0, 30.6}}, // smoldering

		// em fuel type MEAN_VALUES
		{{26.7, 16.8, 14.7, 105., 3365., 4.2, 6.5}, // flaming
			{39.4, 30.0, 25.6, 292., 2994., 15.6, 15.1}}}; // smoldering

	EmissionFuelTypeEnum em_typeGLOBAL[] = {
		MEAN_VALUES, // 0 unused
		MEAN_VALUES, // 1 ice aka barren COLD_BARRENveg
		MIXED_CONIFER_SLASH, // 2 tundra aka alpine TUNDRAveg
		MIXED_CONIFER_SLASH, // 3 taiga-tundra TAIGA_TUNDRAveg
		MIXED_CONIFER_SLASH, // 4 boreal needleleaf forest BOREAL_NEEDLELEAF_FORESTveg
		MIXED_CONIFER_SLASH, // 5 boreal mixed woodland BOREAL_WOODLANDveg
		MIXED_CONIFER_SLASH, // 6 subalpine forest SUBALPINE_FORESTveg
		DOUGLAS_FIR_SLASH, // 7 maritime needleleaf forest MARITIME_EN_FORESTveg
		MIXED_CONIFER_SLASH, // 8 temperate needleleaf forest MESIC_TEMPERATE_NEEDLELEAF_FORESTveg
		HARDWOODS_SLASH, // 9 temperate deciduous broadleaf forest TEMPERATE_DB_FORESTveg
		MIXED_CONIFER_SLASH, // 10 temperate cool mixed forest COOL_MIXED_FORESTveg
		PONDEROSA_LODGEPOLE_PINE_SLASH, // 11 temperate warm mixed forest TEMPERATE_WARM_MIXED_FORESTveg
		PONDEROSA_LODGEPOLE_PINE_SLASH, // 12 temperate needleleaf woodland TEMPERATE_EN_WOODLANDveg
		HARDWOODS_SLASH, // 13 temperate deciduous broadleaf woodland TEMPERATE_DB_WOODLANDveg
		JUNIPER_SLASH, // 14 temperate cool mixed woodland TEMPERATE_COOL_MIXED_WOODLANDveg
		JUNIPER_SLASH, // 15 temperate warm mixed woodland TEMPERATE_WARM_MIXED_WOODLANDveg
		SAGEBRUSH, // 16 C3 shrubland C3SHRUBveg
		SAGEBRUSH, // 17 C3 grassland C3GRASSveg
		SAGEBRUSH, // 18 temperate desert TEMPERATE_DESERTveg
		PONDEROSA_LODGEPOLE_PINE_SLASH, // 19 subtropical needleleaf forest SUBTROPICAL_EN_FORESTveg
		PONDEROSA_LODGEPOLE_PINE_SLASH, // 20 subtropical deciduous broadleaf forest SUBTROPICAL_DB_FORESTveg
		HARDWOODS_SLASH, // 21 subtropical evergreen broadleaf forest WARM_EB_FORESTveg
		PONDEROSA_LODGEPOLE_PINE_SLASH, // 22 subtropical mixed forest SUBTROPICAL_MIXED_FORESTveg
		HARDWOODS_SLASH, // 23 subtropical needleleaf woodland SUBTROPICAL_EN_WOODLANDveg
		HARDWOODS_SLASH, // 24 subtropical deciduous broadleaf woodland SUBTROPICAL_DB_WOODLANDveg
		PONDEROSA_LODGEPOLE_PINE_SLASH, // 25 subtropical evergreen broadleaf woodland SUBTROPICAL_EB_WOODLANDveg
		HARDWOODS_SLASH, // 26 subtropical mixed woodland SUBTROPICAL_MIXED_WOODLANDveg
		SAGEBRUSH, // 27 C4 shrubland C4SHRUBveg
		SAGEBRUSH, // 28 C4 grassland C4GRASSveg
		SAGEBRUSH, // 29 subtropical desert SUBTROPICAL_DESERTveg
		HARDWOODS_SLASH, // 30 tropical evergreen broadleaf forest TROPICAL_EB_FORESTveg
		HARDWOODS_SLASH, // 31 tropical deciduous woodland TROPICAL_DECIDUOUS_WOODLANDveg
		HARDWOODS_SLASH, // 32 tropical savanna TROPICAL_SAVANNAveg
		DOUGLAS_FIR_SLASH, // 33 WCR_SSZveg
		DOUGLAS_FIR_SLASH, // 34 WCR_PSFZveg
		SAGEBRUSH, // 35 tropical desert TROPICAL_DESERTveg
		MIXED_CONIFER_SLASH, // 36 cool moist needleleaf forest MOIST_TEMPERATE_NEEDLELEAF_FORESTveg
		MEAN_VALUES, // 37 unused or Lynx_AgricultureGrazing
		MIXED_CONIFER_SLASH, // 38 subalpine meadow SUBALPINE_MEADOWveg
		MEAN_VALUES, // 39 water and wetlands WATERveg
		MEAN_VALUES, // 40 natural barren NATURAL_BARRENveg
		MEAN_VALUES, // 41 developed DEVELOPEDveg
		MIXED_CONIFER_SLASH, // 42 larch forest 
		DOUGLAS_FIR_SLASH, // 43 SSZ, from 7 maritime needleleaf forest MARITIME_EN_FORESTveg
		DOUGLAS_FIR_SLASH, // 44 WHZ, from 7 maritime needleleaf forest MARITIME_EN_FORESTveg
		DOUGLAS_FIR_SLASH, // 45 PSFZ, from 7 maritime needleleaf forest MARITIME_EN_FORESTveg
		MIXED_CONIFER_SLASH, // 46 MHZ, from 6 subalpine forest SUBALPINE_FORESTveg
		MIXED_CONIFER_SLASH, // 47 SAFZ, from 6 subalpine forest SUBALPINE_FORESTveg
		MIXED_CONIFER_SLASH, // 48 PKLZ, from 5 boreal mixed woodland BOREAL_WOODLANDveg 
		MIXED_CONIFER_SLASH, // 49 cool dry needleleaf forest DRY_TEMPERATE_NEEDLELEAF_FORESTveg
		// PONDEROSA_LODGEPOLE_PINE_SLASH, // 50 XERIC_NEEDLELEAF_WOODLANDveg
	}; // end of em_typeGLOBAL[]
	switch (runParamsP->baseCalibration)
	{
		/*
		   case mc2VEMAP: 
		   case mc2NA8km:
		   case mc2WWETAC:
		   case mc2CA08: 
		   case mc2YOSE:
		   case mc2VINCERA:
		   assert(0);
		   break;
		   case mc2US50km:
		   */
		case mc2WCR:
		case mc2ConUS:
		case mc2GLOBAL:
		case mc2ConUS_LC:
		case mc2California:
			em_fuel_type = em_typeGLOBAL[vtype];
			return(emfacGLOBAL[em_fuel_type][em_source][em_species]);
		default: assert(0); break;
	}
} // end of MC_FireModel::emfac()


void MC_FireModel::ann_effect(float cen_outvars[600], float cen_state[200], float part_burn)
{
	float lgras, lleaf, lwod1, lwod2;
	float dstnd, mlittr, slittr, littr, dwod1, dwod2;
	float tmp[8];
	float bio_sum1, bio_sum2;

	/* make var assigments from CENTURY output */
	lgras   = cen_outvars[3];
	lleaf   = cen_outvars[445];
	lwod1   = cen_outvars[419];
	lwod2   = cen_outvars[452];
	dstnd   = cen_outvars[262];
	mlittr  = cen_outvars[104];
	slittr  = cen_outvars[251];
	littr   = mlittr + slittr;
	dwod1   = cen_outvars[482];
	dwod2   = cen_outvars[483];

	bio_sum1 = lgras + lleaf + lwod1 + lwod2 + dstnd + littr + dwod1 + dwod2;
	bio_sum2 = lleaf + lwod1 + lwod2;

	/* calculate biomass consumed */

	tmp[0] = lleaf *   cen_state[0];
	tmp[1] = lwod1 *   cen_state[1];
	tmp[2] = lwod2 *   cen_state[2];
	tmp[3] = dwod1 *   cen_state[3];
	tmp[4] = dwod2 *   cen_state[4];
	tmp[5] = lgras *   cen_state[5];
	tmp[6] = dstnd *   cen_state[6];
	tmp[7] = littr *   cen_state[7];

	m_consume_totbio = 0.;
	for (int i = 0; i < 8; i++) m_consume_totbio += tmp[i];

	m_consume_live = tmp[0] + tmp[1] + tmp[2] + tmp[5];
	m_consume_dead = tmp[3] + tmp[4] + tmp[6] + tmp[7]; 

	/* calculate black carbon formation */
	tmp[0] = lleaf *  m_blkc[1] * part_burn; 
	tmp[1] = lwod1 *  m_blkc[2] * part_burn;
	tmp[2] = lwod2 *  m_blkc[2] * part_burn; 
	tmp[3] = dwod1 *  m_blkc[3] * part_burn;
	tmp[4] = dwod2 *  m_blkc[5] * part_burn;
	tmp[5] = lgras *  m_blkc[0] * part_burn;
	tmp[6] = dstnd *  m_blkc[3] * part_burn;
	tmp[7] = littr *  m_blkc[3] * part_burn;

	m_blkc_totbio = 0.;
	for (int i = 0; i < 8; i++) m_blkc_totbio += tmp[i];

	/* calculate biomass killed */
	tmp[0] = lleaf * cen_state[8];
	tmp[1] = lwod1 * cen_state[9];
	tmp[2] = lwod2 * cen_state[10];

	m_death_totbio = 0.;
	for (int i = 0; i < 3; i++) m_death_totbio += tmp[i];

	/* assign value to stem death */
	// death_stem = tmp[2];

} // end of MC_FireModel::ann_effect()


