/* MAPSSfcns.cpp */

#ifndef UNIX_MAPSS
#include "StdAfx.h"
#endif

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#ifndef UNIX_MAPSS
  #include "EnvModel.h"
#endif
#include "ScienceFcns.h"
#include "category_bgc.h"
#include "MAPSSvegClasses.h"
#include "ProcessModel.h"
#include "MAPSSbiogeographyModel.h"

#ifndef ASSERT
   #define ASSERT assert
#endif

#ifdef UNIX_MAPSS
  extern Boolean diagsMAPSS;
#else
  extern bool diagsMAPSS;
  extern CString diagMsg;
#endif

/*
 *
 * ----------------------------  GrassPotAtSatisfied
 *
 */
bool MAPSS_BiogeogModel::GrassPotAtSatisfied(float pot_at, State * curr_mo_State, float conductance, float lai[], 
    float * grass_lai, float begin_h2o, float pot_transp)
	/*******
		check of curr_mo_State->transpire vs. pot_transp below is required to
		preclude calculation of negative lai values
			(if x >= 1, log x >= 0, *lai negative)
	********/
{
	float	begin_lai;
	bool	DoLai;
  
  begin_lai = *lai;
  DoLai = m_cats.mclass==Unknown;
	
	if (pot_at == 0.0) 
  {	// all demand was satisfied 
		if ((*lai == *grass_lai) || (curr_mo_State->soil_h2o[SURFACE] <= (m_inadeq_h2o[SURFACE] + m_soil.wilt_pt[SURFACE]))) 			
        return(TRUE);
		else 
    {
      curr_mo_State->transpire[GRASS] += curr_mo_State->soil_h2o[SURFACE] - m_soil.wilt_pt[SURFACE];
      if (curr_mo_State->transpire[GRASS] >= pot_transp) 
      { // revise lai upward and repeat calculations 
        curr_mo_State->soil_h2o[SURFACE] = begin_h2o;
        if (DoLai) 
        {
    			*lai = *grass_lai;
          return(FALSE);
        }
        else return(TRUE);
      }
    }
	}
	else if (begin_h2o == m_soil.wilt_pt[SURFACE]) // ??? should this be <= instead of ==?
  { // no water in upper soil layer for plants 
		curr_mo_State->soil_h2o[SURFACE] = begin_h2o;
		if (DoLai) 
    {
			*lai = 0.0;
			return(FALSE);
		}
		else return(TRUE);
	}
  
	if (DoLai) 
  {
    *lai = ((-cond_lai_max[m_zone][GRASS][m_cats.broadleaf]) *
				log((-(curr_mo_State->transpire[GRASS] / pot_transp)) + 1)) / (conductance * k_transp[m_zone][SURFACE][0][m_cats.broadleaf]);
		if (*lai > *grass_lai) *lai = *grass_lai; // cannot exceed maximum potential lai 
		if (fabs(begin_lai - *lai) < 0.01) 
    {	/* Not worth expend effort */
			*lai = begin_lai;			
			return(TRUE);
		}
	}
  
	curr_mo_State->soil_h2o[SURFACE] = begin_h2o; // preparation for next iteration, if any
	return(!DoLai);
} // end of GrassPotAtSatisfied()




/*
 *
 * ----------------------------  AccumulateTranspiration
 *
 */
float MAPSS_BiogeogModel::AccumulateTranspiration(State * curr_mo_State, float * pot_transp, float water, int soil)
{
	float	mm_transpired;
	
	mm_transpired = MIN(*pot_transp, water);
	*pot_transp -= mm_transpired;
	curr_mo_State->soil_h2o[soil] -= mm_transpired;
	
	return(mm_transpired);
} // end of AccumulateTranspiration()


/*
 *
 * ----------------------------  calc_swp
 *
 */
double Soil_MAPSS::calc_swp(double pct_soil_h2o, int layer)
{
	double psi;	// water potential, kPa					
	double psi_e;	// water potential (kPa) at air entry	
  double h2o_10; // s_h2o at 10kPa (m^3/m^3)				
  double pct_by_vol;

	pct_by_vol = pct_soil_h2o*h2o_sat[layer];

	psi = aa[layer] * pow(pct_by_vol, bb[layer]);
	if (psi <= 10.0) 
  { //	set air entry water potential (psi_e) 
		psi_e = 100.0 * (m + (n * h2o_sat[layer]));	// water potential (kPa) at air entry		
		h2o_10 = exp((2.302 - log(aa[layer])) / bb[layer]); //	set 10kPa soil water content (h2o_10) 
		psi = 10.0 - (pct_by_vol - h2o_10)*(10.0 - psi_e)/(h2o_sat[layer] - h2o_10);
	}

	psi /= 1000.0; // Convert from KPa to MPa.
	return(psi);
} // end of calc_swp()


/*
 *
 * ----------------------------  StomatalConductanceMC2
 *
 */
float MAPSS_BiogeogModel::StomatalConductanceMC2(State * site, float pet, int lifeform, const float conductance[][2], 
    int mo, bool shrubs, bool broadleaf)
// Returns stomatal conductance.
// Also updates site->swp[lifeform] (a double) and site->pct_soil_h2o[lifeform] (also a double).
{
	double pct_soil_h2o,	/* soil h2o as fraction of saturation */
					sum_h2o,
					swp,				/* soil water potential (MPa) */
					pet_coeff,			/* pet effects coefficient */
					cond,				/* calculation of conductance */
					inner_value;
	int				i,
					k,			/* index to previous conductance */
					lform_index =
						(lifeform ? (shrubs ? SHRUB : TREE) : GRASS);
						/* lifeform index into global conductance parameters */	
	/******
		lifeform is either GRASS (0) or WOODY (1).
		lform_index is one of GRASS (0), TREE (1), or SHRUB (2).
	*******/
  
  sum_h2o = site->soil_h2o[SURFACE];
  if (lifeform!=GRASS) sum_h2o += site->soil_h2o[INTERMEDIATE];
   ASSERT(lifeform==GRASS || lifeform==WOODY); // m_total_saturation is only initialized for GRASS and WOODY
	pct_soil_h2o = sum_h2o / m_total_saturation[lifeform];
  
  // Note that, for lifeform==WOODY, while pct_soil_h2o is based on the water in the top 2 soil layers,
  // the calculation of swp in the next statement treats pct_soil_h2o as if it were specific to the
  // second layer alone.  For lifeform==WOODY, the swp for the second layer is used as if it were the
  // average swp for all the roots in both layers.
  swp = lifeform==GRASS ? -m_soil.calc_swp(pct_soil_h2o, SURFACE) : -m_soil.calc_swp(pct_soil_h2o, INTERMEDIATE);
	site->swp[lifeform] = swp;
	site->pct_soil_h2o[lifeform] = pct_soil_h2o;

	pet_coeff = (a_slope[lform_index] * pet);
	if (pet_coeff < 0.0) pet_coeff = 0.0;
	pet_coeff *= swp;

	inner_value = pet_coeff*pet_coeff - b[m_zone][lform_index][broadleaf]*swp*swp + c[m_zone][lform_index][broadleaf];
	if (inner_value < 0.0) return(0.0);
	cond = 0.5 * (sqrt(inner_value) + pet_coeff);

	if (cond < ((double) cond_min[m_zone][lform_index][broadleaf])) 
  {
		for (i = JAN; i < mo; i++) 
    {
			if (conductance[i][lifeform] < cond_min[m_zone][lform_index][broadleaf]) 
      {
				for (k = i + 1; k < mo; k++) if (conductance[k][lifeform] >= cond_min[m_zone][lform_index][broadleaf]) 
            return(cond_min[m_zone][lform_index][broadleaf]);
			}
		}
		for (i = cond_min_months[lform_index]; i; i--) 
    {
			if ((k = mo - i) < JAN) k += MONTHS;	/* wrap around to end of array */
			if (conductance[k][lifeform] >= cond_min[m_zone][lform_index][broadleaf]) return((float) MAX(cond, 0.0));
		}
		return(cond_min[m_zone][lform_index][broadleaf]); // set conductance constraint 
	}
	
	return((float) MAX(cond, 0.0));
}



/*
 * 
 * ----------------------------  MonthMin
 * 
 */
/*
 * for each evaluation, the cumulative one month precipitation value is
 * compared to the minimum for both the start and end of the time period; if
 * the time period straddles the end of a month, the accumulation of that
 * month may be a local minimum as well and must be checked
 */
void MAPSS_BiogeogModel::MonthMin(float current, float next, float prior, float * min, Interval * period)
{
	float	precip_value;

	
	if (period->beg <= 0.5) {
		precip_value = current + (prior - current) * (0.5f - period->beg);
		if (precip_value < *min) {
			*min = precip_value;
		}
		if ((period->end >= 0.5) && (current < *min)) {
			*min = current;
		}
	}
	else {
		precip_value = current + (next - current) * (period->beg - 0.5f);
		if (precip_value < *min) {
			*min = precip_value;
		}
	}

	if (period->end <= 0.5) {
		precip_value = current + (prior - current) * (0.5f - period->end);
		if (precip_value < *min) {
			*min = precip_value;
		}
	}
	else {
		precip_value = current + (next - current) * (period->end - 0.5f);
		if (precip_value < *min) {
			*min = precip_value;
		}
	}
	
}


/*
 * 
 * ----------------------------  OneThreshold
 * 
 */
/*******
	I think I could do this as a macro.

# define ONE_THRESHOLD(_current, _next, _tmp) \
	((_tmp - _curr_tmp) / (_next_tmp - _curr_tmp))
********/
float MAPSS_BiogeogModel::OneThreshold(float curr_tmp, float next_tmp, float frost_line)
{
	
	return((frost_line - curr_tmp) / (next_tmp - curr_tmp));
}



/*
 * 
 * ----------------------------  TimeInterval
 * 
 */
void MAPSS_BiogeogModel::TimeInterval(float curr_tmp, float next_tmp, float frost_line, Interval * period)
{
	
	if (curr_tmp < frost_line) {
		if (next_tmp > frost_line) {
			/******
				(cur < min) && (next > min)
			*******/
			period->end = 1.0;
			period->beg = OneThreshold(curr_tmp, next_tmp, frost_line);
		}
		else {
			/*******
				(cur < min) && (next <= min)
			********/
			period->end = 0.0;	/* season does not occur in this interval */
		}
	}
	else if (next_tmp < frost_line) {
		/******
			(cur >= min) && (next < min)
		*******/
		period->beg = 0.0;
		period->end = OneThreshold(curr_tmp, next_tmp, frost_line);
	}
	else {
		/******
			(cur >= min) && (next >= min)
			The entire time period is in the growing season.
		*******/
		period->beg = 0.0;
		period->end = 1.0;
	}
	
}


/*
 *
 * ----------------------------  PS_CenturyNew
 *
 */
/*******************************************************
	PS biomass canopy-type function based on CENTURY model
*******************************************************/
void MAPSS_BiogeogModel::PS_CenturyNew(float lai[][2], PsPath * biomass)
{
   double x, prodd;
	float			meantmp[MONTHS],
					n_days,
					tmp_diff,
					tmp_delt,
					soil_mean[MONTHS],
					tmp[DAYS_PER_YEAR],
					tmp_sum,
					soil_tmp[DAYS_PER_YEAR],
					live_bio[MONTHS],
					prod,
					max_prod,
					pstype[MONTHS],
					bio_sum,
					weight[MONTHS],
					C3,
					C3_75,
					C3_C4,
					C4_75,
					C4,
					max_weight;
	int				i,
					j,
					k,
					month,
					dayy,
					julday = 15,
					flagx = 0,
					prev_mo1,
					prev_mo2,
					/* next_mo, */
					C4_consec,
					mixed_consec;
          
	
	/*******************************************************
		VARIABLE INITIALIZATION
	********************************************************/

	for (month = JAN; month <= DEC; month++) {
		biomass->C3[month] = 0.0;
		biomass->C3_75[month] = 0.0;
		biomass->C3_C4[month] = 0.0;
		biomass->C4_75[month] = 0.0;
		biomass->C4[month] = 0.0;
		meantmp[month] = m_tmp[month];
	}

	/*******************************************************
		ESTIMATE DAILY MEAN TEMP BY INTERPOLATION
	********************************************************/

	for (month = JAN; month <= DEC; month++) {
		n_days = ((sciFn.days_per_mo[(month + 1) % 12] / 2.0f) +
			(sciFn.days_per_mo[month] / 2.0f));

		n_days = ceil(n_days);
    
		tmp_diff = meantmp[(month + 1) % 12] - meantmp[month];
		tmp_delt = (n_days != 0) ? (tmp_diff / n_days) : 0.0f;
 
		for (dayy = 0; dayy < (int) n_days; dayy++) {
			if (julday == DAYS_PER_YEAR) {
				julday = 0;
				flagx = 1;
			}
			if ((julday == 15) && (flagx == 1)) {
				break;
			}
			tmp[julday]  = meantmp[month] + (dayy * tmp_delt);
			julday += 1;
		}
	}

	/*******************************************************
		ESTIMATE DAILY SOIL TEMP FROM RUNNING AVE OF AIR TEMP
	********************************************************/

	for (i = 0; i < DAYS_PER_YEAR; i++) {
		if (tmp[i] < 0.0) {
			tmp[i] = 0.0;
		}
	}

	for (i = 0; i < DAYS_PER_YEAR; i++) {
		tmp_sum = 0.;
 
		for (j = 13; j >= 0; j--) {
 
			if ((i - j) < 0.0) {
				k = 364 + (i - j);
			}
			else {
				k = i - j;
			}
			tmp_sum += tmp[k];
		}
 
		soil_tmp[i] = tmp_sum / 14.0f;
 
		if (soil_tmp[i] < 0.0) {
			soil_tmp[i] = 0.0;
		}
	}

	/*******************************************************
		ESTIMATE MONTHLY SOIL TEMP 
	********************************************************/

	for (month = JAN; month <= DEC; month++) {
		soil_mean[month] = (soil_tmp[sciFn.month_days[month][MONTH_BEGIN]] + 
			soil_tmp[sciFn.month_days[month][MONTH_END]]) / 2.0f;
	}

	/*******************************************************
		ESTIMATE MONTHLY LIVE GRASS BIOMASS
	********************************************************/

	for (month = JAN; month <= DEC; month++) {
		if (lai[month][GRASS] <= 0.0) {
			live_bio[month] = 0.0;
		}
		else {
			live_bio[month] = (lai[month][GRASS] / 2.3f) / 0.0125f;
		}
	}

	/* beginning of monthly loop */
	for (month = JAN; month <= DEC; month++) {
		x = soil_mean[month];

		/******************************************************** 
			100% C3 DOMINANT 
		********************************************************/

		/* compute relative production of 100% C3 canopy
		using CENTURY model function, constrain value */

		prodd = (0.1788196117882357 + x * (0.03606057742960269+
			x * -0.001302093024389319)) /
			(1.0 + x * (-0.08856743168395311 +
			x * (0.004250159047418467 +
			x * -6.002654039746731E-05)));
      prod = (float)prodd;
		if (prod < 0.0) {
			prod = 0.0;
		}

		if (prod > 1.0) {
			prod = 1.0;
		}
	
		biomass->C3[month] = prod;

		/********************************************************
			75% C3 DOMINANT
		********************************************************/
 
		/* compute relative production of 75% C3 canopy
		using CENTURY model function, constrain value */
 
		prodd = (0.1079635295632225 + x *(0.02432464534358727 +
			x * -0.0007847687532581031)) /
			(1.0 + x * (-0.09067553364319364 +
			x * (0.003606799891237951 +
			x * -4.272256928231842E-05)));
      prod = (float) prodd;
		if (prod < 0.0) {
			prod = 0.0;
		}

		if (prod > 1.0) {
			prod = 1.0;
		}

		biomass->C3_75[month] = prod;

		/********************************************************
			C3-C4 CODOMINANT
		********************************************************/
 
		/* compute relative production of C3-C4 canopy
		using CENTURY model function, constrain value */

		prodd = (0.009993854887257886 + x * (0.008621196063664201 +
			x * (0.001632115117639580 + x * -4.768732294661518E-05))) /
			(1.0 + x * (-0.1011702820305769 +
			x * (0.005432818669449008 +
			x * -8.638320046043276E-05)));
      prod = (float) prodd;

		if (prod < 0.0) {
			prod = 0.0;
		}

		if (prod > 1.0) {
			prod = 1.0;
		}

		biomass->C3_C4[month] = prod;

		/********************************************************
			75% C4 DOMINANT
		********************************************************/
 
		/* compute relative production of 75% C4 canopy
		using CENTURY model function, constrain value */

		prodd = (0.02193052822277199 + x * (0.001399110037600328 +
			x * (0.002205857298167872 + x * -4.986575419838470E-05))) /
			(1.0 + x *(-0.0004880525815289832 +
			x * (-0.002485621630475747 +
			x * (0.0001037711905193199 +
			x * -1.009336575618056E-06))));
      prod = (float) prodd;
 
		if (prod < 0.0) {
			prod = 0.0;
		}
 
		if (prod > 1.0) {
			prod = 1.0;
		}

		biomass->C4_75[month] = prod;

		/********************************************************
			100% C4 DOMINANT
		********************************************************/

		/* compute relative production of 100% C4 canopy 
		using CENTURY model function, constrain value */

		prodd = (0.006239519724182143 + x * (0.006513678347464982 +
			x * (-0.0005775866974182693 + x * (0.0001491401308644086 +
			x * -3.094445566664409E-06)))) /
			(1.0 + x * (-0.0001694826018411467 +
			x * (-1.244515882217747E-05 +
			x * (-1.454616827810395E-05 +
			x * 7.533664205038196E-07))));
      prod = (float) prodd;

		if (prod < 0.0) {
			prod = 0.0;
		}
 
		if (prod > 1.0) {
			prod = 1.0;
		}

		biomass->C4[month] = prod;

		/********************************************************
			MONTHLY CANOPY TYPE DETERMINATION
		********************************************************/

		biomass->ps_type[month] = -1;
		max_prod = -1.0;

		if (live_bio[month] <= 0.) {
			biomass->ps_type[month] = 0;
		}
		else {
			if (biomass->C3[month] > max_prod) { 
				max_prod = biomass->C3[month]; 
				biomass->ps_type[month] = 1;
			}

			if (biomass->C3_75[month] > max_prod) {
				max_prod = biomass->C3_75[month]; 
				biomass->ps_type[month] = 2;
			}

			if (biomass->C3_C4[month] > (max_prod - 0.03)) { 
				max_prod = biomass->C3_C4[month]; 
				biomass->ps_type[month] = 3;
			}

			if (biomass->C4_75[month] > max_prod) { 
				max_prod = biomass->C4_75[month]; 
				biomass->ps_type[month] = 4; 
			}

			if (biomass->C4[month] > max_prod) {  
				max_prod = biomass->C4[month];  
				biomass->ps_type[month] = 5;
			}
		}

	} /* end of monthly loop */

	/********************************************************
		AGGREGATE CANOPY TYPE DETERMINATION
	********************************************************/
	bio_sum = 0.0;

	for (month = JAN; month <= DEC; month++) {
		if (biomass->ps_type[month] == 0) {
			pstype[month] = 0;
		}
		else if (biomass->ps_type[month] == 1) {
			pstype[month] = 1;
		}
		else if (biomass->ps_type[month] == 2) {
			pstype[month] = 2;
		}
		else if (biomass->ps_type[month] == 3) {
			pstype[month] = 3;
		}
		else if (biomass->ps_type[month] == 4) {
			pstype[month] = 4;
		}
		else if (biomass->ps_type[month] == 5) {
			pstype[month] = 5;
		}

		bio_sum += live_bio[month];
	}

	for (month = JAN; month <= DEC; month++) {
		weight[month] = (bio_sum != 0.0) ?
			(live_bio[month] / bio_sum) : 0.0f;
	}

	C3 =
		C3_75 =
		C3_C4 =
		C4_75 =
		C4 = 0.0;
	C4_consec = mixed_consec = 0;

	for (month = JAN; month <= DEC; month++) {

		if (month == JAN) {
			prev_mo1 = DEC;
		}
		else {
			prev_mo1 = month - 1;
		}

		if (month == JAN) {
			prev_mo2 = NOV;
		}
		else if (month == FEB) {
			prev_mo2 = DEC;
		}
		else {
			prev_mo2 = month - 2;
		}
		/***
		if (month == DEC) {
			next_mo = JAN;
		}
		else {
			next_mo = month + 1;
		}
		****/
		if (pstype[month] == 3 && soil_mean[month] >= 22.0) {
			pstype[month] = 4;
		}
 
		if ((pstype[month] >= 4) &&
				(pstype[prev_mo1] >= 4) &&
				(pstype[prev_mo2] >= 3)) {
			C4_consec += 1;
		}
		else if ((pstype[month] >= 3) &&
				(pstype[prev_mo1] >= 2)) {
			mixed_consec += 1;
		}

		if (pstype[month] == 1) {
			C3 += weight[month];
		}
		else if (pstype[month] == 2) {
			C3_75 += weight[month];
		}
		else if (pstype[month] == 3) {
			C3_C4 += weight[month];
		}
		else if (pstype[month] == 4) {
			C4_75 += weight[month];
		}
		else if (pstype[month] == 5) {
			C4 += weight[month];
		}
	}

	if (C4_consec > 0) {
		C4_consec += 1;
	}
	if (mixed_consec > 0) {
		mixed_consec += 1;
	}
	/****
	biomass->canopy = -1;
	****/  
	if (bio_sum <= 0.0) {
		/*****
			This could present a problem!
		******/
		biomass->canopy = C3C4Ignore;
	}
	else if (C4_consec >= 2) {
		biomass->canopy = C4Dominance;
	}
	else if (mixed_consec >=2) {
		biomass->canopy = C3C4Mixed;
	}
	else {
		max_weight = -1.0;

		if (C3 > max_weight) {
			max_weight = C3;
			biomass->canopy = C3Dominance;
		}
 
		if (C3_75 > max_weight) {
			max_weight = C3_75;
			biomass->canopy = C3C4Mixed;
		}

		if (C3_C4 > max_weight) {
			max_weight = C3_C4;
			biomass->canopy = C3C4Mixed;
		}

		if (C4_75 > max_weight) {
			max_weight = C4_75;
			biomass->canopy = C4Dominance;
		}

		if (C4 > max_weight) {
			max_weight = C4;
			biomass->canopy = C4Dominance;
		}
	}
	
} // end of PS_CenturyNew()


/*
 *
 * ----------------------------  GrassAlone
 *
 */
void MAPSS_BiogeogModel::GrassAlone(State site[], float conductance[][2])
{
	State			*curr_mo_State, *prev_mo;
	float			grass_lai[MONTHS], *currlai;
	unsigned		j;
	int	months, mo;
	State last_year[MONTHS];
	float threshold;
  float prev_grass_lai;	/* grass lai this month, prior cycle */
	int cycles;
  bool equilibrium;
  int H2Oiter;

	for (j = JAN; j < MONTHS; j++) 
  {
		grass_lai[j] = 0.0;
		m_lai[j][GRASS] = 0.0;
		m_lai[j][WOODY] = 0.0;
    last_year[j] = site[j];
	}

	LightAttenuate(0.0, grass_lai, m_lai, TRUE, site);
  
  // Set PET.  Setting m_pet_adj to tallgrass_pet_factor seems plausible, but setting m_petP to the SHRUB array rather than the GRASS
  // array seems odd.  But that's the way it is in the original MC1 code.
  m_pet_adj = tallgrass_pet_factor;
  m_petP = m_pet_array[SHRUB];

	InitSoilWater(site);
	InitSnow(site);

  cycles = 0;
  prev_grass_lai = -1.;
  threshold = MAX_DIFF_EQUILIB;
  months = 0;
  H2Oiter = 0;
	do 
  {	/* h2o equilibrium cycle */		
		mo = SetCurrPrev(months, site, &curr_mo_State, &prev_mo, &currlai);
		DistributePpt(curr_mo_State, currlai, mo);
		for (j = SURFACE; j <= DEEP; j++) curr_mo_State->soil_h2o[j] = prev_mo->soil_h2o[j];
		transpiration(curr_mo_State, conductance, currlai, grass_lai, mo);
    H2Oiter++;
    if (0 && diagsMAPSS) 
        printf("MC2 GrassAlone H2Oiter, mo, grass_lai[mo], conductance[mo][0], curr_mo_State->soil_h2o[0] = %d, %d, %f, %f, %f\n",
        H2Oiter, mo, grass_lai[mo], conductance[mo][0], curr_mo_State->soil_h2o[0]);
    { // Test whether equilibrium has been achieved.
      equilibrium = FALSE;
      if (months == MONTHS) 
      {
        cycles = 0;
        threshold = MAX_DIFF_EQUILIB;
      }
      else if (cycles == MAX_CYCLES) 
      { /* adjust equilibrium thresholds for excessive cycles */
        cycles = 0;
        threshold = threshold + MAX_DIFF_INCREMENT;
      }
      if (months >= MONTHS) 
      {	/* define basis for acceptable equilibrium */
        if (CheckStateValues(site, last_year, threshold, mo, 0.0) && (prev_grass_lai == m_lai[mo][0])) 
        { // all fields of site and last_year are equal and same grass lai
          if (threshold > MAX_DIFF_EQUILIB) fprintf(stderr, "GrassAloneEquilibrium: equil fact (unfettered) %5.3f\n",
              threshold);			
          equilibrium = TRUE;
        }
        else ++cycles;	/* unsuccessful equilibrium cycle */
      }
      
      if (!equilibrium)
      {	/* update stored state variables for next annual comparison */
        last_year[mo] = site[mo];
        prev_grass_lai = (mo < DEC) ? m_lai[mo + 1][GRASS] : m_lai[0][GRASS];
        months++;		/* number of months processed */
      }      
    } // end of block to test equilibrium
	} while (!equilibrium);   
} // end of GrassAlone()


/*
 * 
 * ----------------------------  initialize_lai
 * Initialize the LAI for each month for both trees and grass. 
 * Initialize m_offset to initial growing season tree LAI.
 * Return max tree LAI.
 */
float MAPSS_BiogeogModel::initialize_lai()
{
	int				mo;

	for (mo = JAN; mo < MONTHS; mo++) 
  {
		if (m_cats.deciduous && m_cats.broadleaf) 
    {
			if (m_tmp[mo] >= frost) m_lai[mo][WOODY] = m_LaiUpperBoundEnergy[TREE];
			else m_lai[mo][WOODY] = INIT_WOODY_LAI_MIN;
		}
		else m_lai[mo][WOODY] = m_LaiUpperBoundEnergy[TREE];
		m_lai[mo][GRASS] = m_grass_lai[mo] = INIT_GRASS_LAI_MIN;
	}
	m_offset = m_LaiUpperBoundEnergy[TREE];
	return(m_LaiUpperBoundEnergy[TREE]);
}



/*
 *
 * ----------------------------  ForestRules
 *
 */
MAPSSvegClass MAPSS_BiogeogModel::ForestRules(float lai_values[], int m_zone)
{
	MAPSSvegClass		classification = MisClass;

	if (m_cats.deciduous) 
   { // deciduous
		if (m_cats.broadleaf) 
      {
			if (m_zone == mapssSUBTROPICAL) classification = ForestMixedWarm_DEB;
			else if (m_zone == mapssTEMPERATE) 
         {
				if (m_cats.mixed) classification = ForestMixedCool;
				else 
            {
					if (lai_values[TREE] > north_hard_threshold) classification = ForestHardwoodCool;
					else classification = ForestDeciduousBroadleaf;
				}
			}
		}
		else ASSERT(0); // ForestDeciduousNeedle already pulled out before ForestRules() is called.
	}
	else 
   { // evergreen
		if (m_cats.broadleaf) 
      {
			if (m_zone == mapssTROPICAL) classification = ForestEvergreenBroadleafTropical;
		}
		else switch (m_zone) 
      { // evergreen needleleaf
			case mapssBOREAL: classification = ForestEvergreenNeedleTaiga;
				break;
			case mapssTEMPERATE:
				if (m_cats.mixed) classification = ForestMixedCool;
				else 
            {
               if (m_growing_deg_days_zero<=mapss_subalpine_threshold) classification = ForestSubalpine;
					else if ((m_max_tmp - m_min_tmp) > maritime_boundary[m_zone]) classification = ForestEvergreenNeedleContinental;
					else classification = ForestEvergreenNeedleMaritime;
				}
				break;
			case mapssSUBTROPICAL: classification = ForestMixedWarm_EN;
				break;
         default: ASSERT(0);
            break;
		}
	}

  ASSERT(classification!=MisClass);
	return(classification);
} // end of ForestRules()


/*
 *
 * ----------------------------  TreeSavannaRules
 *
 */
MAPSSvegClass MAPSS_BiogeogModel::TreeSavannaRules(float lai_values[], int m_zone)
{
	MAPSSvegClass		classification = MisClass;
	
	if (m_cats.deciduous) 
  {
		if (m_cats.broadleaf) 
    {
			if (m_zone == mapssTEMPERATE) classification = TreeSavannaDeciduousBroadleaf;
			else if (m_zone == mapssSUBTROPICAL) classification = TreeSavannaMixedWarm_DEB;
		}
		else ASSERT(0); // This should never occur. classification = TreeSavannaDeciduousNeedle;
	}
	else 
  { // evergreen
		if (m_cats.broadleaf) 
    {
			if (lai_values[TREE] > dry_trop_threshold) classification = ForestSeasonalTropical_ED;
			else classification = ForestSavannaDryTropical_ED;
		}
		else switch (m_zone) 
    {
			case mapssBOREAL: classification = TreeSavannaMixedCool_EN;
				break;
			case mapssTEMPERATE :
				if (m_growing_deg_days_zero<=mapss_subalpine_threshold) classification = TreeSavannaSubalpine;
				else if ((m_max_tmp - m_min_tmp) > maritime_boundary[m_zone]) 
        {
					if (lai_values[TREE] < pj_max_lai_continental) 
          {
						if (m_events_summer_avg <= pj_xeric_threshold) classification = TreeSavannaPJXericContinental;
						else classification = TreeSavannaPJContinental;
					}
					else classification = TreeSavannaEvergreenNeedleContinental;
				}
				else 
        { // maritime
					if (lai_values[TREE] < pj_max_lai_maritime) classification = TreeSavannaPJMaritime;
					else classification = TreeSavannaEvergreenNeedleMaritime;
				}
				break;
			case mapssSUBTROPICAL: classification = TreeSavannaMixedWarm_EN;
        break;
      default: ASSERT(0); break;
		}
	}

  ASSERT(classification!=MisClass);
	return(classification);
} // end of TreeSavannaRules()


/*
 *
 * ----------------------------  ShrubRules
 *
 */
MAPSSvegClass MAPSS_BiogeogModel::ShrubRules(float lai_values[], int m_zone)
{
	MAPSSvegClass		classification = MisClass;
	
	if (lai_values[GRASS] <= 0.0) classification = OpenShrublandNoGrass;
	else if (m_cats.deciduous) 
  {
		if (m_cats.broadleaf) 
    {
			if (m_zone == mapssTEMPERATE) classification = ShrubSavannaDeciduousBroadleaf;
			else if (m_zone == mapssSUBTROPICAL) classification = ShrubSavannaMixedWarm_DEB;
		}
		else ASSERT(0); // This should never occur.	classification = ShrubSavannaDeciduousMicro;
	}
	else 
  { // evergreen
		if (m_cats.broadleaf) classification = ShrubSavannaTropical_EB;
		else switch (m_zone) 
    {
			case mapssBOREAL: classification = ShrubSavannaMixedCool_EN;
				break;
			case mapssTEMPERATE: classification = ShrubSavannaEvergreenMicro;
				break;
			case mapssSUBTROPICAL :
				if (lai_values[GRASS] > xeric_savanna_threshold) classification = ShrubSavannaSubTropicalMixed;
				else 
        { 
					if (lai_values[SHRUB] > mediterranean_savanna_threshold) classification = ShrublandSubTropicalMediterranean;
					else classification = ShrublandSubTropicalXeromorphic;
				}
				break;
		}
	}

  ASSERT(classification!=MisClass);
	return(classification);
} // end of ShrubRules()


/*
 *
 * ----------------------------  GrassRules
 *
 */
MAPSSvegClass MAPSS_BiogeogModel::GrassRules(float lai_values[], int m_zone, C3C4Dominance canopy_type)
{

	MAPSSvegClass		classification = MisClass;
	
	if (lai_values[GRASS] > semi_desert_threshold) 
  {
		if (lai_values[GRASS] > short_grass_threshold) 
    {
			if (lai_values[GRASS] > tall_grass_threshold) switch (canopy_type)
      { // tall_grass_threshold<=lai_values
        case C3Dominance :
          classification = GrassTallC3;
          break;
        case C3C4Mixed :
          classification = GrassTallC3C4;
          break;
        case C4Dominance :
          classification = GrassTallC4;
          break;
        case C3C4Ignore :
        default : // ASSERT(0);
          break;
			}
			else switch (canopy_type) 
      { // short_grass_threshold<lai_values<=tall_grass_threshold						
        case C3Dominance :
          classification = GrassMidC3;
          break;
        case C3C4Mixed :
          classification = GrassMidC3C4;
          break;
        case C4Dominance :
          classification = GrassMidC4;
          break;
        case C3C4Ignore :
        default : // ASSERT(0);
          break;
			}
		}
		else switch (canopy_type) 
    { // semi_desert_threshold<lai_values<=short_grass_threshold
      case C3Dominance: classification = GrassShortC3;
        break;
      case C3C4Mixed: classification = GrassShortC3C4;
        break;
      case C4Dominance: classification = GrassShortC4;
        break;
      case C3C4Ignore :
      default : // ASSERT(0);
        break;
		}
	} 
	else 
  { // lai_values<=semi_desert_threshold
		if (lai_values[GRASS_SUM] > desert_grass_sum_threshold) switch (canopy_type) 
    {
      case C3Dominance: classification = GrassSemiDesertC3;
        break;
      case C3C4Mixed: classification = GrassSemiDesertC3C4;
        break;
      case C4Dominance: classification = GrassSemiDesertC4;
        break;
      case C3C4Ignore :
      default : // ASSERT(0);
        break;
		}
		else classification = DesertRules(m_zone);
	}
	
  // ASSERT(classification!=MisClass);
	return(classification);
} // end of GrassRules()


/*
 *
 * ----------------------------  DesertRules
 *
 */
MAPSSvegClass MAPSS_BiogeogModel::DesertRules(int m_zone)
{
	MAPSSvegClass		classification = MisClass;
	
	switch (m_zone) 
  {
		case mapssBOREAL: classification = DesertBoreal;
			break;
		case mapssTEMPERATE: classification = DesertTemperate;
			break;
		case mapssSUBTROPICAL: classification = DesertSubtropical;
			break;
		case mapssTROPICAL: classification = DesertTropical;
			break;
    default: ASSERT(0);
      break;
	}

  ASSERT(classification!=MisClass);	
	return(classification);
} // end of DesertRules()


/*
 *
 * ----------------------------  HeatLimitedRules
 *
 */
MAPSSvegClass MAPSS_BiogeogModel::HeatLimitedRules(int m_zone)
{
	MAPSSvegClass classification = MisClass;

	if (m_growing_deg_days_zero <= ice_boundary) classification = Ice;
	else if (m_growing_deg_days_zero <= tundra_boundary[m_zone]) classification = Tundra;
	else if (m_growing_deg_days_zero <= taiga_tundra_boundary[m_zone]) classification = TaigaTundra;

  ASSERT(classification!=MisClass);
	return(classification);
} // end of HeatLimitedRules()


/*
 *
 * ----------------------------  ClassifyStation
 *
 */
MAPSSvegClass MAPSS_BiogeogModel::ClassifyStation(float lai_values[], C3C4Dominance canopy_type)
{
	MAPSSvegClass		classification = MisClass;
	float			tsg;
	
	tsg = lai_values[TREE] + lai_values[SHRUB] + lai_values[GRASS];
  if (m_cats.deciduous == TRUE && m_cats.broadleaf == FALSE)
		classification = ForestDeciduousNeedle;
	else 
  {
		if (m_growing_deg_days_zero <= taiga_tundra_boundary[m_zone]) classification = HeatLimitedRules(m_zone);
		else if (lai_values[TREE] >= forest_threshold) classification = ForestRules(lai_values, m_zone);
		else if (lai_values[TREE] >= min_tree_lai[m_zone][m_cats.broadleaf]) classification = TreeSavannaRules(lai_values, m_zone);
		else if ((m_zone == mapssSUBTROPICAL) || (m_zone == mapssTROPICAL)) 
    {
			if ((m_max_grass >= max_grass_threshold) && (lai_values[GRASS_SUM] < cool_grass_threshold)) 
      {
				if (lai_values[SHRUB] > max_grass_shrub_threshold) classification = ShrubRules(lai_values, m_zone);
				else if (lai_values[GRASS] > desert_grass_threshold) classification = GrassRules(lai_values, m_zone, canopy_type);
				else classification = DesertRules(m_zone);
			}
			else if ((tsg >= tsg_threshold) && (lai_values[SHRUB] > desert_shrub_threshold)) classification = ShrubRules(
          lai_values, m_zone);
			else if (lai_values[GRASS] > desert_grass_threshold) classification = GrassRules(lai_values, m_zone, canopy_type);
			else classification = DesertRules(m_zone);
		}
		else if ((lai_values[SHRUB] > desert_shrub_threshold) && (tsg >= tsg_threshold)) classification = ShrubRules(
        lai_values, m_zone);
		else if (lai_values[GRASS] > desert_grass_threshold) classification = GrassRules(lai_values, m_zone, canopy_type);
		else
			classification = DesertRules(m_zone);
	}

  // ASSERT(classification!=MisClass);
	return(classification);
}




/*
 *
 * ----------------------------  CheckEvergreen
 *
 */
/*******
	This section of code is where things look into the areas where we
	(Ron and Jim really) think that there should be evergreen forest.
	The model is dredicting a deciduous forest.	 We have based a
	change on the growing season aet.  If there is not enough aet to
	support deciduous, based upon a threshold value, then switch to an
	evergreen and redo the calculation.
********/
bool MAPSS_BiogeogModel::CheckEvergreen(State site[], State site_begin[], float offset_begin, float lai_values[], 
    State ** curr_yearP, State ** prev_yearP)
// returns TRUE when switching to evergreen
{
	float			offset;
	Rulecat			new_cats;
	MAPSSvegClass		preclass;
	int				veg;
	int			mo;
	float		aet, prod_lai, ptemp;
   float f1p5 = 1.5;

  ASSERT(m_cats.deciduous);
	new_cats.broadleaf = FALSE;
	new_cats.deciduous = FALSE;
	new_cats.mixed = FALSE;

  // Calculate LAI from AET.
	aet = 0.0;
	for (mo = JAN; mo <= DEC; mo++) if (m_tmp[mo] >= frost) aet += site[mo].transpire[WOODY];
	ptemp = pow(m_maxlai, f1p5);	
   prod_lai = aet / ptemp;

	BroadDecid(&new_cats, prod_lai, (m_pet_adj == tree_pet_factor), TRUE);
  if (new_cats.deciduous) return(FALSE);

  m_cats.mixed = TRUE;

  /**********
    This is rather ugly.  I am going to have to do this again
    (for all veg) in the grasses routine.  I have to do it
    here in order to preclassify.  So, for now, it is done
    twice for those sites which "fail" the AET/Evergreen rule.
    It is really only a waste for those sites which do not
    have the water balance rerun (all the northern hardwood).
  ***********/
  if ((m_pet_adj == tree_pet_factor)) 
  {
    FindMeanLai(&lai_values[TREE_SUM], WOODY);
    lai_values[TREE] = m_maxlai;
  }
  else 
  {
    FindMeanLai(&lai_values[SHRUB_SUM], WOODY);
    lai_values[SHRUB] = m_maxlai;
  }

  preclass = ClassifyStation(lai_values, C3C4Ignore);

  /***********
    This is a very specific test.  The preclassification must
    be done to decide if the entire water balance must be
    rerun.
  ************/
  if (preclass != ForestMixedWarm_DEB) return(FALSE);
  
  fprintf(stderr, "Rerunning site for new (evergreen) water balance\n");
  m_cats.grassfire = FALSE;
  m_cats.deciduous = new_cats.deciduous;
  m_cats.broadleaf = new_cats.broadleaf;
  ASSERT(m_cats.mclass==Unknown);
  m_petP = m_pet_array[TREE]; // Points to the 12 monthly PET values for TREEs
  for (veg = GRASS; veg <= SHRUB; veg++) m_LaiUpperBoundEnergy[veg] = LaiUpperBoundsEnergy[m_zone][veg];
  m_maxlai = initialize_lai();

  LaiCycle(&offset, site_begin, offset_begin, curr_yearP, prev_yearP);
  GrassWoodyPart1(*curr_yearP, *prev_yearP, site, m_offset);
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
      GrassAlone(*prev_yearP, m_conductance);
      m_lai_values[GRASS] = FindMeanLai(&m_lai_values[GRASS_SUM], GRASS);
      m_lai_values[TREE] = m_lai_values[SHRUB] = m_lai_values[TREE_SUM] = m_lai_values[SHRUB_SUM] = 0.0;
    }
    else 
    {
      if ((m_pet_adj == tree_pet_factor)) 
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

  } // end of block for GrassWoodyPart2 

  return(TRUE);
	
} // end of CheckEvergreen()




/*
 *
 * ----------------------------  FindMeanLai
 *
 */
/******
	Return the average LAI for the growing season.	Modify a parm to
	hold the total LAI (LAD).
*******/
float MAPSS_BiogeogModel::FindMeanLai(float * total_lai, int lifeform)
{
	int	mon, growth_months;
  float lai_sum, growing_lai;

  lai_sum = 0.0;
  growth_months = 0;
	for (mon = JAN; mon <= DEC; mon++) 
  {
		lai_sum += m_lai[mon][lifeform];
		if (m_lai[mon][lifeform] > 0.0) ++growth_months;
	}

	*total_lai = lai_sum;
  growing_lai = growth_months>0 ? lai_sum/growth_months : 0.f;
	
	return(growing_lai);
}


/*
 *
 * ----------------------------  NoDeficits
 *
 */
bool MAPSS_BiogeogModel::NoDeficits(float deficit[])
{
	int	mo;

	for (mo = JAN; mo < MONTHS; mo++) if (deficit[mo] > 0.0) return(FALSE);
	
	return(TRUE);
} // end of NoDeficits()




/*
 *
 * ----------------------------  AdjustMaxlai
 *
 */
void MAPSS_BiogeogModel::AdjustMaxlai(float excess, State * site, float lai[][2],  
    float grass_lai[], bool deciduous, State curr_mo_State[])
{
	int				i;
	float			total_transpiration,
					prev_woody_lai,
					change,
					mininum_tree;
	
	mininum_tree = min_tree_lai[m_zone][m_cats.broadleaf];
  
	/* reduce tree lai relative to extent of moisture deficit */
	if (excess < 0.0) 
  {
		total_transpiration = -excess;
		for (i = JAN; i < MONTHS; i++, site++) total_transpiration += site->transpire[1];

		/* if MINIMAL_LAI_REDUCTION < change < 0, increase lai change */
		change = excess / (total_transpiration / m_maxlai);
		if (change > MINIMAL_LAI_REDUCTION) 
    {
			m_maxlai += MINIMAL_LAI_REDUCTION;
			if (m_maxlai < 0.0) m_maxlai = 0.0;
		}
		else 
    {	/* avoid prolonged recovery to roughness threshold */
			prev_woody_lai = m_maxlai;
			m_maxlai += change;
			if ((prev_woody_lai > mininum_tree) && (m_maxlai < mininum_tree)) m_maxlai = mininum_tree + MINIMAL_LAI_REDUCTION;
		}
	}
	else m_maxlai += MINIMAL_LAI;
  
	LightAttenuate(m_maxlai, grass_lai, lai, deciduous, curr_mo_State);
	
} // end of AdjustMaxlai()



/*
 *
 * ----------------------------  GrassWoodyEquilibrium
 *
 */
bool MAPSS_BiogeogModel::GrassWoodyEquilibrium(State site[], int mo,  
    int * months, bool * shrubs, float deficit[], bool deciduous, State curr_mo_State[])
{
	static State	last_year[MONTHS];
  static State * prev_year;
	static float threshold = MAX_DIFF_EQUILIB;
	static int cycles = 0;
  static int search = UNDIRECTED;
	float	excess;	/* water reserve above threshold */
	bool	DoLai = m_cats.mclass == Unknown;

	if (!DoLai) search = SEEK_FINAL;
	prev_year = &last_year[mo];
  
	if (*months >= MONTHS) // Check for equilibrium only every 12 months.
  {	/* define basis for acceptable equilibrium */
		if (*months == MONTHS) 
    {
			cycles = 0;
			threshold = MAX_DIFF_EQUILIB;
		}
		else if (cycles == MAX_CYCLES) 
    {
			/* adjust equilibrium thresholds for excessive cycles */
			cycles = 0;
			threshold = threshold + MAX_DIFF_INCREMENT;
		}
    
		if ((deficit[mo] > 0.0) && (m_maxlai > MINIMAL_LAI)) 
    {
			if (search >= SEEK_CEILING) search = SEEK_FINAL; /* deficit following excess of floor */
			else 
      {
				search = SEEK_FLOOR;	/* may update undirected search */
				if (DoLai) AdjustMaxlai(-deficit[mo], site, m_lai, m_grass_lai, deciduous, curr_mo_State);
				if (!*shrubs) pet_adjust();
				*months = JAN;				
				return(FALSE);
			}
		}

		if ((search != SEEK_FINAL) && (m_maxlai < m_LaiUpperBoundEnergy[*shrubs ? SHRUB : TREE]) &&
				((excess = minimum_reserve(site) - m_excess_h2o[INTERMEDIATE]) > 0.0)) 
    {
		/******
		if ((search != SEEK_FINAL) && (m_maxlai < *lai_ceiling) &&
				((excess = minimum_reserve(site) -
					m_excess_h2o[INTERMEDIATE]) > 0.0)) {
		*******/
			/*
			 * detected excess water reserve in driest month,
			 * increase lai
			 */
			if ((search) >= (SEEK_FLOOR)) search = SEEK_CEILING;		
			if (DoLai) AdjustMaxlai(excess, site, m_lai, m_grass_lai, deciduous, curr_mo_State);
			if (!*shrubs) pet_adjust();
			*months = JAN;			
			return(FALSE);
		}
		else 
    {	/* begin test for equilibrium */
			if (CheckStateValues(site, last_year, threshold, mo, 0.001f)) 
      { // All fields of site and prev_year are equal (or sorta equal) and same grass lai.
				if (!DoLai || NoDeficits(deficit) || (m_maxlai == 0.0)) 
        {
					if (!*shrubs && (m_maxlai < min_tree_lai[m_zone][m_cats.broadleaf])) 
          { // Shrubs above minimum tree lai.
						*shrubs = TRUE;
						pet_adjust();
						search = SEEK_FLOOR;
						if (DoLai) 
            {
							m_maxlai = m_LaiUpperBoundEnergy[SHRUB];
							LightAttenuate(m_maxlai, m_grass_lai, m_lai, deciduous, curr_mo_State);
						}
						*months = JAN;						
						return(FALSE);
					}
					if (threshold > MAX_DIFF_EQUILIB) 
          {
						m_equilibrium_factor = threshold;
						threshold = MAX_DIFF_EQUILIB;
					}
					else m_equilibrium_factor = 0;
					search = UNDIRECTED;	/* for next station */	
          				
					return(TRUE); /* This is the only TRUE for this routine */
				}
				else 
        {	/* deficits occur after ceiling hit */
					if (DoLai) 
          {
						if ((m_maxlai -= MINIMAL_LAI) < 0.0) m_maxlai = 0.0;
					}
					search = (search) < (SEEK_CEILING) ? SEEK_FLOOR : SEEK_FINAL;
					/*
					 * handles special case straddling
					 * roughness threshold
					 */
					if (DoLai) LightAttenuate(m_maxlai, m_grass_lai, m_lai, deciduous, curr_mo_State);
					if (!*shrubs) pet_adjust();
					*months = JAN;
					
					return(FALSE);
				}
			}
			else ++cycles;	/* unsuccessful equilibrium cycle */
		}
	} /* end of if (*months >= YEAR) */
  
	/* update stored state variables for next annual comparison */
	*prev_year = site[mo];
	++*months;		/* number of months processed */
	
	return(FALSE);
} // end of GrassWoodyEquilibrium()


/*
 *
 * ----------------------------  LightAttenuate
 *
 */
void MAPSS_BiogeogModel::LightAttenuate(float maxlai, float grass_lai[], float lai[][2], bool deciduous, State site[])
   {
	float			grasslai;
   double	common_subexp;
	
	if (maxlai <= no_attenuation_lai) grasslai = m_LaiUpperBoundEnergy[GRASS];
	else if (maxlai >= full_attenuation_lai)
      grasslai = 0.0;
	else 
      {
		common_subexp = m_LaiUpperBoundEnergy[GRASS] / (no_attenuation_lai - full_attenuation_lai);
		grasslai = (float)(((maxlai * common_subexp) + m_LaiUpperBoundEnergy[GRASS]) - (no_attenuation_lai * common_subexp));
		grasslai = MIN(grasslai, m_LaiUpperBoundEnergy[GRASS]);
      }
	NewMonthlyLai(maxlai, grass_lai, grasslai, lai, deciduous, site);
	
   } // end of LightAttenuate()


/*
 *
 * ----------------------------  minimum_reserve
 *
 */
float MAPSS_BiogeogModel::minimum_reserve(State * site)
{
	float	low_water, avail;
	int	i;

  low_water = m_soil.swhc[SURFACE] + m_soil.swhc[INTERMEDIATE];
  
	for (i = JAN; i <= DEC; i++) 
  {
		if (m_lai[i][WOODY] == m_maxlai) 
    {	/* growing season water availability */
			avail = site[i].soil_h2o[SURFACE] + site[i].soil_h2o[INTERMEDIATE] - m_unavail;
			if (low_water > avail) low_water = avail;
		}
	}
	
	return(low_water);
} // end of minimum_reserve()



/*
 *
 * ----------------------------  CheckStateValues
 * Compare values of current state variables with their values at the same month one year ago.
 */
bool MAPSS_BiogeogModel::CheckStateValues(State * curr_mo_State, State lastyear[], float threshold, int mo, float base)
{
	float			diff;

	diff = fabs(curr_mo_State->transpire[GRASS] - lastyear[mo].transpire[GRASS]);
	if ((diff > base) &&
			(diff > fabs(threshold * lastyear[mo].transpire[GRASS]))) {
		
		return(FALSE);
	}
	diff = fabs(curr_mo_State->evap - lastyear[mo].evap);
	if ((diff > base) && (diff > fabs(threshold * lastyear[mo].evap))) {
		
		return(FALSE);
	}
	diff = fabs(curr_mo_State->infiltrate - lastyear[mo].infiltrate);
	if ((diff > base) && (diff > fabs(threshold * lastyear[mo].infiltrate))) {
		
		return(FALSE);
	}
	diff = fabs(curr_mo_State->soil_h2o[SURFACE] -
		lastyear[mo].soil_h2o[SURFACE]);
	if ((diff > base) &&
			(diff > fabs(threshold * lastyear[mo].soil_h2o[SURFACE]))) {
		
		return(FALSE);
	}
	diff = fabs(curr_mo_State->soil_h2o[INTERMEDIATE] -
		lastyear[mo].soil_h2o[INTERMEDIATE]);
	if ((diff > base) &&
			(diff > fabs(threshold * lastyear[mo].soil_h2o[INTERMEDIATE]))) {
		
		return(FALSE);
	}
	diff = fabs(curr_mo_State->soil_h2o[DEEP] - lastyear[mo].soil_h2o[DEEP]);
	if ((diff > base) &&
			(diff > fabs(threshold * lastyear[mo].soil_h2o[DEEP]))) {
		
		return(FALSE);
	}
	diff = fabs(curr_mo_State->snow - lastyear[mo].snow);
	if ((diff > base) && (diff > fabs(threshold * lastyear[mo].snow))) {
		
		return(FALSE);
	}
	diff = fabs(curr_mo_State->surf_run - lastyear[mo].surf_run);
	if ((diff > base) && (diff > fabs(threshold * lastyear[mo].surf_run))) {
		
		return(FALSE);
	}
	diff = fabs(curr_mo_State->base_run - lastyear[mo].base_run);
	if ((diff > base) && (diff > fabs(threshold * lastyear[mo].base_run))) {
		
		return(FALSE);
	}
	
	return(TRUE);
}


/*
 *
 * ----------------------------  calc_hycon
 *
 */
double Soil_MAPSS::calc_hycon(double sand, double clay, double h2o, int layer)
{
	double			hycon;

	/********
		I am trying to think of some additional way to do some
		indirection into these arrays.
	*********/
	hycon = 2.778e-6 *
		exp(pp[layer] + (qq[layer] * sand) +
			((rr[layer] +
				(tt[layer] * sand) +
					(uu[layer] * clay) +
				(vv[layer] * clay * clay)) / h2o));
	/***
	hycon_sat = 2.778e-6 *
		exp(p + (q * sand) +
			((r + (t * sand) +(u * clay) + (v * clay * clay)) / h2o_sat));
	****/
	return(hycon);
}


/*
 *
 * ----------------------------  saturated_drainage
 *
 */
void Soil_MAPSS::saturated_drainage(State * site, int layer, float * destination /* h2o in destination layer */ )
{
	double Theta_vol, Theta_vol_fc, h2o, K_theta, QQQ_theta, reserve;
	
	Theta_vol = site->soil_h2o[layer];
	Theta_vol_fc = field_cap[layer];
	h2o = Theta_vol - Theta_vol_fc;
	if (h2o > 0.0) 
  {
    K_theta = KK[SATURATED][layer] * pow((float)(h2o/(swhc[layer] - Theta_vol_fc)), exp_perc[SATURATED][layer]);
    QQQ_theta = h2o * K_theta;
    if (layer < DEEP) 
    {
      reserve = swhc[layer + 1] - *destination;
      if (QQQ_theta > reserve) QQQ_theta = reserve;
    }
			
    *destination += (float)QQQ_theta; /* water flows to the next layer */			
    site->soil_h2o[layer] -= (float)QQQ_theta; /* water is removed from the current layer */
	}	
} // end of saturated_drainage()


/*
 *
 * ----------------------------  unsaturated_drainage
 *
 */
void Soil_MAPSS::unsaturated_drainage(State * site, int layer, float * destination /* h2o in destination layer */ )
{
	double		Theta_vol, K_theta, QQQ_theta, reserve;
	
	Theta_vol = site->soil_h2o[layer];
	if ((Theta_vol > wilt_pt[layer])) 
  {
    K_theta = KK[UNSATURATED][layer] *
				pow((float)((Theta_vol - wilt_pt[layer])/(swhc[layer] - wilt_pt[layer])), exp_perc[UNSATURATED][layer]);
    QQQ_theta = (site->soil_h2o[layer] - wilt_pt[layer]) * K_theta;
    
		if (layer < DEEP) 
    {
			reserve = swhc[layer + 1] - *destination;
			if (QQQ_theta > reserve) QQQ_theta = reserve;
		}
		
		*destination += (float)QQQ_theta; /* water flows to the next layer */		
		site->soil_h2o[layer] -= (float)QQQ_theta; /* water is removed from the current layer */

		if (site->soil_h2o[layer] < wilt_pt[layer]) 
        site->soil_h2o[layer] = wilt_pt[layer]; // ??? This seems to violate mass balance.
	}
} // end of unsaturated_drainage()

/*
 *
 * ----------------------------  CalcTrbxfr
 *
 *
*/ 
void MAPSS_BiogeogModel::CalcTrbxfr(int lifeform)
{
	int mo;
   float t;
   // CString msg;

	for (mo = JAN; mo <= DEC; mo++) 
   {
       if (m_vp_sat[mo]>m_vpr[mo])
       {
         t = (float)sciFn.trbxfr(c_z_all, modelParamsP->z0[m_zone][lifeform], m_tmp[mo], m_tmp[mo], 
               m_vpr[mo], m_vp_sat[mo], m_wnd[mo], m_elevation); 
         // trbxfr() returns a positive number when an error occurs.
       }
       else t = 0;
       if (t>=0.0 || m_vp_sat[mo]<m_vpr[mo]) 
       {
#ifdef UNIX_MAPSS
           // printf("*** mc2_functions.c/CalcTrbxfr(): mo, m_vp_sat[mo], m_vpr[mo] t = %d, %f, %f, %f\n",
               // mo, m_vp_sat[mo], m_vpr[mo], t);
#else
           CString msg;
           msg.Format("*** mc2_functions.c/CalcTrbxfr(): mo, m_vp_sat[mo], m_vpr[mo] t = %d, %f, %f, %f\n",
               mo, m_vp_sat[mo], m_vpr[mo], t);
           Report::InfoMsg(msg);
#endif
            t = 0.;
            // ASSERT(0);
       }
	    m_pet_array[lifeform][mo] = t * SECONDS_PER_DAY * sciFn.days_per_mo[mo] * -1.0f;
       ASSERT(m_pet_array[lifeform][mo]>=0.0);
   }
} // end of CalcTrbxfr()


/*
 *
 * ----------------------------  InitSoilWater
 * Start off the first month of the water year (START_MONTH, e.g. October)
 * with the top 2 soil layers at wilting point, the deep layer at
 * field capacity, and no snow on the ground.
 */
void MAPSS_BiogeogModel::InitSoilWater(State * curr_year)
{

	State	* start_month;

  ASSERT((START_MONTH - 1)>=0);
	start_month = &curr_year[START_MONTH - 1];

	/* add unavailable water */
	start_month->soil_h2o[SURFACE] = m_soil.wilt_pt[SURFACE];
	start_month->soil_h2o[INTERMEDIATE] = m_soil.wilt_pt[INTERMEDIATE];

	/* start subsoil "water table" at field capacity */
	start_month->soil_h2o[DEEP] = (m_soil.field_cap[DEEP] > m_soil.wilt_pt[DEEP]) ? 
      m_soil.field_cap[DEEP] : m_soil.swhc[DEEP]; // ??? shouldn't this be min(wilt_pt, swhc)?

	start_month->snow = 0.0;
	
} // end of MAPSS_BiogeogModel::InitSoilWater()


/*
 *
 * ----------------------------  InitSnow
 * Starting with the warmest month of the year, calculate the snow depth for each month
 * for 12 months. 
 */
void MAPSS_BiogeogModel::InitSnow(State * curr_year)
{
	int	mo;
	
	mo = m_max_tmp_mo;
  ASSERT(PREV_MO(mo)>=0);
	curr_year[PREV_MO(mo)].snow = 0.0;
	do 
  {
    if (m_months[mo].melt_rate > (curr_year[PREV_MO(mo)].snow + m_months[mo].snow))
    { // If there is any snow to start with, it all melts off this month
			curr_year[mo].snow = 0.0;
			curr_year[mo].melt = curr_year[PREV_MO(mo)].snow + m_months[mo].snow;
    }
    else 
    { // Some of the snow melts this month
      curr_year[mo].snow = curr_year[PREV_MO(mo)].snow + m_months[mo].snow - m_months[mo].melt_rate;
      curr_year[mo].melt = m_months[mo].melt_rate;
    }

		mo = NEXT_MO(mo);
	} while (mo != m_max_tmp_mo);
  
	
} // end of MAPSS_BiogeogModel::InitSnow(State * curr_year);


/*
 *
 * ----------------------------  LaiCycle
 *
 */
void MAPSS_BiogeogModel::LaiCycle(float * offsetP, State site_begin[MONTHS*2], float offset_begin, 
    State ** curr_yearP, State ** prev_yearP)
{
	unsigned sl;
	int	months, mo, lf, this_mo;
	int	inc_mode;
	float	* currlai;
  State * curr_year;
  State * prev_year;
	State * curr_mo_State;
	State *	prev_mo;
  bool unbalanced;
  int LAIiter, H2Oiter;
  float prev_maxlai;

  ASSERT(m_cats.mclass==Unknown);
  
	/* set lai iteration cycle site pointers, search mode */
  inc_mode = BINARY_MODE;
	curr_year = site_begin;
	prev_year = &site_begin[MONTHS];
  InitConductance(m_conductance, cond_max[m_zone][GRASS][m_cats.broadleaf], cond_max[m_zone][TREE][m_cats.broadleaf]);
  LAIiter = 0;
  do
  { // Loop to bring LAI into balance with water
    months = 0;			/* start each iteration with 0 months */
    InitSoilWater(curr_year); // Initialize soil water in the first month of the water year, START_MONTH, e.g. Oct.
    InitSnow(curr_year); // Calculate the snow depth for the 12 months starting with the warmest month.
    
    pet_adjust();
    if (diagsMAPSS && LAIiter<0) 
    {
      printf("MC2 m_petP - m_pet_array = %ld\n", m_petP - &(m_pet_array[0][0]));
      for (mo = 0; mo<12; mo++) 
      {
        printf("MC2 *(m_petP+%d) = %f\n", mo, *(m_petP+mo));
        // printf("MC2 m_pet_array[][%d] = %f, %f, %f\n", mo, m_pet_array[0][mo], m_pet_array[1][mo], m_pet_array[2][mo]);
      }
    }

    lf = (m_pet_adj == tree_pet_factor) ? TREE : SHRUB;
    H2Oiter = 0; 
    do 
    {	// Loop to establish annual cycle equilibrium at the given LAI
      mo = SetCurrPrev(months, curr_year, &curr_mo_State, &prev_mo, &currlai); // Establish pointers to current and previous
          // months, based on the number of iterations so far. On return, curr_mo_State points to curr_year[mo].
      DistributePpt(curr_mo_State, currlai, mo); // Preliminary partition of this month's rain between evap and infiltrate.
      for (sl = SURFACE; sl <= DEEP; sl++) curr_mo_State->soil_h2o[sl] = prev_mo->soil_h2o[sl];
      if (diagsMAPSS && LAIiter==0 && H2Oiter<=112) 
      {
#ifndef UNIX_MAPSS
        diagMsg.Format("MC2 LAIiter, H2Oiter, mo, currlai[0], currlai[1], m_grass_lai[mo] = %d, %d, %d, %f, %f, %f\n",
            LAIiter, H2Oiter, mo, *currlai, *(currlai+1), m_grass_lai[mo]);
        Report::LogMsg(diagMsg);
#endif
      }
      transpiration(curr_mo_State, m_conductance, currlai, m_grass_lai, mo);
      H2Oiter++;
      if (diagsMAPSS && LAIiter==0 && H2Oiter<=112)
      {
#ifndef UNIX_MAPSS
        diagMsg.Format("MC2 LAIiter, H2Oiter, *(m_petP+mo), m_maxlai, curr_mo_State->soil_h2o[] = %d, %d, %f, %f, %f, %f, %f\n", 
            LAIiter, H2Oiter, *(m_petP+mo), m_maxlai, curr_mo_State->soil_h2o[0], curr_mo_State->soil_h2o[1], curr_mo_State->soil_h2o[2]);
        Report::LogMsg(diagMsg);
#endif
      }
    } while (!AnnualCycleEquilibrium(curr_mo_State, mo, &months));
    
    prev_maxlai = m_maxlai;
    unbalanced = AdjustLaiIfUnbalanced(&curr_year, &prev_year, &inc_mode, lf); // Puts the adjusted LAI into m_maxlai.
    LAIiter++;
    if (1 && diagsMAPSS && m_maxlai!=prev_maxlai)
    { 
#ifndef UNIX_MAPSS
      diagMsg.Format("MC2 LAIiter, months, prev_maxlai, m_maxlai, curr_mo_State->soil_h2o[] = %d, %d, %f, %f, %f, %f, %f\n", 
          LAIiter, months, prev_maxlai, m_maxlai, curr_mo_State->soil_h2o[0], curr_mo_State->soil_h2o[1], curr_mo_State->soil_h2o[2]);
      Report::LogMsg(diagMsg);
#endif
    }
    if (unbalanced)
    { // Prepare for next iteration.
      /* maintain pointer to previous end State */
      SwapYears(&curr_year, &prev_year);
      
      // Set the woody LAI to m_maxlai when the leaves are on, 0 otherwise.  
      // Grass LAI should be 0 until GrassWoody() is called.
      // ??? This is a little different from the original MAPSS code in MC1. It used to call NewMonthlyLai()      
      for (this_mo = JAN; this_mo <= DEC; this_mo++)
      {
        /*
        if (m_tmp[this_mo] >= frost) m_lai[this_mo][WOODY] = m_maxlai; // growing season 
        else m_lai[this_mo][WOODY] = m_cats.deciduous ? 0 : m_maxlai; // outside of growing season
        */
        if (!m_cats.deciduous || m_tmp[this_mo]>=frost) m_lai[this_mo][WOODY] = m_maxlai;
        else if (m_lai[this_mo][WOODY]>m_maxlai) m_lai[this_mo][WOODY] = 0.0;
        ASSERT((m_lai[this_mo][GRASS]==0.) && (m_grass_lai[this_mo]==0.));
      }
    }
	} while (unbalanced);
  // Now the LAI is in balance w/ water.
  *curr_yearP = curr_year;
  *prev_yearP = prev_year;
  
  if (diagsMAPSS) 
  {
#ifdef UNIX_MAPSS
     printf("MC2 LAIiter, H2Oiter, mo, currlai[0], currlai[1], m_grass_lai[mo] = %d, %d, %d, %f, %f, %f\n",
         LAIiter, H2Oiter, mo, *currlai, *(currlai+1), m_grass_lai[mo]);
     for (int i=0; i<12; i++) printf("MC2 LaiCycle: m_grass_lai[%d] = %f\n", i, m_grass_lai[i]);
#else
     diagMsg.Format("MC2 LaiCycle: LAIiter, H2Oiter = %d, %d\n", LAIiter, H2Oiter);
     Report::LogMsg(diagMsg);
#endif
  }

} // end of LaiCycle(...)



/*
 * 
 * ----------------------------  BroadDecid
 * Distinguish between EB, DB, and EN
 */
void MAPSS_BiogeogModel::BroadDecid(Rulecat * cats, float productivity, bool IsTree, bool SecondTime)
{
	float	broad_ppt;
	
  ProcessSeason(frost, &broad_ppt); // calculate minimum monthly precip during growing season

	if (m_min_tmp < n_decid_bound) cats->broadleaf = cats->deciduous = FALSE; // evergreen needleleaf
	
  else if (m_thaw == NEVER_FREEZES) 
  { // evergreen broadleaf
    cats->broadleaf = TRUE;
    cats->deciduous = FALSE;	
	}
	
  else if ((m_mix_ratio >= evergreen_gdd_ratio) || (broad_ppt < broad_ppt_min) || (SecondTime
      && ((m_growing_deg_days_frost <= evergreen_gdd) 
      || (productivity <= (IsTree ? evergreen_productivity_tree : evergreen_productivity_shrub)))))
      cats->broadleaf = cats->deciduous = FALSE; // again, evergreen needleleaf
  
  else cats->broadleaf = cats->deciduous = TRUE; // deciduous broadleaf

} // end of BroadDecid()


/*
 *
 * ----------------------------  SetCurrPrev
 * Advance the pointers by one month, i.e. "Set Current and Previous Months".
 * Incoming value of "months" is a month counter, i.e. number of months processed so far, and it starts w/ 0.
 * Returned value is a month index, i.e. Jan=0, Feb=1, ...
 * The first month processed is the month which begins the water year, specified as START_MONTH.
 */
int MAPSS_BiogeogModel::SetCurrPrev(int months, State site[], State * curr_mo_State[], State * prev_mo[], float * * currlai)
{
  int mo;
	
	mo = (months + START_MONTH) % MONTHS;
	*curr_mo_State = &site[mo];
	*currlai = m_lai[mo]; // m_lai is dimensioned [MONTHS][2], so currlai points to a pair of values, one for GRASS
      // and one for WOODY

	*prev_mo = &site[PREV_MO(mo)];
	
	return(mo);
} // end of SetCurrPrev()




/*
 *
 * ----------------------------  DistributePpt
 *
 */
void MAPSS_BiogeogModel::DistributePpt(State * curr_mo_State, float * lai, int mo)
{
	float	intercept, interc_event, rain;
  float	lai_sum = lai[WOODY] + lai[GRASS];
	
	rain = m_ppt[mo] - m_months[mo].snow;
	interc_event = interc_lai * (1.0f - exp(-lai_sum));
	intercept = MIN(interc_event * m_months[mo].events, rain);
	curr_mo_State->evap = MIN(intercept, *(m_petP + mo));

	/* infiltration estimate prior to soil water uptake */
	curr_mo_State->infiltrate = curr_mo_State->melt + rain - curr_mo_State->evap;
	curr_mo_State->surf_run = curr_mo_State->base_run = 0.0;
	
} // end of DistributePpt()



/*
 *
 * ----------------------------  transpiration
 *
 */ 
float MAPSS_BiogeogModel::transpiration(State * curr_mo_State, float conductance[][2], float lai[], float grass_lai[], int mo)
{
  int sl;
	float	tot_infiltrate; 
  float excess_pot_at;	/* unsatisfied tree pot_at */
	float	water;

  tot_infiltrate = curr_mo_State->infiltrate;
  curr_mo_State->surf_run = curr_mo_State->infiltrate * pow((curr_mo_State->soil_h2o[SURFACE] - m_soil.wilt_pt[SURFACE])/m_soil.sat2mat[SURFACE],
			m_soil.k_surfrun);
	curr_mo_State->infiltrate -= curr_mo_State->surf_run;

	infiltrate(curr_mo_State); // First infiltration

  excess_pot_at = TranspireStep(curr_mo_State, conductance, lai, &grass_lai[mo], mo);
  
	/***
		1st saturated flow
	****/
	m_soil.saturated_drainage(curr_mo_State, DEEP, &curr_mo_State->base_run);
	m_soil.saturated_drainage(curr_mo_State, INTERMEDIATE, &curr_mo_State->soil_h2o[DEEP]);
	m_soil.saturated_drainage(curr_mo_State, SURFACE, &curr_mo_State->soil_h2o[INTERMEDIATE]);

	/***
		1st unsaturated flow
	****/
	m_soil.unsaturated_drainage(curr_mo_State, DEEP, &curr_mo_State->base_run);
	m_soil.unsaturated_drainage(curr_mo_State, INTERMEDIATE, &curr_mo_State->soil_h2o[DEEP]);
	m_soil.unsaturated_drainage(curr_mo_State, SURFACE, &curr_mo_State->soil_h2o[INTERMEDIATE]);

	/***
		reset soil water reservoirs with "post-transpiration"
		infiltration 2nd infiltrate
	****/
	if (curr_mo_State->infiltrate > 0.0) infiltrate(curr_mo_State);

	/***
		2nd saturated flow
	****/
	m_soil.saturated_drainage(curr_mo_State, DEEP, &curr_mo_State->base_run);
	m_soil.saturated_drainage(curr_mo_State, INTERMEDIATE, &curr_mo_State->soil_h2o[DEEP]);
	m_soil.saturated_drainage(curr_mo_State, SURFACE, &curr_mo_State->soil_h2o[INTERMEDIATE]);

	/***
		2nd unsaturated flow
	****/
	m_soil.unsaturated_drainage(curr_mo_State, DEEP, &curr_mo_State->base_run);
	m_soil.unsaturated_drainage(curr_mo_State, INTERMEDIATE, &curr_mo_State->soil_h2o[DEEP]);
	m_soil.unsaturated_drainage(curr_mo_State, SURFACE, &curr_mo_State->soil_h2o[INTERMEDIATE]);

	/***
		3rd infiltrate -- post-drainage, final infiltration
	****/
	if (curr_mo_State->infiltrate > 0.0) infiltrate(curr_mo_State);

	/***
		3rd saturated flow
	****/
	m_soil.saturated_drainage(curr_mo_State, DEEP, &curr_mo_State->base_run);
	m_soil.saturated_drainage(curr_mo_State, INTERMEDIATE, &curr_mo_State->soil_h2o[DEEP]);
	m_soil.saturated_drainage(curr_mo_State, SURFACE, &curr_mo_State->soil_h2o[INTERMEDIATE]);

	/***
		3rd unsaturated flow
	****/
	m_soil.unsaturated_drainage(curr_mo_State, DEEP, &curr_mo_State->base_run);
	m_soil.unsaturated_drainage(curr_mo_State, INTERMEDIATE, &curr_mo_State->soil_h2o[DEEP]);
	m_soil.unsaturated_drainage(curr_mo_State, SURFACE, &curr_mo_State->soil_h2o[INTERMEDIATE]);

	/***
		remaining infiltrate is surface runoff
	****/
	curr_mo_State->surf_run += curr_mo_State->infiltrate;

	/***
		reset infiltrate for display of simulation
	****/
	curr_mo_State->infiltrate = tot_infiltrate;

	/***
		Hack for rounding errors that may cause sqrt DOMAIN error in the next simulation month.
		Note that this can "create" water in the system. ???
	****/
	for (sl = SURFACE; sl <= DEEP; sl++) if (curr_mo_State->soil_h2o[sl] < m_soil.wilt_pt[sl]) 
			curr_mo_State->soil_h2o[sl] = m_soil.wilt_pt[sl];
		
	if (excess_pot_at != 0.0) return(excess_pot_at + m_inadeq_h2o[INTERMEDIATE]);

	water = curr_mo_State->soil_h2o[SURFACE] + curr_mo_State->soil_h2o[INTERMEDIATE] - m_unavail;
	if (water < m_inadeq_h2o[INTERMEDIATE]) return(m_inadeq_h2o[INTERMEDIATE] - water);
  
  return(0.0);	/* indicates balance, precludes change in tree lai */
} // end of transpiration()




/*
 *
 * ----------------------------  AnnualCycleEquilibrium
 *
 */
bool MAPSS_BiogeogModel::AnnualCycleEquilibrium(State * curr_mo_State, int mo, int * months)
{
	static State	last_year[MONTHS];
	static float	threshold = MAX_DIFF_EQUILIB;
	static int		cycles = 0;
	
	/* adjust equilibrium thresholds for excessive cycles */
	if (cycles == MAX_CYCLES) 
  {
		cycles = 0; // Start the cycle count over.
		threshold = threshold + MAX_DIFF_INCREMENT; // Increase the threshold.
	}
  
	if (*months >= MONTHS) // Check for equilibrium only after the end of the first year.
  {	/* define basis for acceptable equilibrium */
		if (CheckStateValues(curr_mo_State, last_year, threshold, mo, 0.0)) 
    { // Equilibrium has been achieved.
      // Reset static variables for next time, and return TRUE.
			cycles = 0;
			if (threshold > MAX_DIFF_EQUILIB) 
      { // The threshold had been raised.
				m_equilibrium_factor = threshold; // Remember what the threshold was raised to.
				threshold = MAX_DIFF_EQUILIB; // Reset the threshold to its starting value.
			} 
      else m_equilibrium_factor = 0; // The threshold was at its nominal value.
			return(TRUE);
		}
		else ++cycles;	/* unsuccessful equilibrium cycle */
	}
  
	/* update stored state variables for next comparison */
	last_year[mo] = *curr_mo_State;
	++*months;		/* number of months processed */

	return(FALSE);
} // end of AnnualCycleEquilibrium()




/*
 *
 * ----------------------------  AdjustLaiIfUnbalanced
 * Tests whether LAI is in balance with water, and if not, adjusts LAI.
 */
/****
	Determine if iteration is to continue, and by which function
	(binary offset or increment) return TRUE to continue lai seek,
	return FALSE to end process.

*****/
bool MAPSS_BiogeogModel::AdjustLaiIfUnbalanced(State * * curr_year, State * * prev_year, int * inc_mode, int lifeform)
/****
	inc_mode:	1 BINARY_MODE -> binary search mode
            0 INCREMENT_MODE aka CANNOT_EXCEED_EXCESS_H2O -> increment mode, h2o cannot exceed m_excess_h2o
*****/
{
	float	low_water;
	float	minimum_tree;
	
	minimum_tree = min_tree_lai[m_zone][m_cats.broadleaf];
  low_water = minimum_reserve(*curr_year); // Find the least soil water in the top 2 layers over the year.

	if (*inc_mode == BINARY_MODE) 
  {
		m_offset *= 0.5;
		if (low_water > m_excess_h2o[INTERMEDIATE]) 
    { // low_water is > EXCESS_H2O fraction (e.g. 0.05) of awc
			if (INCREMENT > (m_offset + m_offset)) return(FALSE); // Stop iterating without going into increment mode.
			else 
      { // Increase the LAI and continue the iteration.
				m_maxlai += m_offset;
				m_maxlai = MIN(m_maxlai, m_LaiUpperBoundEnergy[lifeform]);
			}
		}
		else if (low_water > m_inadeq_h2o[INTERMEDIATE]) 
    { // low_water is between INADEQ_H2O fraction and EXCESS_H2O fraction (e.g. between 0.025 and 0.05) of awc
			// Switch to increment mode, increase the LAI by a small increment, and continue the iteration. 
			*inc_mode = CANNOT_EXCEED_EXCESS_H2O;
			m_maxlai += INCREMENT; // e.g. 0.1 LAI unit
			m_maxlai = MIN(m_maxlai, m_LaiUpperBoundEnergy[lifeform]);
		}
		else 
    { // low_water is <= INADEQ_H2O fraction (e.g. 0.025) of awc
			if (m_maxlai < MINIMAL_LAI) return(FALSE); // lai approximately 0.0, so stop iterating. 
      // Decrease the LAI and continue the iteration.
			if (INCREMENT > (m_offset + m_offset)) m_offset += m_offset;
			if ((m_maxlai > minimum_tree) && (lifeform == TREE)) 
      {
				if (m_maxlai <= (minimum_tree + INCREMENT)) m_maxlai -= m_offset;
				else if ((m_maxlai - m_offset) < minimum_tree) m_maxlai -= MIN(m_offset, (m_maxlai - minimum_tree) - INCREMENT);
        else m_maxlai -= m_offset;
			}
			else m_maxlai -= m_offset;
      ASSERT(m_maxlai>=0.);
		}	    	
		return(TRUE); // Continue the iteration.
	}
	else 
  {		/* increment mode, not binary mode.  */
    ASSERT(*inc_mode==INCREMENT_MODE);
		if (low_water <= m_inadeq_h2o[INTERMEDIATE]) 
    {	// low_water is <= INADEQ_H2O fraction (e.g. 0.025) of awc
      /* end iteration, backup one cycle to previous lai level */
			m_maxlai -= INCREMENT;
			SwapYears(curr_year, prev_year); /* point to previous State */			
			return(FALSE); // Stop iterating because all the water has been used up.
		}
    m_maxlai += INCREMENT;
		if (m_maxlai > m_LaiUpperBoundEnergy[lifeform]) 
    { /* do not exceed maximum lai */
			m_maxlai = MIN(m_LaiUpperBoundEnergy[lifeform], m_maxlai - INCREMENT);			
			return(FALSE); // Stop iterating because the LAI has hit the maximum.
		}		
		return(TRUE); // Continue the iteration.
	}
} // end of AdjustLaiIfUnbalanced()


/*
 *
 * ----------------------------  InitConductance
 *
 */
void MAPSS_BiogeogModel::InitConductance(float conductance[][2], float grass, float woody)
{
	int			mo;
	
	for (mo = JAN; mo <= DEC; mo++) 
  {
		conductance[mo][GRASS] = grass;
		conductance[mo][WOODY] = woody;
	}	
} // end of InitConductance()


/*
 *
 * ----------------------------  GrassWoodyPart1
 *
 */
void MAPSS_BiogeogModel::GrassWoodyPart1(State curr_year[], State prev_year[], State site_begin[], float offset_begin)
{
	float	* currlai,	/* at end of cycle, points to output lai array
								during cycle, points to next month's lai */
					deficit[MONTHS];	/* water available to trees */
	State			*site,
					*curr_mo_State,
					*prev_mo;	/* current, prev month within site array */
	unsigned		j;
	int				months = JAN,	/* tree lai can change if months >= YEAR */
					mo;
	bool			shrubs = FALSE,
						/* if TRUE, shrubs > 3.0 lai, no pet_adjust */
					deciduous = m_cats.deciduous &&
						m_cats.broadleaf ? TRUE : FALSE,
					broadleaf = m_cats.broadleaf;
	// float			c3c4_ratio = 0.0;
	float			delta_lai = 0.0;
	float			winter_ppt_avg;
   float pet_ann;

  ASSERT(m_cats.mclass==Unknown);	
	m_lai_values[GRASS] = m_lai_values[TREE] = m_lai_values[SHRUB] = m_lai_values[GRASS_SUM] = m_lai_values[TREE_SUM] =
      m_lai_values[SHRUB_SUM] = 0.0;
  // c3c4_ratio = 0.;
  m_at_woody_ppt = 0.0;
  m_at_woody_ppt_norm = 0.0;
  m_woody_alone_at = 0.;
  m_k_factor = 0.;
	site = curr_year;
  shrubs = !(m_pet_adj == tree_pet_factor);
	
	m_old_tree_lai = m_maxlai;

	if ((!shrubs) && (k_factor_pet_boundary[m_zone][broadleaf ? 1 : 0] > 0.0) && (!m_cats.grassfire)) 
  {
		pet_ann = 0.0;
		for (mo = JAN; mo <= DEC; mo++) 
    {
      if (0 && diagsMAPSS) printf("MC2 site[%d].transpire[], pet[%d] = %f, %f, %f\n", 
          mo, mo, site[mo].transpire[0], site[mo].transpire[1], *(m_petP + mo));
      m_woody_alone_at += site[mo].transpire[WOODY];
			pet_ann += *(m_petP + mo);
		}
    if (m_ppt_winter > 0.0) m_at_woody_ppt = m_woody_alone_at / m_ppt_winter;
		else m_at_woody_ppt = 0.0;

		if ((m_ppt_winter > 0.0) && (m_length_winter > 0)) m_at_woody_ppt_norm = 
        (m_woody_alone_at /((float) (MONTHS - m_length_winter))) 
        / (m_ppt_winter /((float) m_length_winter));
		else m_at_woody_ppt_norm = 0.0;

		if (m_length_winter > 0) winter_ppt_avg = m_ppt_winter / ((float) m_length_winter);
		else winter_ppt_avg = 0.0;

		if (k_factor_slope[m_zone][broadleaf ? 1 : 0] > 0.0) m_k_factor = k_factor_slope[m_zone][broadleaf ? 1 : 0] * ((m_at_woody_ppt_norm / k_factor_constraint) - 1.0f);
		else m_k_factor = k_factor_slope[m_zone][broadleaf ? 1 : 0] * ((winter_ppt_avg / k_factor_constraint) - 1.0f);

		if ((m_length_winter > k_factor_winter_boundary_lower) && (pet_ann > k_factor_pet_boundary[m_zone][broadleaf ? 1 : 0])) 
    {
			delta_lai = m_k_factor * m_at_woody_ppt_norm;
			m_LaiUpperBoundEnergy[TREE] = MAX(m_maxlai - delta_lai, 0.0f);
			/*******
				This NEEEEEEEDS to be done here just in case delta_lai
				was a large negative number.  We could have set the
				max tree lai somewhere in the range of MAXINT.
			********/
			m_LaiUpperBoundEnergy[TREE] = MIN(m_LaiUpperBoundEnergy[TREE], LaiUpperBoundsEnergy[m_zone][TREE]);
			/*******
				I have to do this here, now.  If the tree lai was
				beaten way back to something less than min_tree_lai[],
				then I have to redo all the tree vs shrub stuff.
			********/
			pet_adjust();
			shrubs = !(m_pet_adj == tree_pet_factor);
			m_maxlai = m_LaiUpperBoundEnergy[TREE];
			for (mo = JAN; mo <= DEC; mo++) m_lai[mo][WOODY] = MIN(m_lai[mo][WOODY], m_maxlai);
		}
	}
  
	m_new_tree_lai = m_maxlai;

	if ((m_maxlai < full_attenuation_lai) || (m_maxlai < min_tree_lai[m_zone][broadleaf])) 
  {
    int H2Oiter;
    
		LightAttenuate(m_maxlai, m_grass_lai, m_lai, deciduous, site);
		InitSoilWater(site);
		InitSnow(curr_year);
    H2Oiter = 0;
		do 
    {	/* h2o equilibrium cycle */
			mo = SetCurrPrev(months, site, &curr_mo_State, &prev_mo, &currlai);
			DistributePpt(curr_mo_State, currlai, mo);
			for (j = SURFACE; j <= DEEP; j++) curr_mo_State->soil_h2o[j] = prev_mo->soil_h2o[j];
			deficit[mo] = transpiration(curr_mo_State, m_conductance, currlai, m_grass_lai, mo);
      H2Oiter++;
      if (0 && diagsMAPSS) printf("MC2 GrassWoodyPart1 H2Oiter, mo, deficit[mo], m_maxlai, m_grass_lai[mo] = %d, %d, %f, %f, %f\n", 
          H2Oiter, mo, deficit[mo], m_maxlai, m_grass_lai[mo]);
		} while (!GrassWoodyEquilibrium(site, mo, &months, &shrubs, deficit, deciduous, site));
	}
  
	if (m_maxlai < full_attenuation_lai) 
      m_lai_values[GRASS] = FindMeanLai(&m_lai_values[GRASS_SUM], GRASS); /* calculate mean grass lai */
} // end of GrassWoodyPart1()


/*
 *
 * ----------------------------  SwapYears
 *
 */
void MAPSS_BiogeogModel::SwapYears(State * * curr_year, State * * prev_year)
{
	State			*tmp;

	tmp = *curr_year;
	*curr_year = *prev_year;
	*prev_year = tmp;
	
} // end of SwapYears()



/*
 *
 * ----------------------------  NewMonthlyLai
 *
 */
void MAPSS_BiogeogModel::NewMonthlyLai(float maxlai, float grass_lai[], float grasslai, float lai[][2],
    bool deciduous, State site[])
{
	int  			mo;
	float			spring_lai = MIN(grasslai, spring_grass);
	bool			had_spring = FALSE;
	
	/* assign lai for woody lifeform based on temperature condition */
	for (mo = JAN; mo <= DEC; mo++) 
  {
		if (m_growing_deg_days_zero <= taiga_tundra_boundary[m_zone]) 
    { // Treat grass differently if it a thermal bound m_zone.
			if ((m_tmp[mo] >= frost) && (site[mo].snow <= grass_snowpack)) 
      {
				if (!had_spring) 
        {
					lai[mo][GRASS] = grass_lai[mo] = spring_lai;
					had_spring = TRUE;
				}
				else lai[mo][GRASS] = grass_lai[mo] = grasslai;
			}
			else lai[mo][GRASS] = grass_lai[mo] = 0.0;
		}
		else 
    {
			if (m_tmp[mo] >= frost) 
      {
				if (m_tmp[PREV_MO(mo)] < frost) lai[mo][GRASS] = grass_lai[mo] = spring_lai; // If the prev month tmp is 
            // below frost and this month is above frost, apply grass.
				else lai[mo][GRASS] = grass_lai[mo] = grasslai;
			}
			else lai[mo][GRASS] = grass_lai[mo] = 0.0;
		}
		if (!deciduous || (m_tmp[mo] >= frost)) lai[mo][WOODY] = maxlai;
		else if (lai[mo][WOODY] > maxlai) lai[mo][WOODY] = 0.0; // Reduces non-growing season lai.
	}	
} // end of NewMonthlyLai()


/*
 * 
 * ----------------------------  ProcessSeason
 * 
 */
/***
	the min calculated here is not instantaneous, but is the
	min for all possible one month windows ending between the begin
	time and end time of the season, inclusive
****/
void MAPSS_BiogeogModel::ProcessSeason(float frost_line, float * ppt_ptr)
{
	float			min = 1000.0;	/* minimum precip rate during season */
	int				mo;
	Interval		period;

	
	for (mo = JAN; mo <= DEC; mo++) 
  {
		TimeInterval(m_tmp[mo], m_tmp[NEXT_MO(mo)], frost_line, &period);
		if (period.end != 0.0) MonthMin(m_ppt[mo], m_ppt[NEXT_MO(mo)], m_ppt[PREV_MO(mo)], &min, &period);
	}

	if (min == 1000.0) min = -1.0;			/* season does not exist */

	*ppt_ptr = min;
	
} // end of ProcessSeason()


/*
 *
 * ----------------------------  infiltrate
 *
 */
void MAPSS_BiogeogModel::infiltrate(State * site)
{
	float			deficit;	/* saturation - current soil h2o */
	int				i;
	
	for (i = SURFACE; i <= DEEP; i++) 
  {
    deficit = m_soil.swhc[i] - site->soil_h2o[i];
    if (deficit<=0.0) continue; // No room for more water in this soil layer.
    if (deficit<=site->infiltrate)
    { // Fill up this layer completely, and reduce the infiltrate accordingly.
      site->soil_h2o[i] = m_soil.swhc[i];
      site->infiltrate -= deficit;
      if (site->infiltrate<0.) site->infiltrate = 0.; // Deal with potential roundoff.
    }
    else 
    { // Add all the remaining infiltrate to this layer.
      site->soil_h2o[i] += site->infiltrate;
      site->infiltrate = 0.;
      break; // The infiltrate is all gone, so there is no need to look at lower levels.
    }
  }	
} // end of infiltrate()


/*
 *
 * ----------------------------  TranspireStep
 *
 */
float MAPSS_BiogeogModel::TranspireStep(State * curr_mo_State, float conductance[][2], float lai[], float * grass_lai, int mo)
	/******
	lifeform GRASS transpiration precedes lifeform WOODY
	because lifeform WOODY will draw from both strata
	*******/
{
	int	veg, lifeform, i;
	float	pot_at[2];
	float	* curr_conduct; // pointer to this month's conductance pair 
	float	begin_h2o; // h2o in upper soil at start 
  float pot_transp; // potential transpiration 
  float remaining = 0.0; // water remaining in upper soil layer 
  float h2o_fraction[2]; // fraction upper h2o available by lifeform 
  float mm_transpired;
	bool	shrub, broadleaf;
	float	curr_canopy_cond_max;
  float potential_water_use[2];
  float sum;
  float pet;
  double swp_orig, swp_MC2, pct_soil_h2o_orig, pct_soil_h2o_MC2;
	
  pet = *(m_petP + mo);
  shrub = !(m_pet_adj == tree_pet_factor);
  broadleaf = m_cats.broadleaf;
	curr_conduct = conductance[mo];
  begin_h2o = curr_mo_State->soil_h2o[SURFACE];
  pot_transp = pet - curr_mo_State->evap;
  
	// Select curr_canopy_cond_max from cond_surface_max[][][] of largest lifeform.
	for (veg = GRASS; veg <= WOODY; veg++) if (lai[veg] > 0.0) curr_canopy_cond_max = cond_surface_max[m_zone][LIFEFORM(veg)];
	for (veg = GRASS; veg <= WOODY; veg++) 
  {
    swp_orig = curr_mo_State->swp[veg];
    pct_soil_h2o_orig = curr_mo_State->pct_soil_h2o[veg];
    curr_conduct[veg] = StomatalConductanceMC2(curr_mo_State, pet, veg, conductance, mo, shrub, broadleaf);
    swp_MC2 = curr_mo_State->swp[veg];
    pct_soil_h2o_MC2 = curr_mo_State->pct_soil_h2o[veg];
    curr_mo_State->swp[veg] = swp_orig;
    curr_mo_State->pct_soil_h2o[veg] = pct_soil_h2o_orig;
    curr_mo_State->swp[veg] = swp_MC2;
    curr_mo_State->pct_soil_h2o[veg] = pct_soil_h2o_MC2;
  }
  // printf("MC2 stomatal conductance = %.8f, %.8f\n", curr_conduct[GRASS], curr_conduct[WOODY]);

	do 	/* while (!GrassPotAtSatisfied(... */
	{	
		for (veg = GRASS; veg <= WOODY; veg++) 
      { float kTransp, laiVeg, currConductVeg, condLAImax;
         lifeform = LIFEFORM(veg);
         //      pot_at[veg] = lai[veg]>0.f ? 
         //          pot_transp*(1.0f - exp(-k_transp[m_zone][lifeform][0][broadleaf]*lai[veg]*curr_conduct[veg]
         //          /cond_lai_max[m_zone][lifeform][m_cats.broadleaf]))
         //          : 0.f;
         kTransp = k_transp[m_zone][lifeform][0][broadleaf];
         laiVeg = lai[veg];
         currConductVeg = curr_conduct[veg];
         condLAImax = cond_lai_max[m_zone][lifeform][m_cats.broadleaf];
         pot_at[veg] = laiVeg>0.f ? 
            pot_transp*(1.0f - exp(-kTransp*laiVeg*currConductVeg
            /condLAImax))
            : 0.0f;
      }
      /* Allocate water between grass and woody lifeforms. */
      sum = 0.;
      for (i = GRASS; i <= WOODY; i++) 
      {
         lifeform = LIFEFORM(i);
         potential_water_use[i] = lai[i]>0.f ? lai[i]*curr_conduct[i]/cond_lai_max[m_zone][lifeform][m_cats.broadleaf] : 0.f;
         sum += potential_water_use[i];
      }
      if (lai[GRASS]<=0. || sum<=0.)
      { // usually winter 
         h2o_fraction[GRASS] = 0.0;
         h2o_fraction[WOODY] = 1.0;
      }
      else for (i = GRASS; i <= WOODY; i++) h2o_fraction[i] = potential_water_use[i]/sum;	

      /* If necessary, adjust transpiration to take into account available water. */
		if ((pot_at[GRASS] + pot_at[WOODY]) > pot_transp) for (veg = GRASS; veg <= WOODY; veg++) 
				pot_at[veg] = pot_transp * h2o_fraction[veg];

		/* transpiration associated with upper soil layer */
		for (veg = GRASS; veg <= WOODY; veg++) curr_mo_State->transpire[veg] = pot_at[veg]>0. ? 
         AccumulateTranspiration(curr_mo_State, &pot_at[veg], (begin_h2o - m_soil.wilt_pt[SURFACE])*h2o_fraction[veg], SURFACE)
         : 0.f;

		/* check if grass could utilize remaining h2o in upper soil */
		remaining = curr_mo_State->soil_h2o[SURFACE] - m_soil.wilt_pt[SURFACE];
		if ((remaining > 0.01) && (pot_at[GRASS] > 0.0)) curr_mo_State->transpire[GRASS] += 
         AccumulateTranspiration(curr_mo_State, &pot_at[GRASS], remaining, SURFACE);

		if ((curr_mo_State->infiltrate > 0.01) && (pot_at[GRASS] > 0.0)) 
      { // Remove water directly from the infiltration pool.
		   mm_transpired = MIN(pot_at[GRASS], curr_mo_State->infiltrate);
		   curr_mo_State->transpire[GRASS] += mm_transpired;
		   curr_mo_State->infiltrate -= mm_transpired;
		}
	} while (!GrassPotAtSatisfied(pot_at[GRASS], curr_mo_State, curr_conduct[GRASS], lai, grass_lai, begin_h2o, pot_transp));

	/* check if grass could utilize remaining h2o in upper soil */
	if ((remaining > 0.0) && (pot_at[GRASS] > 0.0)) curr_mo_State->transpire[GRASS] +=
			AccumulateTranspiration(curr_mo_State, &pot_at[GRASS], remaining, SURFACE);
        
	/* Transpire h2o from middle soil layer if woody plants need more water. */
	if (pot_at[WOODY] > 0.0) curr_mo_State->transpire[WOODY] += AccumulateTranspiration(curr_mo_State, &pot_at[WOODY], 
        curr_mo_State->soil_h2o[INTERMEDIATE] - m_soil.wilt_pt[INTERMEDIATE], INTERMEDIATE);

	return(pot_at[WOODY]);	/* return any excess pot_at of trees */
} // end of TranspireStep()


void MAPSS_BiogeogModel::pet_adjust()
{
  if (m_maxlai >= min_tree_lai[m_zone][m_cats.broadleaf ? 1 : 0])
  {
    m_petP = m_pet_array[TREE];
    m_pet_adj = tree_pet_factor;
  }
  else 
  {
    m_petP = m_pet_array[SHRUB];
    m_pet_adj = tallgrass_pet_factor;
  }
} // end of pet_adjust()


