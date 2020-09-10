/*
 *  ScienceFcns.h
 *
 */

#ifndef SCIENCEFCNS_H
#define SCIENCEFCNS_H

#ifdef UNIX_MAPSS
#ifndef Boolean
#define Boolean int
#endif
#endif

#ifndef FALSE
#define FALSE false
#endif

#ifndef TRUE
#define TRUE true
#endif

# ifndef MIN
#  define MIN(x, y)	(((x) < (y)) ? (x) : (y))
# endif				/* MIN */

# ifndef MAX
#  define MAX(x, y)	(((x) > (y)) ? (x) : (y))
# endif				/* MAX */

#define PI 3.14159265

#define SEA_LEVEL 1.013246e5 /* standard sea level pressure (Pa) */
#define GRAV 9.80665 /* gravitational acceleration at reference latitude 45d 32m 33s */
#define MOL_AIR 28.9644 /* molecular weight of air (g mole-1) */
#define MOL_H2O         18.0153 /* molecular weight of water vapor (kg / kmole) */
#define	TB 288.0
#define	STDLR	-0.0065
#define RGAS            8.31432e3 /* gas constant (J / kmole / deg) */
#define AH		1.0	/* ratio sensible/momentum phi func	*/
#define AV		1.0	/* ratio latent/momentum phi func	*/
#define CP_AIR 1.005e3 /* specific heat of air at constant pressure (J / kg / deg) */
#define DALR ( GRAV / CP_AIR ) /*  dry adiabatic lapse rate (deg / m) */
#define FREEZE 2.7316e2 /* freezing point of water at standard pressure (deg K) */
#define BOIL   3.7315e2 /* boiling point of water at standard pressure (deg K) */
#define VON_KARMAN 3.5e-1 /* Von Karman constant */
#define PAESCHKE	7.35	/* Paeschke's const	*/
#define THRESH		1.e-5	/* convergence threshold		*/
#define ITMAX		50	/* max # iterations allowed		*/
#define SECONDS_PER_DAY		(60.0f * 60.0f * 24.0f)	/* SEC * MIN * HRS */
#define BETA_S		5.2
#define BETA_U		16
#define DAYS_PER_YEAR 365
#define DAYS_PER_MONTH {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}

# define JAN			0
# define FEB			1
# define MAR			2
# define APR			3
# define MAY			4
# define JUN			5
# define JUL			6
# define AUG			7
# define SEP			8
# define OCT			9
# define NOV			10
# define DEC			11
#define MONTHS 12
#define MONTH_BEGIN     	0
#define MONTH_END		1

# define NEXT_MO(_mo)			((_mo + 1) % MONTHS)
# define PREV_MO(_mo)			((_mo + DEC) % MONTHS)

#define SM		0 /* psi function code: momentum */
#define SH		1 /* psi function code: sensible heat flux */
#define SV		2 /* psi function code: latent heat flux */

# define mapssBOREAL			0
# define mapssTEMPERATE		1
# define mapssSUBTROPICAL	        2
# define mapssTROPICAL		3


/*
 *  integral of hydrostatic equation over layer with linear temperature
 *  variation
 *
 *      pb = base level pressure, tb = base level temp (K),
 *      L = lapse rate (deg/km), h = layer thickness (km),
 *      g = grav accel (m/s^2), m = molec wt (kg/kmole),
 *
 *      (the factors 1.e-3 and 1.e3 are for units conversion)
 */
#define HYSTAT(pb,tb,L,h,g,m)           ((pb) * (((L)==0.) ?\
			exp(-(g)*(m)*(h)*1.e3/(RGAS*(tb))) :\
			pow((tb)/((tb)+(L)*(h)),(g)*(m)/(RGAS*(L)*1.e-3))))
			/*
			 *  specific humidity, as function of e & P
			 *      e = vapor pressure, P = pressure (same units)
			 */
#define SPEC_HUM(e,P)           ((e)*MOL_H2O/(MOL_AIR*(P)+(e)*(MOL_H2O-MOL_AIR)))

			/*
			 *  give density of a gas (kg/m^3) as a function of mol. wt, pressure, and temperature
			 *
			 *      p = pressure (Pa), m = molecular weight (kg/kmole),
			 *      t = temperature (K), rho = density (kg/m^3)
			 */
#define GAS_DEN(p,m,t)          ((p)*(m)/(RGAS*(t)))

			/*
			 *  virtual temperature, i.e. the fictitious temperature that air must have at
			 *  pressure p to have the same density, as a water vapor-air mixture at
			 *  pressure P, temperature t, and vapor pressure e
			 *
			 *      t = temperature (K),
			 *      e = vapor pressure, P = pressure (e and P in same units),
			 */
#define VIRT(t,e,P)             ((t)/(1.-(1.-MOL_H2O/MOL_AIR)*((e)/(P))))

			/*
			 *  latent heat of vaporization at temperature t (deg K)
			 */
#define LH_VAP(t)               (2.5e6 - 2.95573e3 *((t) - FREEZE))

			/*
			 *  latent heat of fusion at temperature t (deg K)
			 */
#define LH_FUS(t)               (3.336e5 + 1.6667e2 * (FREEZE - (t)))


			typedef enum {UNKNOWNzone = 0, ARCTICzone, BOREALzone, TEMPERATEzone, SUBTROPICALzone, TROPICALzone} ClimateZone;


			class ScienceFcns
{ // functions which communicate only thru the calling sequence
	public:
		ScienceFcns();
		~ScienceFcns();

		float efold(float tau, float raw, float previous);
		void FindSlopeYIntercept(float X1, float Y1, float X2, float Y2, float * slopeP, float * y_interceptP); // Find the slope and y-intercept 
		// of a line passing between 2 points.
		float C3prodPct(float meantmp[]); // Estimate C3 production as a % of total production, from monthly air temps, 
		// without weighting by month-to-month LAIs    
		float normalizedC3prod(float soil_tmp);
		float normalizedC4prod(float soil_tmp);
		float estimate_max_soil_tmp(float meantmp[MONTHS]);
		float DewpointTemperature(float satvp);
		float FindRatio(float x, float slope, float y_intercept);
		float GrowingDegreeDays(float tmp[], float base_temp);
		float GrowingDegreeDays(float tmp[], float base_temp, bool SouthernHemisphereFlag);
		double sati(double tk); // saturated vapor pressure over ice
		double satw(double tk); // saturated vapor pressure over water
		int hle1(
				double	press,	/* in: air pressure (Pa)			*/
				double	ta,	/* in: air temperature (K) at height za	*/
				double	ts,	/* in: surface temperature (K)		*/
				double	za,	/* in: height of air temp measurement (m)	*/
				double	ea,	/* in: vapor pressure (Pa) at height zq	*/
				double	es,	/* in: vapor pressure (Pa) at surface	*/
				double	zq,	/* in: height of spec hum measurement (m)	*/
				double	u,	/* in: wind speed (m/s) at height zu	*/
				double	zu,	/* in: height of wind speed measurement (m)	*/
				double	z0,	/* in: roughness length (m)			*/
				double *h,	/* out: sens heat flux (+ to surf) (W/m^2)	*/
				double *le,	/* out: latent heat flux (+ to surf) (W/m^2)	*/
				double *e);	/* out: mass flux (+ to surf) (kg/m^2/s)	*/
		double trbxfr(double z, double z0, double t, double t0, double e, double e0, double u, double elev);
		double satvp(double tmp);
		double psi(
				double zeta, /* z/lo				*/
				int code); /* code identifies which psi function */
		float annual_max(float monthly_vals[]);
		float annual_min(float monthly_vals[]);
		float annual_average(float monthly_vals[], bool shift_flag);
		float annual_sum(float monthly_vals[], bool shift_flag);
		bool close_enough(double a, double b, double tolerance);
		bool close_enough(double a, double b, double tol_ratio, double tol_diff);
		void MixIndex(float ppt[], float tmp[], float p_hi_mult, float * tmp_indexP, float * ppt_indexP, float * mix_indexP);
		float FrostIndex(float tmp[]);
		int ClimateZone4MAPSS(float min_tmp, float n_decid_bound, float s_decid_bound, float frost);
		ClimateZone ClimateZone4Biogeography(float gdd_zero, float min_tmp, float az_thres, float bz_thres, float tz_thres, float stz_thres);
		void PhytomorphologicalIndices(float tmp_index, float ppt_index, float * evergP, float * needlP);
		void allometry_DavidKing(float abovegr_live_wood, Boolean db_flag, float * tree_htP, float * dbhP);
		void tree_dim(float tree_lai, float ltree, bool deciduous, float * tree_htP, float * dbhP);
		void live_wood(float dbh, bool deciduous, float lwood, float * lbranchP, float * lstemP);
		float bed_depth(float tot_fuel_bed_bio, float depth_ratio);
		int doy_of_dom0(int tgt_month, bool southernHemiFlag);
		float clip(float val, float minval, float maxval);
		float crown_kill(float ht, float cl, float hk);

		// unit conversions
		float CtoF(float degC) { return(1.8f*degC + 32.0f); }
		float FtoC(float degF) { return((degF - 32.f)*5.f/9.f); }
		float mm_to_in(float mm) { return(mm*0.0393700787f); }
		float in_to_mm(float inches) { return(inches/0.0393700787f); }
		float ft_to_m(float ft) { return(ft*0.3048f); }
		float m_to_ft(float m) { return(m*3.28084f); }
		float btu_per_sec_per_ft_to_kW_per_m(float btu_per_sec_per_ft); 
		float g_per_m2_to_tons_per_acre(float g_per_m2) { return(g_per_m2*.0044609f); } 
		float lbs_per_ac_to_g_per_m2(float lbs_per_ac) { return(lbs_per_ac*0.112); }
		float cm_to_in(float cm) { return(cm*0.393700787); }
		float kJ_per_day_to_W(float kJ_per_day) { return(kJ_per_day/86.4f); }
		float m_per_sec_to_mph(float m_per_sec) { return(2.23694f*m_per_sec); }
		float degC_to_degK(float degC) { return(degC + 273.15); }
		float degK_to_degC(float degK) { return(degK - 273.15); }

		int days_per_mo[MONTHS];  
		int month_days[MONTHS][2];

}; // end of class ScienceFcns

#endif /* SCIENCEFCNS_H */
