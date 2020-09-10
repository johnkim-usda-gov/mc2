/*
 *  commandFile.cpp
 *  MC2
 */

#include <iostream>
#include <vector>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <assert.h> 
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include "netcdf.h"
#include "commandFile.h"
#include "ScienceFcns.h"
#include "category_bgc.h"
#include "MAPSSvegClasses.h"
#include "ProcessModel.h"
#include "MAPSSbiogeographyModel.h"
#include "MC2.h"

#define streq(s1, s2)		(!strcmp(s1,s2))
#define isminus(_ch)		(_ch == '-')
#define isdecimal(_ch)		(_ch == '.')


void getRunData(int argc, char * argv[]);
void makeUpperCase(char s[]);
void read_run_param(int num_valid_occur[], char * name, char * buf);
void read_run_specs(int num_valid_occur[], int argc, char * argv[]);
char * tokncpy(char * destinationP, char * sourceP, size_t maxsize, char const * delimitersP);


/* External variables local to this source file */
static char const commaOrWhitespace[]=", \t\n";

// typedef enum {unknownBaseCalibration, mc2GLOBAL, mc2ConUS, mc2ConUS_LC, mc2WCR, mc2California} BaseCalibrationEnum;
const char * BaseCalS[] = {"unknownBaseCalibration", "mc2GLOBAL", "mc2ConUS", "mc2ConUS_LC", "mc2WCR", "mc2California"}; 


const ParamStruct RunParamsList[] = 
{  						/* alphabetical order */
  {"climate_data_directory", rpString},
  {"CO2_file", rpString},
  {"col_begin", rpInt},
  {"col_end", rpInt},
  {"dummy_wind_switch", rpString},
  {"earth_data_directory", rpString},
  {"fire_mode_switch", rpString},
  {"fire_set_interval", rpInt},
  {"fire_set_jday", rpInt},
  {"fire_suppression_erc_threshold", rpFloat},
  {"fire_suppression_first_year", rpInt},
  {"fire_suppression_fli_threshold", rpFloat},
  {"fire_suppression_ros_threshold", rpFloat},
  {"fire_suppression_switch", rpString},
  {"first_calendar_year_of_run", rpInt},
  {"grid_name", rpString},
  {"mask_file", rpString},
  // {"model_parameter_set", rpString},
  {"monthly_output_switch", rpString},
  {"multiyr_len", rpInt},
  {"multiyr_output_switch", rpString},
  {"multiyr_start", rpInt},
  {"output_file_prefix", rpString},
  {"output_variable_file", rpString},
  {"output_variable", rpString},
  {"row_begin", rpInt},
  {"row_end", rpInt},
  {"run_mode", rpString},
  {"soil_bulk_density_file", rpString},
  {"soil_data_file", rpString},
  {"space_before_time_switch", rpString},
  {"suppressed_fire_cell_fraction", rpFloat},
  {"time_invariant_output_switch", rpString},
  {"warmstart_file", rpString},
  {"warmstart_output_switch", rpString},
  {"yearly_output_switch", rpString},
  {"years_offset_into_input_data", rpInt},
  {"years_to_run", rpInt}
};
#define NUM_OF_RUN_PARAMS (sizeof(RunParamsList)/sizeof(ParamStruct)) 

#define SWITCH_STR_LEN 4   
char on[SWITCH_STR_LEN + 1]  = "on  ";
char off[SWITCH_STR_LEN + 1] = "off ";

void Simulation::writeRunParameters(ncid_type * ncidP)
   {
   int rtnval;
   // typedef enum {unknownRunMode, MAPSS_EQ, CENTURY_EQ, SPINUP, TRANSIENT} RunModeEnum;				
   const char * RM[] = {"unknownRunMode", "MAPSS_EQ", "CENTURY_EQ", "SPINUP", "TRANSIENT"};
   
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "climate_data_directory", strlen(runParams.climate_data_directory), runParams.climate_data_directory); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "CO2_file", strlen(runParams.CO2_file), runParams.CO2_file); assert(chk_nc(rtnval));
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "col_begin", NC_INT, 1, &runParams.col_begin); assert(chk_nc(rtnval));
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "col_end", NC_INT, 1, &runParams.col_end); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "dummy_wind_switch", SWITCH_STR_LEN, runParams.dummyWindFlag ? on : off); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "earth_data_directory", strlen(runParams.earth_data_directory), runParams.earth_data_directory); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "fire_mode_switch", SWITCH_STR_LEN, runParams.fireModelFlag ? on : off); assert(chk_nc(rtnval));
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "fire_set_interval", NC_INT, 1, &runParams.fire_set_interval); assert(chk_nc(rtnval));
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "fire_set_jday", NC_INT, 1, &runParams.fire_set_jday); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "fire_suppression_erc_threshold", NC_FLOAT, 1, &runParams.fire_suppression_erc_threshold); assert(chk_nc(rtnval));
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "fire_suppression_first_year", NC_INT, 1, &runParams.fire_suppression_first_year); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "fire_suppression_fli_threshold", NC_FLOAT, 1, &runParams.fire_suppression_fli_threshold); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "fire_suppression_ros_threshold", NC_FLOAT, 1, &runParams.fire_suppression_ros_threshold); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "fire_suppression_switch", SWITCH_STR_LEN, runParams.fireSuppressionFlag ? on : off); assert(chk_nc(rtnval));
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "first_calendar_year_of_run", NC_INT, 1, &runParams.first_calendar_year_of_run); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "grid_name", strlen(runParams.grid_name), runParams.grid_name); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "mask_file", strlen(runParams.mask_file), runParams.mask_file); assert(chk_nc(rtnval));
   
   assert((unsigned int)m_baseCalibration<sizeof(BaseCalS)/sizeof(char *));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "m_baseCalibration", strlen(BaseCalS[m_baseCalibration]), BaseCalS[m_baseCalibration]); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "monthly_output_switch", SWITCH_STR_LEN, runParams.monthlyOutputFlag ? on : off); assert(chk_nc(rtnval));
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "multiyr_len", NC_INT, 1, &runParams.multiyr_len); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "multiyr_output_switch", SWITCH_STR_LEN, runParams.multiyrOutputFlag ? on : off); assert(chk_nc(rtnval));
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "multiyr_start", NC_INT, 1, &runParams.multiyr_start); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "output_file_prefix", strlen(runParams.output_file_prefix), runParams.output_file_prefix); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "output_variable_file", strlen(runParams.output_variable_file), runParams.output_variable_file); assert(chk_nc(rtnval));
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "row_begin", NC_INT, 1, &runParams.row_begin); assert(chk_nc(rtnval));
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "row_end", NC_INT, 1, &runParams.row_end); assert(chk_nc(rtnval));
   
   assert((unsigned int)runParams.runMode<sizeof(RM)/sizeof(char *));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "run_mode", strlen(RM[runParams.runMode]), RM[runParams.runMode]); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "soil_bulk_density_file", strlen(runParams.soil_bulk_density_file), runParams.soil_bulk_density_file); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "soil_data_file", strlen(runParams.soil_data_file), runParams.soil_data_file); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "suppressed_fire_cell_fraction", NC_FLOAT, 1, &runParams.suppressed_fire_cell_fraction); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "time_invariant_output_switch", SWITCH_STR_LEN, runParams.timeInvariantOutputFlag ? on : off); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "warmstart_file", strlen(runParams.warmstart_file), runParams.warmstart_file); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "warmstart_output_switch", SWITCH_STR_LEN, runParams.warmstartOutputFlag ? on : off); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "yearly_output_switch", SWITCH_STR_LEN, runParams.yearlyOutputFlag ? on : off); assert(chk_nc(rtnval));
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "years_offset_into_input_data", NC_INT, 1, &runParams.years_offset_into_input_data); assert(chk_nc(rtnval));
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "years_to_run", NC_INT, 1, &runParams.years_to_run); assert(chk_nc(rtnval));   
   
   rtnval = nc_put_att_double(ncidP->fileid, NC_GLOBAL, "latitude0", NC_DOUBLE, 1, &runParams.latitude0); assert(chk_nc(rtnval));
   rtnval = nc_put_att_double(ncidP->fileid, NC_GLOBAL, "longitude0", NC_DOUBLE, 1, &runParams.longitude0); assert(chk_nc(rtnval));
   rtnval = nc_put_att_double(ncidP->fileid, NC_GLOBAL, "cell_spacing", NC_DOUBLE, 1, &runParams.cell_spacing); assert(chk_nc(rtnval));
        
   } // end of writeRunParameters()


const ParamStruct ModelParamsList[] = 
{
  {"alt_fuel_load_calculation", rpString},
  {"alt_tree_allometry_calculation", rpString},
  {"bui_threshold", rpFloat},
  {"bui_threshold_by_vtype", rpFloat},
  {"bz_thres", rpFloat},
  {"c3_threshold", rpFloat},
  {"century_path", rpString},
  {"century_runoff_slope", rpFloat},
  {"century_runoff_x_intercept", rpFloat},
  {"code_flag", rpCodeFlag},
  {"crown_fire_mortality_pct", rpFloat},
  {"desert_grass_C_max", rpFloat},
  {"double_CO2_grass_AT_mult", rpFloat},
  {"double_CO2_grass_NPP_mult", rpFloat},
  {"double_CO2_tree_AT_mult", rpFloat},
  {"double_CO2_tree_NPP_mult", rpFloat},
  {"dry_temperate_threshold", rpFloat},
  {"efold_t_max", rpFloat},
  {"ffmc_threshold", rpFloat},
  {"ffmc_threshold_by_vtype", rpFloat},
  {"fire_min", rpFloat},
  {"forest_thres_C", rpFloat},
  {"grassfrac_thres", rpFloat},
  {"lait_lower_limit", rpFloat},
  {"MAPSS_parameter_set", rpString},
  {"MAPSS_subalpine_threshold", rpFloat},
  {"maritime_threshold", rpFloat},
  {"max_LAI", rpFloat},
  {"max_grass_NPP", rpFloat},
  {"max_tree_NPP", rpFloat},
  {"moist_temperate_threshold", rpFloat},
  {"p_hi_mult", rpFloat}, 
  {"PKLZ_lower_bound_coefficients", rpFloat},
  {"pprdwc_grass", rpFloat}, 
  {"pprdwc_tree", rpFloat}, 
  {"PSFZ_lower_bound_coefficients", rpFloat},
  {"SAFZ_lower_bound_coefficients", rpFloat},
  {"southern_hemisphere_switch", rpString},
  {"SSZ_upper_bound_coefficients", rpFloat},
  {"subalpine_threshold", rpFloat},
  {"tmmin_threshold", rpFloat},
  {"unlimited_N_switch", rpString},
  {"vtype_fire_return_interval_range", rpFloat},
  {"woodl_thres_C", rpFloat}
};
#define NUM_OF_MODEL_PARAMS (sizeof(ModelParamsList)/sizeof(ParamStruct))

#define NUM_OF_PARAMS (NUM_OF_RUN_PARAMS+NUM_OF_MODEL_PARAMS)


void Simulation::writeModelParameters(ncid_type * ncidP)
   {
   int rtnval;
   char cf[NUM_CODE_FLAGS*SWITCH_STR_LEN + 1];
   
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "alt_fuel_load_calculation", SWITCH_STR_LEN, modelParams.altFuelLoadFlag ? on : off); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "alt_tree_allometry_calculation", SWITCH_STR_LEN, modelParams.altTreeAllometryFlag ? on : off); assert(chk_nc(rtnval));
   for (int zone = 1; zone<=5; zone++)
      {
      char att_name[] = "bui_threshold_zone_x_1_to_8";
      att_name[19] = '0' + zone;
      // printf("strlen(att_name), att_name = %d, %s\n", (int)strlen(att_name), att_name);
      rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, att_name, NC_FLOAT, 8, modelParams.bui_threshold[zone] + 1); 
      assert(chk_nc(rtnval));
      };
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "bui_threshold_by_vtype_1_to_MAX_VTYPE", NC_FLOAT, MAX_VTYPE, modelParams.bui_threshold_by_vtype + 1); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "bz_thres", NC_FLOAT, 1, &modelParams.bz_thres); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "c3_threshold", NC_FLOAT, 1, &modelParams.c3_threshold); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "century_path", strlen(modelParams.century_path), modelParams.century_path); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "century_runoff_slope", NC_FLOAT, 1, &modelParams.m_century_runoff_slope); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "century_runoff_x_intercept", NC_FLOAT, 1, &modelParams.m_century_runoff_x_intercept); assert(chk_nc(rtnval));
   
   modelParams.code_flags[ALT_FUEL_LOAD_CODE_FLAG] = modelParams.altFuelLoadFlag;
   modelParams.code_flags[ALT_TREE_ALLOMETRY_CODE_FLAG] = modelParams.altTreeAllometryFlag;
   cf[0] = 0;
   for (int k = 0; k<NUM_CODE_FLAGS; k++) strcat(cf, modelParams.code_flags[k] ? on : off);
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "code_flags", strlen(cf), cf); assert(chk_nc(rtnval));
   
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "moist_temperate_threshold", NC_FLOAT, 1, &modelParams.moist_temperate_threshold); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "crown_fire_mortality_pct", NC_FLOAT, 1, &modelParams.crown_fire_mortality_pct); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "desert_grass_C_max", NC_FLOAT, 1, &modelParams.desert_grass_C_max); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "double_CO2_grass_AT_mult", NC_FLOAT, 1, &modelParams.double_CO2_grass_AT_mult); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "double_CO2_grass_NPP_mult", NC_FLOAT, 1, &modelParams.double_CO2_grass_NPP_mult); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "double_CO2_tree_AT_mult", NC_FLOAT, 1, &modelParams.double_CO2_tree_AT_mult); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "double_CO2_tree_NPP_mult", NC_FLOAT, 1, &modelParams.double_CO2_tree_NPP_mult); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "dry_temperate_threshold", NC_FLOAT, 1, &modelParams.dry_temperate_threshold); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "efold_t_max", NC_FLOAT, 1, &modelParams.efold_t_max); assert(chk_nc(rtnval));
   for (int zone = 1; zone<=5; zone++)
      {
      char att_name[] = "ffmc_threshold_zone_x_1_to_8";
      att_name[20] = '0' + zone;
      rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, att_name, NC_FLOAT, 8, modelParams.ffmc_threshold[zone] + 1); 
      assert(chk_nc(rtnval));
      };
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "ffmc_threshold_by_vtype_1_to_MAX_VTYPE", NC_FLOAT, MAX_VTYPE, modelParams.ffmc_threshold_by_vtype + 1); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "fire_min", NC_FLOAT, 1, &modelParams.fire_min); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "forest_thres_C", NC_FLOAT, 1, &modelParams.m_forest_thres_C); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "grassfrac_thres", NC_FLOAT, 1, &modelParams.grassfrac_thres); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "lait_lower_limit", NC_FLOAT, 1, &modelParams.m_lait_lower_limit); assert(chk_nc(rtnval));
   
   assert((unsigned int)modelParams.MAPSSparameterSet<sizeof(MpSN)/sizeof(char *));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "MAPSS_parameter_set", strlen(MpSN[modelParams.MAPSSparameterSet]), MpSN[modelParams.MAPSSparameterSet]); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "MAPSS_subalpine_threshold", NC_FLOAT, 1, &modelParams.MAPSS_subalpine_threshold); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "maritime_threshold", NC_FLOAT, 1, &modelParams.maritime_threshold); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "max_LAI", NC_FLOAT, 5, modelParams.max_LAI); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "max_grass_NPP", NC_FLOAT, 3, modelParams.max_grass_NPP); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "max_tree_NPP", NC_FLOAT, 5, modelParams.max_tree_NPP); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "p_hi_mult", NC_FLOAT, 1, &modelParams.p_hi_mult); assert(chk_nc(rtnval)); 
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "PKLZ_lower_bound_coefficients", NC_FLOAT, 9, modelParams.PKLZ_lb_coeffs); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "pprdwc_grass", NC_FLOAT, 3, modelParams.pprdwc_grass); assert(chk_nc(rtnval)); 
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "pprdwc_tree", NC_FLOAT, 3, modelParams.pprdwc_tree); assert(chk_nc(rtnval)); 
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "PSFZ_lower_bound_coefficients", NC_FLOAT, 9, modelParams.PSFZ_lb_coeffs); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "SAFZ_lower_bound_coefficients", NC_FLOAT, 9, modelParams.SAFZ_lb_coeffs); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "southern_hemisphere_switch", SWITCH_STR_LEN, modelParams.southernHemisphereFlag ? on : off); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "SSZ_upper_bound_coefficients", NC_FLOAT, 9, modelParams.SSZ_ub_coeffs); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "subalpine_threshold", NC_FLOAT, 1, &modelParams.subalpine_threshold); assert(chk_nc(rtnval));
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "tmmin_threshold", NC_FLOAT, 1, &modelParams.tmmin_threshold); assert(chk_nc(rtnval));
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "unlimited_N_switch", SWITCH_STR_LEN, modelParams.unlimitedNflag ? on : off); assert(chk_nc(rtnval));
   for (int vtypeNdx = 0; vtypeNdx<(MAX_VTYPE+1); vtypeNdx++)
      {
      char att_name[] = "vtypeXX_fire_return_interval_range";
      assert(vtypeNdx<=99);
      att_name[5] = '0' + vtypeNdx/10;
      att_name[6] = '0' + vtypeNdx%10;
      
      rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, att_name, NC_FLOAT, 2, modelParams.vveg2mfri[vtypeNdx] + 1); 
      assert(chk_nc(rtnval));
      };
   rtnval = nc_put_att_float(ncidP->fileid, NC_GLOBAL, "woodl_thres_C", NC_FLOAT, 1, &modelParams.m_woodl_thres_C); assert(chk_nc(rtnval));

   } // end of writeModelParameters()


void Simulation::writeCommands(ncid_type * ncidP, char * cmd_line, char * file_name_and_path)
   {
   int rtnval;
   char buf[STRLEN+1];
   FILE * input;
   int line_num;
   char cmdline_att_name[] = "cmdfile_linexxx";
   char line_num_str[6];
   int buflen;
   bool end_flag;
   
   rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, "cmdline", strlen(cmd_line), cmd_line); assert(chk_nc(rtnval));
   
   access_file(&input, file_name_and_path, "r");
   line_num = 1;
   end_flag = false;
   while (!feof(input) && !end_flag)
      {
      sprintf(line_num_str, "%d", line_num);
      while (line_num_str[0]==' ') 
         { 
         int len = strlen(line_num_str);
         for (int i = 0; i<len; i++) line_num_str[i] = line_num_str[i+1];
         }
      cmdline_att_name[strlen("cmdfile_line")] = 0;
      strcat(cmdline_att_name, line_num_str);
      
      buf[0] = 0;
      fgets(buf, sizeof(buf), input);    
      buflen = strlen(buf);
      if (buf[buflen - 1]=='\n') buf[buflen - 1] = 0;
      assert(buflen<=STRLEN);
      
      rtnval = nc_put_att_text(ncidP->fileid, NC_GLOBAL, cmdline_att_name, strlen(buf), buf); assert(chk_nc(rtnval)); 
      line_num++;  
      
      // end_flag = strcmp(buf, "end_of_commands")==0;
      end_flag = line_num>=m_cmd_file_len;
      }
   fclose(input);
   
   } // end of writeCommands()


bool Simulation::getRunParameters(char * command_file_name_and_path)
  {
  int num_valid_occur[NUM_OF_PARAMS];
  const char * charP;
  
  runParams.baseCalibration = m_baseCalibration;

  printf("\n\n*** Setting default values of parameters for calibration %s now.\n", BaseCalS[m_baseCalibration]);
  for (unsigned int i=0; i<NUM_OF_RUN_PARAMS; i++) if (RunParamsList[i].fcnType!=rpUnknown)
    { /* Set default values. */
    printf("default ");
    switch (RunParamsList[i].fcnType)
      {
      case rpString: charP = rpStringFcn(RunParamsList[i].keyword, 0, NULL); break;
      case rpInt: charP = rpIntFcn(RunParamsList[i].keyword, 0, NULL); break;
      case rpFloat: charP = rpFloatFcn(RunParamsList[i].keyword, 0, NULL); break;
      case rpCodeFlag: charP = rpCodeFlagFcn(RunParamsList[i].keyword, 0, NULL); break;
      default: assert(0); break;
      } // end of switch on fcnType
    printf("\n");
    // if (charP!=NULL) { printf("\n\n%s\n", charP); err_exit("getRunParameters internal error #1"); }
    } // end of loop thru all the run parameters
  for (unsigned int i=0; i<NUM_OF_MODEL_PARAMS; i++) if (ModelParamsList[i].fcnType!=rpUnknown)
    { /* Set default values. */
    printf("default ");
    switch (ModelParamsList[i].fcnType)
      {
      case rpString: charP = rpStringFcn(ModelParamsList[i].keyword, 0, NULL); break;
      case rpInt: charP = rpIntFcn(ModelParamsList[i].keyword, 0, NULL); break;
      case rpFloat: charP = rpFloatFcn(ModelParamsList[i].keyword, 0, NULL); break;
      case rpCodeFlag: charP = rpCodeFlagFcn(ModelParamsList[i].keyword, 0, NULL); break;
      default: assert(0); break;
      } // end of switch on fcnType
    printf("\n");
    // if (charP!=NULL) { printf("\n\n%s\n", charP); err_exit("getRunParameters internal error #2"); }
    } // end of loop thru all the model parameters
  
   // Initialize the discretionary output variable dictionary and output lists.
   printf("\nLAST_MONTHLY_VAR, LAST_YEARLY_VAR, LAST_MULTIYR_VAR, LAST_SINGLE_VAR, NUM_OPTIONAL_OUTVARS ="
         " %d, %d, %d, %d, %d\n", 
         LAST_MONTHLY_VAR, LAST_YEARLY_VAR, LAST_MULTIYR_VAR, LAST_SINGLE_VAR, NUM_OPTIONAL_OUTVARS);
   initializeVarDict();
   for (int index = 0; index<NUM_OPTIONAL_OUTVARS; index++)
         m_moOutList[index].show = m_yrOutList[index].show = m_multiyrOutList[index].show  = false;             
   m_num_monthly_vars = 0;
   m_num_yearly_vars = 0;
   m_num_multiyr_vars = 0;
   m_num_single_vars = 0;
   
   // Load up a few default output variables.
   addVarToYrList(VTYPE);
   addVarToMultiyrList(VTYPE);
   addVarToYrList(PHYSIOGNOMIC_CLASS);
   addVarToMultiyrList(PHYSIOGNOMIC_CLASS);
   addVarToYrList(BIOME);
   addVarToMultiyrList(BIOME);
   addVarToYrList(C_VEG);
   addVarToMultiyrList(C_VEG);
   addVarToMoList(NPP);
   addVarToYrList(NPP);
   addVarToMultiyrList(NPP);
   addVarToYrList(PART_BURN);
   addVarToMultiyrList(PART_BURN);

  /* Now get the specific run parameters for this simulation. */
  printf("\n\n*** Getting specific values of run parameters for this simulation now. ***");
  for (unsigned int rpi=0; rpi<NUM_OF_PARAMS; rpi++) num_valid_occur[rpi] = 0; 
  m_cmd_file_len = process_command_file(command_file_name_and_path, num_valid_occur);     
  
  printf("\n");

  runParams.runMode = runParams.interpretRunMode(runParams.run_mode);
  runParams.gridName = runParams.interpretGridName(runParams.grid_name);
  runParams.latlonFlag = (0.<runParams.cell_spacing && runParams.cell_spacing<=180.)
	    && (-90.<runParams.latitude0 && runParams.latitude0<=90.)
      && (-180<=runParams.longitude0 && runParams.longitude0<180.);
  runParams.singleCellFlag = runParams.diags = (runParams.col_begin==runParams.col_end) && (runParams.row_begin==runParams.row_end);

  return(true);
  } // end of Simulation::getRunParameters()
  
  
bool Simulation::addVarToList(char * name)
   {
   /* look through VarDict for match */
   bool found = FALSE;
   int i = 0;
   while (!found && i<NUM_OPTIONAL_OUTVARS) 
      {
      found = strcmp(name, VarDict[i].name)==0;
      if (!found) i++;
      }
      
   if (found)
      { // Add the variable to the appropriate output file lists.
      VarDict[i].show = true; VarDict[i].m_pS = this;
      switch (VarDict[i].interval)
         {
         case MONTH_INTERVAL: addVarToMoList(i);
         case YEAR_INTERVAL: addVarToYrList(i);
         case MULTIYR_INTERVAL: addVarToMultiyrList(i); 
            break;
         case SINGLE_INTERVAL: addVarToSingleList(i);
            break;
         default: assert(0); break;
         }
      }
   else printf("*** Could not find %s in the variable dictionary.\n", name);
   
   return(found);
   } // end of Simulation::addVarToList()
   

void Simulation::addVarToSingleList(int varNdx)
   {
   if (VarDict[varNdx].singleOutNdx>=0) return;
   VarDict[varNdx].singleOutNdx = m_num_single_vars; 
   VarDict[varNdx].show = true; VarDict[varNdx].m_pS = this;
   m_singleOutList[m_num_single_vars] = VarDict[varNdx];
   m_num_single_vars++;
   } // end of addVarToSingleList()
   
   
void Simulation::addVarToMoList(int varNdx)
   {
   if (VarDict[varNdx].moOutNdx>=0) return;
   VarDict[varNdx].moOutNdx = m_num_monthly_vars; 
   VarDict[varNdx].show = true; VarDict[varNdx].m_pS = this;
   m_moOutList[m_num_monthly_vars] = VarDict[varNdx];
   m_num_monthly_vars++;
   } // end of addVarToMoList()
   
   
void Simulation::addVarToYrList(int varNdx)
   {
   if (VarDict[varNdx].yrOutNdx>=0) return;
   VarDict[varNdx].yrOutNdx = m_num_yearly_vars; 
   VarDict[varNdx].show = true; VarDict[varNdx].m_pS = this;
   m_yrOutList[m_num_yearly_vars] = VarDict[varNdx];
   m_num_yearly_vars++;
   } // end of addVarToYrList()
   
   
void Simulation::addVarToMultiyrList(int varNdx)
   {
   if (VarDict[varNdx].multiyrOutNdx>=0) return;
   VarDict[varNdx].multiyrOutNdx = m_num_multiyr_vars; 
   VarDict[varNdx].show = true; VarDict[varNdx].m_pS = this;
   m_multiyrOutList[m_num_multiyr_vars] = VarDict[varNdx];
   m_num_multiyr_vars++;
   } // end of addVarToMultiyrList()



int Simulation::process_command_file(char * file_name_and_path, int num_valid_occur[])
/* Read values specified in command file. */	
   {
   FILE * input;
   const char equalsOrWhitespace[]="= \t";
   char buf[STRLEN+1], name[STRLEN+1], value[STRLEN+1];
   char * charP;
   int buflen, namelen;
   bool skip_flag;
   int cmd_file_len = 0;
	
   fprintf(stdout, "\n");
   
   access_file(&input, file_name_and_path, "r");
   while (!feof(input))
      {
      cmd_file_len++;
      buf[0] = 0;
      skip_flag = FALSE;
      
      fgets(buf, sizeof(buf), input );
      
      buflen = strlen(buf);
      skip_flag = buflen<=0 || ispunct(buf[0]) || isspace(buf[0]);
      if (!skip_flag)
         {
         /* Copy parameter name into name[]. */
         charP = strtok(buf, equalsOrWhitespace);
         strncpy(name, charP, sizeof(name)-1); name[sizeof(name)-1] = 0;
         if ((strlen(name)==(strlen("end_of_commands")+1) || strlen(name)==(strlen("end_of_commands")+2))
                && strncmp(name, "end_of_commands", strlen("end_of_commands"))==0) 
                break;        
         /* Copy remainder of line or argument string into value[]. */
         namelen = strlen(name);
         if (namelen<buflen)
            {
            charP += namelen + 1;
            charP += strspn(charP, equalsOrWhitespace);
            strncpy(value, charP, sizeof(value)-1); value[sizeof(value)-1] = 0;
            }
         else value[0] = 0;
         
         read_param(num_valid_occur, name, value);
         }
         
      }
   fclose(input);
   
   return(cmd_file_len);
   } /* end of read_command_file() */


void Simulation::read_param(int num_valid_occur[], char * name, char * buf)
{
  const char * keyword;
  const char * msg;
   unsigned int i;
  int j;
  
  for (i=0; i<(NUM_OF_RUN_PARAMS + NUM_OF_MODEL_PARAMS); i++)
  {  
    j = i - NUM_OF_RUN_PARAMS;
    if ((i<NUM_OF_RUN_PARAMS ? RunParamsList[i].fcnType : ModelParamsList[j].fcnType)==rpUnknown) continue;
    keyword = i<NUM_OF_RUN_PARAMS ? RunParamsList[i].keyword : ModelParamsList[j].keyword;
    if (strcmp(keyword, name)==0) 
    {
      printf("for this simulation, ");
      switch (i<NUM_OF_RUN_PARAMS ? RunParamsList[i].fcnType : ModelParamsList[j].fcnType)
      {
        case rpString: msg = rpStringFcn(keyword, num_valid_occur[i], buf); break;
        case rpInt: msg = rpIntFcn(keyword, num_valid_occur[i], buf); break;
        case rpFloat: msg = rpFloatFcn(keyword, num_valid_occur[i], buf); break;
        case rpCodeFlag: msg = rpCodeFlagFcn(keyword, num_valid_occur[i], buf); break;
        default: assert(0); break;
      }
      printf("\n");
      if (msg==NULL) num_valid_occur[i]++;
      else 
      {
        printf("\n%s incorrectly specified.\n%s\n", keyword, msg);
        err_exit("read_param(): Couldn't interpret parameter value");
      }
      break;
    }
  }
  if (i>=(NUM_OF_RUN_PARAMS + NUM_OF_MODEL_PARAMS)) 
  {
    printf("\n\n*****************************************"
        "\nUnrecognized line in command file:"
		    "\n%s %s", name, buf);
    err_exit("read_param: Didn't recognize parameter name");
  }
} /* end of read_param() */


/* rpFcn - function for processing a run parameter
   char * rpFcn(keyword, num_valid_occur, input_buf)
   if (input_buf==NULL) put default value(s) for this parameter into RunData
   else 
   {
     try to interpret input_buf as value(s) of the parameter
     if unsuccessful, return with address of descriptive message
     put the values into RunData
   }
   print the values and the keyword
   return NULL
*/   


/*
typedef struct
{
  char * name;
  int max_num;
  double default_value;
  float lower_limit;
  float upper_limit;
  Boolean inclusive_lower_limit;
  Boolean inclusive_upper_limit;
  char * msg;  
} FloatParamInfoStruct;
*/

const char * Simulation::rpFloatFcn(const char * keyword, int times_called, char * buf)
{
  static const FloatParamInfoStruct info[] = 
  {
/* 0 */ {"fire_suppression_erc_threshold", 1, 60.0, 0.0, 1000., TRUE, TRUE,
        "fire_suppression_erc_threshold <ERC value> default = 60"},
/* 1 */ {"fire_suppression_fli_threshold", 1, 900.0, 0.0, 10000.0, TRUE, TRUE,
        "fire_suppression_fli_threshold <FLI value> default = 900"},
/* 2 */ {"fire_suppression_ros_threshold", 1, 100.0, 0.0, 1000.0, TRUE, TRUE, 
        "fire_suppression_ros_threshold <ROS value> default = 100"},
/* 3 */ {"suppressed_fire_cell_fraction", 1, 0.0006, 0.0, 1.0, TRUE, TRUE,
        "suppressed_fire_cell_fraction <fraction> default = 0.0006"},
/* 4 */ {"spare4", 1, 0., 0., 0., FALSE, FALSE, "spare"},
/* 5 */ {"spare5", 1, 0., 0., 0., FALSE, FALSE, "spare"},
/* 6 */ {"century_runoff_slope", 1, 0.15, 0.0, 1.0, TRUE, TRUE, 
        "century_runoff_slope <century runoff slope> default = 0.15"},
/* 7 */ {"century_runoff_x_intercept", 1, 8.0, 0.0, 100.0, TRUE, TRUE, 
        "century_runoff_x_intercept <century runoff X intercept> default = 8 cm/mo "},
/* 8 */ {"efold_t_max", 1, 10.0, 1.0, 100.0, FALSE, TRUE,
        "efold_t_max <max Efold value> default = 10.0"},
/* 9 */ {"spare9", 1, 0., 0., 0., FALSE, FALSE, "spare"},
/*10 */ {"spare10", 1, 0., 0., 0., FALSE, FALSE, "spare"},
/*11 */ {"fire_min", 1, -9999., 0.0, 1.0, TRUE, TRUE,
        "fire_min <minimum fire effect fraction> default = 0.0 (VINCERA), 0.33 (otherwise)  minimum fraction of live C pools killed by fire; "
        "applies only to the part of the cell affected by the fire"},
/*12 */ {"lait_lower_limit", 1, 0.2, 0.0, 100.0, TRUE, TRUE,
        "lait_lower_limit     <lait lower limit> default = 0.2"},
/*13 */ {"maritime_threshold", 1, -9999., 10.0, 25.0, TRUE, TRUE,
        "maritime_threshold <upper limit of efolded continentality index for maritime EN forests> default depends on biogeography option"},
/*14 */ {"p_hi_mult", 1, 1.0, 0.0, 10.0, TRUE, TRUE,
        "p_hi_mult <multiplier for ppt_avg in calculation of p_hi in lifeform_rules()> default = 1.0"},
/*15 */ {"subalpine_threshold", 1, -9999., 0.0, 10000.0, TRUE, TRUE,
        "subalpine_threshold <upper GDD (ref. 0 C) limit for subalpine zone> default depends on biogeography option"},
/*16 */ {"tmmin_threshold", 1, 0, -50.0, 50.0, FALSE, FALSE, 
        "tmmin_threshold <threshold for mean temp of coldest month for distinguishing maritime needleleaf from cool needleleaf forests> default = 0 C"},
/*17 */ {"forest_thres_C", 1, 3000., 0.0, 1e7, TRUE, FALSE, 
        "forest_thres_C <lower limit of live tree carbon for forests, gC m-2> default depends on biogeography option"},
/*18 */ {"woodl_thres_C", 1, 1150., 0.0, 1e6, TRUE, FALSE, 
        "woodl_thres_C <lower limit of live tree carbon for woodlands and savannas, gC m-2> default depends on biogeography option"},
/*19 */ {"MAPSS_subalpine_threshold", 1, 1900., -1., 10000.0, TRUE, TRUE,
        "MAPSS_subalpine_threshold <upper GDD (ref. 0 C) limit for subalpine forests and tree savanna> default = 1900., but may be affected by MAPSS parameter set option"},
/*20 */ {"desert_grass_C_max", 1, 385., 0.0, 1e6, FALSE, FALSE,
         "desert_grass_C_max <upper limit of grass C for deserts and semi-deserts,"
         " g C m-2> default = 385. g C m-2"},
/*21 */ {"c3_threshold", 1, 55.0, 0.0, 100.0, TRUE, TRUE,
         "c3_threshold <C3% of production> default = 55%; threshold for distinguishing C3 veg types from C4 veg types"},
/*22 */ {"grassfrac_thres", 1, 0.63, 0.0, 100., TRUE, TRUE,
         "grassfrac_thres <min grass fraction for grassland and savanna> default = 0.63; threshold for distinguishing shrubland from grassland/savanna"},
/*23 */ {"crown_fire_mortality_pct", 1, 98., 0., 100., TRUE, TRUE,
         "crown_fire_mortality_pct <mortality %> default = 98%"},
/*24 */ {"bz_thres", 1, -15., -50.0, 50.0, TRUE, TRUE, 
         "bz_thres <threshold for mean temp of coldest month for distinguishing boreal zone "
         "from temperate zone, deg C> default = -15. C"}, 
/*25 */ {"double_CO2_tree_AT_mult", 1, NC_FILL_FLOAT, 0., 1., FALSE, TRUE,
         "double_CO2_tree_AT_mult < multiplier between 0 and 1 > effect of doubled CO2 on transpiration rate of woody plants (orig. CO2ITR in tree.100)"},
/*26 */ {"double_CO2_tree_NPP_mult", 1, NC_FILL_FLOAT, 1., 2., TRUE, TRUE,
         "double_CO2_tree_NPP_mult < multiplier between 1. and 2. > effect of doubled CO2 on NPP of woody plants (orig. CO2IPR in tree.100)"},
/*27 */ {"max_LAI", 5, NC_FILL_FLOAT, 0., 100., FALSE, FALSE,
         "max_LAI < maximum leaf area index for SUPRT	-DN-	-EN-	-DB-	-EB-, projected leaf area m2 m-2 > (orig. MAXLAI in tree.100)"},
/*28 */ {"max_tree_NPP", 5, NC_FILL_FLOAT, 0., 1000., FALSE, FALSE,
         "max_tree_NPP < maximum monthly NPP for SUPRT	-DN-	-EN-	-DB-	-EB-, g C m-2 month-1 > (orig. PRDX(4) in tree.100)"},
/*29 */ {"double_CO2_grass_AT_mult", 1, NC_FILL_FLOAT, 0., 1., FALSE, TRUE,
         "double_CO2_grass_AT_mult < multiplier between 0 and 1 > effect of doubled CO2 on transpiration rate of herbaceous plants (orig. CO2ITR in crop.100)"},
/*30 */ {"double_CO2_grass_NPP_mult", 1, NC_FILL_FLOAT, 1., 2., TRUE, TRUE,
         "double_CO2_grass_NPP_mult < multiplier between 1. and 2. > effect of doubled CO2 on NPP of herbaceous plants (orig. CO2IPR in crop.100)"},
/*31 */ {"max_grass_NPP", 3, NC_FILL_FLOAT, 0., 1000., FALSE, FALSE,
         "max_grass_NPP < maximum monthly NPP for SUPRG	-C3-	-C4-, g C m-2 month-1 > (orig. PRDX(1) in crop.100)"},
/*32 */ {"pprdwc_grass", 3, NC_FILL_FLOAT, 0., 10., TRUE, TRUE,
         "pprdwc_grass < A, B, C parameters to pprdwc() for herbaceous plants > determines how production is limited by water availability relative to PET"},
/*33 */ {"pprdwc_tree", 3, NC_FILL_FLOAT, 0., 10., TRUE, TRUE,
         "pprdwc_tree < A, B, C parameters to pprdwc() for woody plants > determines how production is limited by water availability relative to PET"},
/*34 */ {"vtype_fire_return_interval_range", 3, NC_FILL_FLOAT, 1., 9999., TRUE, TRUE,
         "vtype_fire_return_interval_range < vtype minFRI maxFRI > fire return interval range for one potential vegetation type"},
/*35 */ {"SSZ_upper_bound_coefficients", 9, NC_FILL_FLOAT, -10000., 10000., FALSE, FALSE,
         "intercept and slope for Sitka spruce zone upper bound calculation"},
/*36 */ {"PSFZ_lower_bound_coefficients", 9, NC_FILL_FLOAT, -10000., 10000., FALSE, FALSE,
         "intercept and slope for Pacific silver fir zone lower bound calculation"},
/*37 */ {"SAFZ_lower_bound_coefficients", 9, NC_FILL_FLOAT, -10000., 10000., FALSE, FALSE,
         "intercept and slope for subalpine fir zone lower bound calculation"},
/*38 */ {"PKLZ_lower_bound_coefficients", 9, NC_FILL_FLOAT, -10000., 10000., FALSE, FALSE,
         "intercept and slope for subalpine parkland zone lower bound calculation"},
/*39 */ {"ffmc_threshold", 3, -9999., 0.0, 1000.0, TRUE, TRUE,
        "ffmc_threshold <climate_zone tree_type FFMC_fire_threshold>"},
/*40 */ {"bui_threshold", 3, -9999., 0.0, 1000.0, TRUE, TRUE,
        "bui_threshold <climate_zone tree_type BUI_fire_threshold>"},
        // {"spare", 1, 0., 0., 0., FALSE, FALSE, "spare"},
/*41 */ {"moist_temperate_threshold", 1, 635., 0.0, 10000.0, TRUE, TRUE,
        "moist_temperate_threshold <threshold between moist temperate needleleaf forests and mesic temperate needleleaf forests> default 635 mm H2O yr-1 (= 25 in H2O yr-1)"},
/*42 */ {"dry_temperate_threshold", 1, 432., 0.0, 10000.0, TRUE, TRUE,
        "dry_temperate_threshold <threshold between dry temperate needleleaf forests and mesic temperate needleleaf forests> default 432 mm H2O yr-1 (= 17 in H2O yr-1)"},
/*43 */ {"ffmc_threshold_by_vtype", 2, -9999., -1., 1000.0, TRUE, TRUE,
        "ffmc_threshold_by_vtype <vtype FFMC_fire_threshold>"},
/*44 */ {"bui_threshold_by_vtype", 2, -9999., -1., 1000.0, TRUE, TRUE,
        "bui_threshold_by_vtype <vtype BUI_fire_threshold>"},
  };	
  int num_of_params = sizeof(info)/sizeof(FloatParamInfoStruct);
  float * param_valP;
  int ip; /* index to info[] */
  int jval; /* index to particular value: param_val */
  int nval; /* number of values */
  float * valP; /* pointer to destination in MP or RunData */
  char valbuf[STRLEN+1];
  int nread;
  char * val_charP;
  int rtnval;
  extern CenParamStruct treeParams[];
  extern CenParamStruct cropParams[];
  extern CenParamStruct pprdwcParams[];
  
  ip = 0;
  while (ip<num_of_params && strcmp(keyword, info[ip].name)!=0) ip++;
  assert(ip<num_of_params);
  if (ip>=num_of_params) return "Oops! rpFloat() internal error #1. <<<<<<<<<<<<<<<<";

  nval = info[ip].max_num;
  param_valP = (float *)malloc(nval*sizeof(double));
  if (param_valP==NULL) { printf("*** Error allocating for param_val array in rpFloat().\n"); err_exit("MALLOC_FAILED"); }
  switch (ip)
  { // set default value(s)
    case 11: /* fire_min */
      switch (m_baseCalibration) 
      {
        default: *(param_valP) = 0.33; break;
        // case mc2VINCERA: *(param_valP) = 0.0; break;
      }
      break;
    case 13: /* maritime_threshold */
      switch (m_baseCalibration) 
      {
        default: *(param_valP) = 18.; break;
        case mc2WCR:
        case mc2ConUS: *(param_valP) = 16.; break;
        case mc2ConUS_LC: *(param_valP) = 16.; break;
        case mc2California: *(param_valP) = 16.; break;
        // case mc2CA08: case mc2YOSE: *(param_valP) = 14.; break;
        // case mc2VINCERA: *(param_valP) = 16.; break;
      }
      break;
    case 15: /* subalpine_threshold */
      switch (m_baseCalibration) 
      {
        default: *(param_valP) = 1900.; break;
        // case mc2NA8km: *(param_valP) = 1600.; break;
        // case mc2CA08: *(param_valP) = 2000.; break;
      }
      break;
    case 16: /* tmmin_threshold */
      switch (m_baseCalibration)
      {
         default: *(param_valP) = 0.0; break;
         case mc2WCR:
         case mc2ConUS: *(param_valP) = 1.5; break;
         case mc2ConUS_LC: *(param_valP) = 1.5; break;
         case mc2California: *(param_valP) = 1.5; break;
      }
      break;
    case 17: /* forest_thres_C */
      switch (m_baseCalibration) 
      {
        default: *(param_valP) = 3000.; break;
        // case mc2VINCERA: *(param_valP) = 3500.; break;
      }
      break;
    case 18: /* woodl_thres_C */
      switch (m_baseCalibration) 
      {
        default: *(param_valP) = 1150.; break;
        // case mc2VINCERA: *(param_valP) = 1250.; break;
        // case mc2US50km: *(param_valP) = 1450.; break;
      }
      break;
   case 25: /* double_CO2_tree_AT_mult */
      rtnval = getCenturyParameter(treeParams, "CO2ITR", 0, param_valP); assert(rtnval==1);
      break;
   case 26: /* double_CO2_tree_NPP_mult */
      rtnval = getCenturyParameter(treeParams, "CO2IPR", 0, param_valP); assert(rtnval==1);
      break;
   case 27: /* max_LAI[5] */
      for (jval = 0; jval<nval; jval++) 
         {
         rtnval = getCenturyParameter(treeParams, "MAXLAI", jval, param_valP + jval); assert(rtnval==1);
         }
      break;
   case 28: /* max_tree_NPP[5] */
      for (jval = 0; jval<nval; jval++) 
         {
         rtnval = getCenturyParameter(treeParams, "PRDX(4)", jval, param_valP + jval); assert(rtnval==1);
         }
      break;
   case 29: /* double_CO2_grass_AT_mult */
      rtnval = getCenturyParameter(cropParams, "CO2ITR", 0, param_valP); assert(rtnval==1);
      break;
   case 30: /* double_CO2_grass_NPP_mult */
      rtnval = getCenturyParameter(cropParams, "CO2IPR", 0, param_valP); assert(rtnval==1);
      break;
   case 31: /* max_grass_NPP[3] */
      for (jval = 0; jval<nval; jval++) 
         {
         rtnval = getCenturyParameter(cropParams, "PRDX(1)", jval, param_valP + jval); assert(rtnval==1);
         }
      break;
   case 32: /* pprdwc_grass[3] originally 0.0, 1.0, 0.6 */
      for (jval = 0; jval<nval; jval++) 
         {
         rtnval = getCenturyParameter(pprdwcParams, "pprdwc_grass", jval, param_valP + jval); assert(rtnval==1);
         }
      break;
   case 33: /* pprdwc_tree[3] originally 0.5, 1.0, 0.9 */
      for (jval = 0; jval<nval; jval++) 
         {
         rtnval = getCenturyParameter(pprdwcParams, "pprdwc_tree", jval, param_valP + jval); assert(rtnval==1);
         }
      break;
   case 34: /* vtype_fire_return_interval_range < vtype minFRI maxFRI > */
      *(param_valP) = 0.;
      *(param_valP + 1) = -9999.;
      *(param_valP + 2) = -9999.;
      break;
      
   // Next 4 are coefficients A-H and P from PNW-GTR-841 Table 3 p.25
   case 35: /* SSZ upper bound coefficients, equation E2 */
      /* A */ *(param_valP) = -647.f; 
      /* B */ *(param_valP + 1) = 0.204f;
      /* C */ *(param_valP + 2) = 45.7f;
      /* D */ *(param_valP + 3) = -29.63f;
      /* E */ *(param_valP + 4) = -10.7f;
      /* F */ *(param_valP + 5) = NC_FILL_FLOAT;
      /* G */ *(param_valP + 6) = 2.640f;
      /* H */ *(param_valP + 7) = 1.0f;
      /* P */ *(param_valP + 8) = 0.5f;
      break;
   case 36: /* PSFZ lower bound coefficients, equation E5 */
      /* A */ *(param_valP) = 1750.f; 
      /* B */ *(param_valP + 1) = -0.300f;
      /* C */ *(param_valP + 2) = -76.2f;
      /* D */ *(param_valP + 3) = 10.97f;
      /* E */ *(param_valP + 4) = 97.5f;
      /* F */ *(param_valP + 5) = -0.0132;
      /* G */ *(param_valP + 6) = -0.18f;
      /* H */ *(param_valP + 7) = -0.8f;
      /* P */ *(param_valP + 8) = 1.0f;
      break;
   case 37: /* SAFZ lower bound coefficients, equation E1 */
      /* A */ *(param_valP) = 1402.f; 
      /* B */ *(param_valP + 1) = 0.168f;
      /* C */ *(param_valP + 2) = -13.7f;
      /* D */ *(param_valP + 3) = 0.0f;
      /* E */ *(param_valP + 4) = -42.7f;
      /* F */ *(param_valP + 5) = NC_FILL_FLOAT;
      /* G */ *(param_valP + 6) = NC_FILL_FLOAT;
      /* H */ *(param_valP + 7) = NC_FILL_FLOAT;
      /* P */ *(param_valP + 8) = 0.5f;
      break;
   case 38: /* PKLZ lower bound coefficients, equation E1 */
      /* A */ *(param_valP) = 2847.f; 
      /* B */ *(param_valP + 1) = -0.540f;
      /* C */ *(param_valP + 2) = -6.1f;
      /* D */ *(param_valP + 3) = 10.97f;
      /* E */ *(param_valP + 4) = 18.3f;
      /* F */ *(param_valP + 5) = NC_FILL_FLOAT;
      /* G */ *(param_valP + 6) = NC_FILL_FLOAT;
      /* H */ *(param_valP + 7) = NC_FILL_FLOAT;
      /* P */ *(param_valP + 8) = 0.3f;
      break;
   case 39: /* ffmc_threshold <climate_zone tree_type FFMC_fire_threshold> */
   case 40: /* bui_threshold <climate_zone tree_type bui_fire_threshold> */
      *(param_valP) = 0.;
      *(param_valP + 1) = 0;
      *(param_valP + 2) = -9999.;
      break;
   case 43: /* ffmc_threshold_by_vtype <vtype ffmc_fire_threshold> */
   case 44: /* bui_threshold_by_vtype <vtype bui_fire_threshold> */
      *(param_valP) = 0.;
      *(param_valP + 1) = -9999.;
      break;
      
    default: for (jval=0; jval<nval; jval++) *(param_valP+jval) = info[ip].default_value; break;
  }
  
  if (buf!=NULL) 
  {
    strncpy(valbuf, buf, sizeof(valbuf)-1); /* Make a local, modifiable copy of the incoming string. */
    valbuf[sizeof(valbuf)-1] = 0;
    val_charP = strtok(valbuf, commaOrWhitespace);
    jval = 0;
    nread = 1;
    while (jval<nval && val_charP!=NULL && nread>0)
    {
      nread = sscanf(val_charP, "%f", param_valP+jval);
      if (nread>0)
      {
        if (*(param_valP+jval)<info[ip].lower_limit 
            || (!info[ip].inclusive_lower_limit && *(param_valP+jval)==info[ip].lower_limit)
            || *(param_valP+jval)>info[ip].upper_limit 
            || (!info[ip].inclusive_upper_limit && *(param_valP+jval)==info[ip].upper_limit))
        {
          printf("\n\nfrom command file: %s = %s", keyword, buf);
          printf("from internal parameter specs: default value = %f, lower limit = %f, upper limit = %f",
               info[ip].default_value, info[ip].lower_limit, info[ip].upper_limit);
          free(param_valP); // malloc in rpFloat()
          return (char *)info[ip].msg;
        }
        jval++;
        val_charP = strtok(NULL, commaOrWhitespace);
      }
    }
    if (jval>0) while (jval<nval) 
    {
      *(param_valP+jval) = *(param_valP+jval-1);
      jval++;
    }
  }

  int jval_of_first_actual_value = 0;
  switch (ip)
  {
    case 0: valP = &runParams.fire_suppression_erc_threshold; break;
    case 1: valP = &runParams.fire_suppression_fli_threshold; break;
    case 2: valP = &runParams.fire_suppression_ros_threshold; break;
    case 3: valP = &runParams.suppressed_fire_cell_fraction; break;
    case 6: valP = &modelParams.m_century_runoff_slope; break;
    case 7: valP = &modelParams.m_century_runoff_x_intercept; break;
    case 8: valP = &modelParams.efold_t_max; break;
    case 11: valP = &modelParams.fire_min; break;
    case 12: valP = &modelParams.m_lait_lower_limit; break;
    case 13: valP = &modelParams.maritime_threshold; break;
    case 14: valP = &modelParams.p_hi_mult; break; 
    case 15: valP = &modelParams.subalpine_threshold; break;
    case 16: valP = &modelParams.tmmin_threshold; break;
    case 17: valP = &modelParams.m_forest_thres_C; break;
    case 18: valP = &modelParams.m_woodl_thres_C; break;
    case 19: valP = &modelParams.MAPSS_subalpine_threshold; break;
    case 20: valP = &modelParams.desert_grass_C_max; break;
    case 21: valP = &modelParams.c3_threshold; break;
    case 22: valP = &modelParams.grassfrac_thres; break;
    case 23: valP = &modelParams.crown_fire_mortality_pct; break;
    case 24: valP = &modelParams.bz_thres; break;
    case 25: valP = &modelParams.double_CO2_tree_AT_mult; break;
    case 26: valP = &modelParams.double_CO2_tree_NPP_mult; break;
    case 27: valP = modelParams.max_LAI; break;
    case 28: valP = modelParams.max_tree_NPP; break;
    case 29: valP = &modelParams.double_CO2_grass_AT_mult; break;
    case 30: valP = &modelParams.double_CO2_grass_NPP_mult; break;
    case 31: valP = modelParams.max_grass_NPP; break;
    case 32: valP = modelParams.pprdwc_grass; break;
    case 33: valP = modelParams.pprdwc_tree; break;
    case 34: // vtype_fire_return_interval_range < vtype minFRI maxFRI >
      if (*(param_valP)>=0 && *(param_valP)<=(MAX_VTYPE+1))
            valP = &(modelParams.vveg2mfri[(int)(*(param_valP))][1]); 
      else valP = &(modelParams.vveg2mfri[0][1]);
      jval_of_first_actual_value = 1;
      break;
    case 35: valP = modelParams.SSZ_ub_coeffs; break;  
    case 36: valP = modelParams.PSFZ_lb_coeffs; break;  
    case 37: valP = modelParams.SAFZ_lb_coeffs; break;  
    case 38: valP = modelParams.PKLZ_lb_coeffs; break;  
    case 39: // ffmc_threshold <climate_zone tree_type ffmc_threshold>
      if (*(param_valP)>=0 && *(param_valP)<=5 && *(param_valP+1)>=0 && *(param_valP+1)<=8)
            valP = &(modelParams.ffmc_threshold[(int)(*(param_valP))][(int)(*(param_valP+1))]); 
      else valP = &(modelParams.ffmc_threshold[0][0]);
      jval_of_first_actual_value = 2;
      break;
    case 40: // bui_threshold <climate_zone tree_type bui_threshold>
      if (*(param_valP)>=0 && *(param_valP)<=5 && *(param_valP+1)>=0 && *(param_valP+1)<=8)
            valP = &(modelParams.bui_threshold[(int)(*(param_valP))][(int)(*(param_valP+1))]); 
      else valP = &(modelParams.bui_threshold[0][0]);
      jval_of_first_actual_value = 2;
      break;
    case 41: valP = &modelParams.moist_temperate_threshold; break;
    case 42: valP = &modelParams.dry_temperate_threshold; break;
    case 43: // ffmc_threshold_by_vtype <vtype ffmc_threshold>
      if (*(param_valP)>=1 && *(param_valP)<=MAX_VTYPE)
            valP = &(modelParams.ffmc_threshold_by_vtype[(int)(*(param_valP))]); 
      else valP = &(modelParams.ffmc_threshold_by_vtype[0]);
      jval_of_first_actual_value = 1;
      break;
    case 44: // bui_alt_threshold <climate_zone vtype bui_threshold>
      if (*(param_valP)>=1 && *(param_valP)<=MAX_VTYPE)
            valP = &(modelParams.bui_threshold_by_vtype[(int)(*(param_valP))]); 
      else valP = &(modelParams.bui_threshold_by_vtype[0]);
      jval_of_first_actual_value = 1;
      break;
    default: 
      free(param_valP); // malloc in rpFloat()
      return "Oops! rpFloat() internal error #2. <<<<<<<<<<<<<<<<";
  }
  for (jval=jval_of_first_actual_value; jval<nval; jval++) 
      *(valP+(jval-jval_of_first_actual_value)) = *(param_valP+jval);
  
  printf("value for %s", keyword);
  if (jval_of_first_actual_value>0) 
      for (jval=0; jval<jval_of_first_actual_value; jval++) printf(" %d", int(*(param_valP+jval)));
  printf(" is");
  for (jval=jval_of_first_actual_value; jval<nval; jval++) 
      printf(" %.6lf ", *(valP+(jval-jval_of_first_actual_value)));
  
  if (times_called>0 && jval_of_first_actual_value==0) 
      printf("\n\n*** Warning: %s specified more than once. <<<<<<<<<<<<<<<<\n", keyword);
  free(param_valP); // malloc in rpFloat()

   if (buf==NULL) // no argument string was passed in, so return the help string
         return (char *)info[ip].msg;
   else 
         return NULL;
} /* end of rpFloat() */


/*
typedef struct
{
  char * name;
  int max_num;
  double default_value;
  float lower_limit;
  float upper_limit;
  char * msg;  
} IntParamInfoStruct;
*/

const char * Simulation::rpIntFcn(const char * keyword, int times_called, char * buf)
{
  static const IntParamInfoStruct info[] = 
  {
/* 0 */ {"col_begin", 1, 0, 0, 9999, 
        "col_begin <first column in range of cells to simulate>"},
/* 1 */ {"col_end", 1, 0, 0, 9999,
        "col_end <final column in range of cells to simulate>"},
/* 2 */ {"fire_set_interval", 1, -1, 0, 9999,
        "fire_set_interval    <number of years between deliberately set fires> default = -1, meaning no deliberately set fires"},
/* 3 */ {"fire_set_jday", 1, 274, 1, 365,
        "fire_set_jday        <julian day (Jan 1 = 1) when fire is deliberately set> default = 274, meaning Oct 1"},
/* 4 */ {"fire_suppression_first_year", 1, 1951, 0, 9999,
        "fire_suppression_first_year <calendar year> default = 1951"},
/* 5 */ {"first_calendar_year_of_run", 1, 0, 0, 9999,
        "first_calendar_year_of_run <calendar year of first simulation year>"},
/* 6 */ {"multiyr_len", 1, 10, 1, 9999,
        "multiyr_len <length of multiyear interval for time-aggregated output variables> default = 10"},
/* 7 */ {"multiyr_start", 1, 1901, 0, 9999,
        "multiyr_start <first year of multiyear interval> default = 1901"},
/* 8 */ {"row_begin", 1, 0, 0, 9999,
        "row_begin <first row in range of cells to simulate>"}, 
/* 9 */ {"row_end", 1, 0, 0, 9999,
        "row_end <final row in range of cells to simulate>"},
/*10 */ {"years_offset_into_input_data", 1, 0, 0, 9999,
        "years_offset_into_input_data <year of first climate data to use relative to beginning of climate files>"},	
/*11 */ {"years_to_run", 1, 1, 1, 3000,
        "years_to_run <maximum years to run>"}	
  };	
  int num_of_params = sizeof(info)/sizeof(IntParamInfoStruct);
  int * param_valP;	
  int ip; /* index to info[] */
  int jval; /* index to particular value: param_val */
  int nval; /* number of values */
  int * valP; /* pointer to destination in MP or RunData */
  char valbuf[STRLEN+1];
  int nread;
  char * val_charP;
  
  ip = 0;
  while (ip<num_of_params && strcmp(keyword, info[ip].name)!=0) ip++;
  assert(ip<num_of_params);
  if (ip>=num_of_params) return "Oops! rpInt() internal error #1. <<<<<<<<<<<<<<<<";

  nval = info[ip].max_num;
  param_valP = (int *)malloc(nval*sizeof(int));
  if (param_valP==NULL) { printf("*** Error allocating for param_val array in rpInt().\n"); err_exit("MALLOC_FAILED"); }
  for (jval=0; jval<nval; jval++) *(param_valP+jval) = info[ip].default_value;
  
  if (buf!=NULL) 
  {
    strncpy(valbuf, buf, sizeof(valbuf)-1); /* Make a local, modifiable copy of the incoming string. */
    valbuf[sizeof(valbuf)-1] = 0;
    val_charP = strtok(valbuf, commaOrWhitespace);
    jval = 0;
    nread = 1;
    while (jval<nval && val_charP!=NULL && nread>0)
    {
      nread = sscanf(val_charP, "%d", param_valP+jval);
      if (nread>0)
      {
        if (*(param_valP+jval)<info[ip].lower_limit || *(param_valP+jval)>info[ip].upper_limit)
        {
          printf("\n\nfrom command file: %s = %s", keyword, buf);
          printf("from internal parameter specs: default value = %d, lower limit = %d, upper limit = %d",
               (int)info[ip].default_value, (int)info[ip].lower_limit, (int)info[ip].upper_limit);
          free(param_valP); // malloc in rpInt()
                return (char *)info[ip].msg;
        }
        jval++;
        val_charP = strtok(NULL, commaOrWhitespace);
      }
    }
    if (jval>0) while (jval<nval) 
    {
      *(param_valP+jval) = *(param_valP+jval-1);
      jval++;
    }
    /* special case default values would be set here, if there were any
    switch (ip)
    {
    }
    */
  }
      
  switch (ip)
  {
    case 0: valP = &runParams.col_begin; break;
    case 1: valP = &runParams.col_end; break;
    case 2: valP = &runParams.fire_set_interval; break;
    case 3: valP = &runParams.fire_set_jday; break;
    case 4: valP = &runParams.fire_suppression_first_year; break;
    case 5: valP = &runParams.first_calendar_year_of_run; break;
    case 6: valP = &runParams.multiyr_len; break;
    case 7: valP = &runParams.multiyr_start; break;
    case 8: valP = &runParams.row_begin; break;
    case 9: valP = &runParams.row_end; break;
    case 10: valP = &runParams.years_offset_into_input_data; break;
    case 11: valP = &runParams.years_to_run; break;
    default: { free(param_valP); return "Oops! rpInt() internal error #2. <<<<<<<<<<<<<<<<"; }   
  }
  for (jval=0; jval<nval; jval++) *(valP+jval) = *(param_valP+jval);
  
  printf("value for %s is", keyword);
  for (jval=0; jval<nval; jval++) printf(" %d", *(valP+jval));
  
  if (times_called>0) printf("\n\n*** Warning: %s specified more than once. <<<<<<<<<<<<<<<<\n", keyword); 
  free(param_valP); // malloc in rpInt()

   if (buf==NULL) // no argument string was passed in, so return the help string
         return (char *)info[ip].msg;
   else 
         return NULL;
} /* end of rpInt() */


const char * Simulation::rpCodeFlagFcn(const char * keyword, int n, char * buf)
/* code_flag <flag#> ON|OFF */   
{
  bool matchFlag;
  static const char * msg="Specify ON or OFF to set the flag TRUE or FALSE.";
  static const char * msg2="Index is out of range, or this code flag is overwritten by"
      "another parameter (e.g. alt_fuel_load_calculation, alt_tree_allometry_calculation)"
      ", and so cannot be set or cleared explicitly.";
  int i;
  int tgt_size;
  bool * tgt_adr;
  const char * tgt_name;
  char value[128];
  char valbuf[STRLEN+1];
  char * val_charP;
  
  tgt_size = NUM_CODE_FLAGS;
  tgt_adr = modelParams.code_flags;
  tgt_name = "code_flags";
  
  if (buf==NULL) 
  {
    for (i=0; i<NUM_CODE_FLAGS; i++) modelParams.code_flags[i] = false; 
    modelParams.code_flags[ALT_FUEL_LOAD_CODE_FLAG] = modelParams.altFuelLoadFlag;
    modelParams.code_flags[ALT_TREE_ALLOMETRY_CODE_FLAG] = modelParams.altTreeAllometryFlag;
    return("code_flag <flag#> < ON | OFF >");
  }
  
  /* Make a local, modifiable copy of the incoming string. */
  strncpy(valbuf, buf, sizeof(valbuf)-1); valbuf[sizeof(valbuf)-1] = 0;
  
  val_charP = strtok(valbuf, commaOrWhitespace);
  sscanf(val_charP, "%d", &i);
  if (i<0 || i>=tgt_size || i==ALT_FUEL_LOAD_CODE_FLAG || i==ALT_TREE_ALLOMETRY_CODE_FLAG) 
      return (char *)msg2;
  
  val_charP = strtok(NULL, commaOrWhitespace);
  strncpy(value, val_charP, sizeof(value)-1); value[sizeof(value)-1] = 0;    
  if (strcmp("ON", value)==0) matchFlag = true;
  else if (strcmp("OFF", value)==0) matchFlag = false;
  else return (char *)msg;

  *(tgt_adr+i) = matchFlag;
  printf("value for %s[%d] is %s", tgt_name, i, matchFlag ? "ON" : "OFF");
  return NULL;
} /* end of rpCodeFlagFcn() */


/*
typedef enum { sptArbitrary, sptOnOff, sptPath, sptFileName } StringParamType;
typedef struct
{
  char * name;
  int max_length;
  char * default_value;
  StringParamType param_type;
  char * msg;  
} StringParamInfoStruct;
*/


const char * Simulation::rpStringFcn(const char * keyword, int times_called, char * buf)
{
  static const StringParamInfoStruct info[] = 
  {
/* 0 */ {"climate_data_directory", PATHNAME_LENGTH, "./", sptPath, 
        "climate_data_directory  <path to input climate data root directory>"}, 
/* 1 */ {"CO2_file", PATHNAME_LENGTH, "NONE", sptFileName,
        "CO2_file <CO2 ramp path and file name>"},	 
/* 2 */ {"earth_data_directory", PATHNAME_LENGTH, "./", sptPath,
        "earth_data_directory <path to directory containing soil and elevation data>"},
/* 3 */ {"fire_mode_switch", 5, "ON", sptOnOff,
        "fire_mode_switch < ON | OFF >"},
/* 4 */	{"fire_suppression_switch", 5, "OFF", sptOnOff,		
        "fire_suppression_switch < ON | OFF >"},
/* 5 */ {"grid_name", PATHNAME_LENGTH, "NONE", sptArbitrary,	
        "grid_name < Global | NA8km | US50km | VEMAP | US12km | US10kmAlbers | US4km | US800m | BLM_PSW4kmAlbers | CA10kmAlbers | ...>"},
/* 6 */ {"mask_file", PATHNAME_LENGTH, "NONE", sptFileName,
        "mask_file <name of mask_file>"},
/* 7 obsolete */ {"model_parameter_set_obsolete", PATHNAME_LENGTH, "NONE", sptArbitrary, 			
        "model_parameter_set < VEMAP | NA8km | CA08 | YOSE | WIES_YOSE | WWETAC | VINCERA | US10km | US50km | GLOBAL >"},
/* 8 */ {"output_file_prefix", PATHNAME_LENGTH, "mc2", sptArbitrary,
        "output_file_prefix <prefix for output files>"},
/* 9 */ {"output_variable_file", PATHNAME_LENGTH, "NONE", sptFileName,
        "output_variable_file <name of file containing list of desired output variables>"},
/*10 */ {"run_mode", PATHNAME_LENGTH, "NONE", sptArbitrary,
        "run_mode < MAPSS_EQ | CENTURY_EQ | SPINUP | TRANSIENT >"},
/*11 */ {"soil_data_file", PATHNAME_LENGTH, "soils_scs.nc", sptFileName,
        "soil_data_file <file name and path of soil_data_file, relative to earth data directory>; defaults to 'soils_scs.nc'"},	
/*12 */ {"warmstart_file", PATHNAME_LENGTH, "NONE", sptFileName,
        "warmstart_file <path and file name of warmstart file, relative to working directory>"},
/*13 */ {"alt_fuel_load_calculation", 5, "ON", sptOnOff,
        "alt_fuel_load_calculation < ON | OFF >" },
/*14 */ {"MAPSS_parameter_set", PATHNAME_LENGTH, "NONE", sptArbitrary,
        "MAPSS_parameter_set < ORNLparams | US10kmSCSparams | US10kmFAOparams | WiesYoseParams"},	
/*15 */ {"alt_tree_allometry_calculation", 5, "ON", sptOnOff,					
        "alt_tree_allometry_calculation < ON | OFF >"},
/*16 */ {"southern_hemisphere_switch", 5, "OFF", sptOnOff,
        "southern_hemisphere_switch < ON | OFF > default = OFF"},
/*17 */ {"unlimited_N_switch", 5, "ON", sptOnOff,
        "unlimited_N_switch < ON | OFF >"},
/*18 */ {"space_before_time_switch", 5, "OFF", sptOnOff,
        "space_before_time_switch < ON | OFF >"},
/*19 */ {"century_path", PATHNAME_LENGTH, "Input/ModelParameters_MC2/", sptFileName, 
         "century_path <path and directory name of the CENTURY data directory, relative to working directory>"},
/*20 */ {"output_variable", PATHNAME_LENGTH, "no default value", sptVarName, 
         "output_variable <output variable name> all caps, letters, digits, and underscores only; variable must exist"
         " in discretionary output variable dictionary"},
/*21 */ {"monthly_output_switch", 5, "OFF", sptOnOff, 
         "monthly_output_switch < ON | OFF >"},
/*22 */ {"yearly_output_switch", 5, "ON", sptOnOff, 
         "yearly_output_switch < ON | OFF >"},
/*23 */ {"multiyr_output_switch", 5, "ON", sptOnOff, 
         "multiyr_output_switch < ON | OFF >"},
/*24 */ {"time_invariant_output_switch", 5, "ON", sptOnOff, 
         "time_invariant_output_switch < ON | OFF >"},
/*25 */ {"warmstart_output_switch", 5, "OFF", sptOnOff, 
         "warmstart_output_switch < ON | OFF >"},
/*26 */ {"dummy_wind_switch", 5, "OFF", sptOnOff, 
         "dummy_wind_switch < ON | OFF >; default is OFF; ON means use dummy wind speed values instead"
         " of reading a wind data file"},
/*27 */ {"soil_bulk_density_file", PATHNAME_LENGTH, "bd.nc", sptFileName,
        "soil_bulk_density_file <file name and path of soil_bulk_density_file, relative to earth data directory>; defaults to 'bd.nc'"}	
  };	
  int num_of_params = sizeof(info)/sizeof(StringParamInfoStruct);
  char ** param_strPP = NULL; // pointer to destination in runParams or modelParams
  char * param_valP;	
  int ip; 			/* index to info[] */
  int max_length; 		/* max string length */
  int len;
  bool * bool_valP = NULL; // pointer to destination in runParams or modelParams
  char valbuf[STRLEN+1];
  char * valbuf_charP;
  char * charP;
  const char whitespace[]=" \t\n";
  bool errflag=FALSE;
  
  ip = 0;
  while (ip<num_of_params && strcmp(keyword, info[ip].name)!=0) ip++;
  if (ip>=num_of_params) return "Oops! rpString() internal error #1. <<<<<<<<<<<<<<<<";
  assert(strlen(info[ip].default_value)>0);

  max_length = info[ip].max_length;
  param_valP = (char *)malloc(max_length+1);
  if (param_valP==NULL) { printf("*** Error allocating for param_val array in rpString().\n"); err_exit("MALLOC_FAILED"); }
  switch (ip) 
  { // default values
    case 19: // century_path
      switch (m_baseCalibration)
      {
        case mc2California:
          strncpy(param_valP, "Input/ModelParameters_MC2_California", max_length); *(param_valP+max_length) = 0;
          break;
        default: strncpy(param_valP, info[ip].default_value, max_length); *(param_valP+max_length) = 0;
          break;
      }
      break;
    default: strncpy(param_valP, info[ip].default_value, max_length); *(param_valP+max_length) = 0; 
      break;
  }
  
  if (buf!=NULL) 
  {
    strncpy(valbuf, buf, sizeof(valbuf)-1); /* Make a local, modifiable copy of the incoming string. */
    valbuf[sizeof(valbuf)-1] = 0;
    valbuf_charP = strtok(valbuf, whitespace);
    errflag = (valbuf_charP==NULL);
    if (!errflag) switch (info[ip].param_type)
    {
      case sptOnOff: errflag = (strcmp(valbuf_charP, "ON")!=0) && (strcmp(valbuf_charP, "OFF")!=0); break;
      case sptPath: break;
      case sptArbitrary: break;
      case sptFileName: break;
      case sptVarName: errflag = !addVarToList(valbuf_charP); break;
      default: { free(param_valP); return "Oops! rpString() internal error #2. <<<<<<<<<<<<<<<<"; }
    }
    if (errflag) 
    { 
          printf("\n\nfrom command file: %s = %s", keyword, buf);
          printf("from internal parameter specs: default value = %s, max string length = %d",
               info[ip].default_value, info[ip].max_length);
      free(param_valP); // malloc in rpString()
      return (char *)info[ip].msg; 
    }
    strncpy(param_valP, valbuf_charP, max_length); *(param_valP+max_length) = 0; 
  }   
      
  switch (ip)
  {
    case 0: param_strPP = &runParams.climate_data_directory; break;
    case 1: param_strPP = &runParams.CO2_file; break;
    case 2: param_strPP = &runParams.earth_data_directory; break;
    case 3: bool_valP = &runParams.fireModelFlag; break;
    case 4: bool_valP = &runParams.fireSuppressionFlag; break;
    case 5: param_strPP = &runParams.grid_name; break;
    case 6: param_strPP = &runParams.mask_file; break;
    // case 7: param_strPP = &runParams.model_parameter_set; break;
    case 8: param_strPP = &runParams.output_file_prefix; break;
    case 9: param_strPP = &runParams.output_variable_file; break;
    case 10: param_strPP = &runParams.run_mode; break;
    case 11: param_strPP = &runParams.soil_data_file; break;
    case 12: param_strPP = &runParams.warmstart_file; break;
    case 13: // bool_valP = &modelParams.code_flags[ALT_FUEL_LOAD_CODE_FLAG]; break; 
             bool_valP = &modelParams.altFuelLoadFlag; break; 
    case 14: param_strPP = &modelParams.MAPSS_parameter_set; break;
    case 15: // bool_valP = &modelParams.code_flags[ALT_TREE_ALLOMETRY_CODE_FLAG]; break;
             bool_valP = &modelParams.altTreeAllometryFlag; break; 
    case 16: bool_valP = &modelParams.southernHemisphereFlag; break;
    case 17: bool_valP = &modelParams.unlimitedNflag; break;
    case 18: bool_valP = &runParams.spaceBeforeTimeFlag; break;
    case 19: param_strPP = &modelParams.century_path; break;
    case 20: // "output_variable"  The work gets done in the addVarToList() call above.
      times_called = 0; // Using output_variable more than once is OK.
      break;
    case 21: bool_valP = &runParams.monthlyOutputFlag; runParams.monthlyOutputSwitchCount++; break;
    case 22: bool_valP = &runParams.yearlyOutputFlag; runParams.yearlyOutputSwitchCount++; break;
    case 23: bool_valP = &runParams.multiyrOutputFlag; runParams.multiyrOutputSwitchCount++; break;
    case 24: bool_valP = &runParams.timeInvariantOutputFlag; runParams.timeInvariantOutputSwitchCount++; break;
    case 25: bool_valP = &runParams.warmstartOutputFlag; runParams.warmstartOutputSwitchCount++; break;    
    case 26: bool_valP = &runParams.dummyWindFlag; break;    
    case 27: param_strPP = &runParams.soil_bulk_density_file; break;
    
    default: { free(param_valP); return "Oops! rpString() internal error #3. <<<<<<<<<<<<<<<<"; }
  }
  switch (info[ip].param_type)
      {
      case sptOnOff: 
         *bool_valP = (strcmp(param_valP, "ON")==0);
         break;
      case sptArbitrary:
      case sptPath:
      case sptFileName:
         /* First copy the string to the destination. */ 
         if (*param_strPP!=NULL) free(*param_strPP);
         *param_strPP = (char *)malloc(strlen(param_valP) + 2); assert(*param_strPP!=NULL); // never free'd
         strncpy(*param_strPP, param_valP, strlen(param_valP) + 1);
         *(*param_strPP + (strlen(param_valP) + 1)) = 0;

         /* Now, if it's a path or a file name, replace NONE with a null string. */
         if ((info[ip].param_type==sptPath || info[ip].param_type==sptFileName) && strcmp("NONE", param_valP)==0)
           *(*param_strPP) = 0;

         /* Finally, if it's a path, make sure that the path ends with a '/'. */
         if (info[ip].param_type==sptPath)
         { 
         len = strlen(*param_strPP);
         if (len>0 && *(*param_strPP+len-1)!='/')
            {
            charP = *param_strPP + (len<max_length ? len : len - 1);
            *charP = '/';
            *(charP+1) = 0;
            }
         }
         break;
      case sptVarName:
         break;
      default:
         assert(false);
         break;
      }

  assert(param_valP!=NULL);
  printf("value for %s is %s", keyword, param_valP);
  
  if (times_called>0) printf("\n\n*** Warning: %s specified more than once. <<<<<<<<<<<<<<<<\n", keyword); 
  free(param_valP);  // malloc in rpString()
  
   if (buf==NULL) // no argument string was passed in, so return the help string
         return (char *)info[ip].msg;
   else 
         return NULL;
} /* end of rpString() */


char * tokncpy(char * destinationP, char * sourceP, size_t maxsize, char const * delimitersP)
/* Sort of a cross between strtok(sourceP,delimitersP) and strncpy(destinationP,sourceP,maxsize).
   Copies up to maxsize characters of the first token in source to destination.
   Also puts a zero at destinationP+maxsize.
   Returns the address of the next token in source. */
{
  size_t sourceSize;
  char * local_sourceP;
  char * tokP;
  char * return_val;
  
  *destinationP = 0;
  if (sourceP==NULL) return(NULL);
  
  sourceSize = strlen(sourceP);
  if (sourceSize<=0) return(NULL);
  
  local_sourceP = (char *)malloc(sourceSize + 1);
  if (local_sourceP==NULL) { printf("*** Error allocating for local_source in tokncpy().\n"); err_exit("MALLOC_FAILED"); }
  
  strcpy(local_sourceP, sourceP);
  tokP = strtok(local_sourceP, delimitersP);
  
  if (tokP==NULL) return_val = NULL;
  else
  {
    strncpy(destinationP, tokP, maxsize); *(destinationP+maxsize) = 0;
    tokP = strtok(NULL, delimitersP);
    return_val = tokP==NULL ? NULL : sourceP + (tokP - local_sourceP);
  }
  
  free(local_sourceP); // malloc in tokncpy()
  return(return_val);  
  
} /* end of tokncpy() */  


void Simulation::access_file(FILE ** file, char * filename, const char * status)
   {
	*file = fopen(filename, status);
   if (*file!=NULL) return; 
   (void) printf("access_file: %s\n", filename);
   err_exit("Cannot access file");	
   } // end of access_file()


void Simulation::writeOutputVarList(const char * file_name, OutVarDescription List[], int list_len)
/* This procedure creates a text file  with information
  about the output variables described in the list. */
{
  FILE * fileP;
  int i, k;
  
  fileP = fopen(file_name, "w");
  
  fprintf(fileP, "## %s\n\n", file_name);
  fprintf(fileP, "##%-30s\t%-20s\t%-13s\t%-6s\t%-s\n",
      "Variable Name", "Units", "Interval", "Type", "Description");

  k = 0;
  for (i=0; i<list_len; i++) if (strcmp(List[i].name, "spareOutVar")!=0)
      {
      fprintf(fileP, "#output_variable %-30s\t%-20s\t%-6s\t%-s\n",
            List[i].name, List[i].units, List[i].type==NC_SHORT ? "int" : "float", List[i].descrip);
      k++;
      }
  fprintf(fileP, "\nnumber of variables = %d\n", k);   
   
  fclose(fileP);
  
} // end of writeOutputVarList()


void instructions();

void Simulation::MC2instructions()
   {
   const char * charP;
  
   printf("\n*** MC2 help text\n");
   
   printf("\nRun parameters which may be set in the command file, and their default values, are:\n\n");
   for (unsigned int i=0; i<NUM_OF_RUN_PARAMS; i++) if (RunParamsList[i].fcnType!=rpUnknown)
      { 
      printf("%s: default ", RunParamsList[i].keyword); 
      switch (RunParamsList[i].fcnType)
         {
         case rpString: charP = rpStringFcn(RunParamsList[i].keyword, 0, NULL); break;
         case rpInt: charP = rpIntFcn(RunParamsList[i].keyword, 0, NULL); break;
         case rpFloat: charP = rpFloatFcn(RunParamsList[i].keyword, 0, NULL); break;
         case rpCodeFlag: charP = rpCodeFlagFcn(RunParamsList[i].keyword, 0, NULL); break;
         default: assert(0); break;
         } // end of switch on fcnType
      printf("\n%s\n\n", charP);
      } // end of loop thru all the run parameters
   printf("\n");
   
   printf("\nModel parameters which may be set in the command file, and their default values, are:\n\n");
   for (unsigned int i=0; i<NUM_OF_MODEL_PARAMS; i++) if (ModelParamsList[i].fcnType!=rpUnknown)
      {
      printf("%s: default ", ModelParamsList[i].keyword); 
      switch (ModelParamsList[i].fcnType)
         {
         case rpString: charP = rpStringFcn(ModelParamsList[i].keyword, 0, NULL); break;
         case rpInt: charP = rpIntFcn(ModelParamsList[i].keyword, 0, NULL); break;
         case rpFloat: charP = rpFloatFcn(ModelParamsList[i].keyword, 0, NULL); break;
         case rpCodeFlag: charP = rpCodeFlagFcn(ModelParamsList[i].keyword, 0, NULL); break;
         default: assert(0); break;
         } // end of switch on fcnType
      printf("\n%s\n\n", charP);
      } // end of loop thru all the model parameters
   printf("\n\n");

   // Initialize the discretionary output variable dictionary and output lists.
   initializeVarDict();
   for (int index = 0; index<NUM_OPTIONAL_OUTVARS; index++)
         m_moOutList[index].show = m_yrOutList[index].show = m_multiyrOutList[index].show  = false;             
   m_num_monthly_vars = 0;
   m_num_yearly_vars = 0;
   m_num_multiyr_vars = 0;
   m_num_single_vars = 0;

   // Load up a few default output variables.
   addVarToYrList(VTYPE);
   addVarToMultiyrList(VTYPE);
   addVarToYrList(PHYSIOGNOMIC_CLASS);
   addVarToMultiyrList(PHYSIOGNOMIC_CLASS);
   addVarToYrList(BIOME);
   addVarToMultiyrList(BIOME);
   addVarToYrList(C_VEG);
   addVarToMultiyrList(C_VEG);
   addVarToYrList(NPP);
   addVarToMultiyrList(NPP);
   addVarToYrList(PART_BURN);
   addVarToMultiyrList(PART_BURN);
   
   printf("\n\nDefault Output Variables\n\n");
   printf("%-30s\t%-20s\t%-13s\t%-6s\t%-s\n",
         "Variable Name", "Units", "Interval", "Type", "Description");
   for (int i=0; i<m_num_yearly_vars; i++) printf("%-30s\t%-20s\t%-13s\t%-6s\t%-s\n",
            m_yrOutList[i].name, m_yrOutList[i].units, "yearly", m_yrOutList[i].type==NC_SHORT ? "int" : "float", m_yrOutList[i].descrip);
   for (int i=0; i<m_num_multiyr_vars; i++) printf("%-30s\t%-20s\t%-13s\t%-6s\t%-s\n",
            m_multiyrOutList[i].name, m_multiyrOutList[i].units, "multiyr", m_multiyrOutList[i].type==NC_SHORT ? "int" : "float", m_multiyrOutList[i].descrip);
   
   printf("\n");
   instructions();
   
   } // end of Simulation::MC2instructions()


 
  

