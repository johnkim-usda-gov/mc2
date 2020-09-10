/*
 *  MC2.cpp
 */

#include <stdlib.h>	// Required for malloc(). JK.
#include <string.h>	// JK.
#include <stdio.h>	// JK.
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <vector>

#include "netcdf.h"
#include "gridspecs.h"
#include "category_bgc.h"
#include "assert.h"

#include "MAPSSvegClasses.h"
#include "ScienceFcns.h"

#include "ProcessModel.h"

#include "MC2.h"
// #include "ScienceFcns.h"
#include "MAPSSbiogeographyModel.h"
#include "MCfire.h"
#include "MCbiogeog.h"
#include "CENTURY.h"

#define COMMAND_NAME "mc2"

using namespace std;

ScienceFcns sciFn;
bool interpretSwitch(char * switchStr);


int main(int argc, char * argv[])
  {
  Simulation theRun(argc, argv);

  if (theRun.m_baseCalibration==unknownBaseCalibration) err_exit("did not recognize calibration name");

  theRun.run();
  
  } // end of main()


Simulation::Simulation(int argc, char * argv[])
  {  
  if (argc!=3) MC2instructions();
  
  // m_cmdLine is free'd in the destructor
  m_cmdLine = (char *)malloc(strlen(argv[0]) + strlen(argv[1]) + strlen(argv[2]) + 3); assert(m_cmdLine!=NULL);
  sprintf(m_cmdLine, "%s %s %s", argv[0], argv[1], argv[2]);
  printf ("mc2 invocation: %s\n", m_cmdLine);
  
  char * baseCalibrationName = argv[1];
  char * cmdFileName = argv[2];
  mask_ncid.fileid = -1;
  soilData_ncid.fileid = -1;
  bd_ncid.fileid = -1;
  ppt_ncid.fileid = -1;
  tmin_ncid.fileid = -1;
  tmax_ncid.fileid = -1;
  tmp_ncid.fileid = -1;
  vpr_ncid.fileid = -1;
  tdmean_ncid.fileid = -1;
  elev_ncid.fileid = -1;
  wnd_ncid.fileid = -1;
  
  fog_ncid.fileid = -1;
  topomoist_ncid.fileid = -1;
  deltaTsl_ncid.fileid = -1;
  aspect_ncid.fileid = -1;
  sw_ncid.fileid = -1;
  cad_ncid.fileid = -1;
  
  MAPSSdata_ncid.fileid = -1;
  EQdata_ncid.fileid = -1;
  WSdata_ncid.fileid = -1;
  m_moOutFile_ncid.fileid = -1; m_num_monthly_vars = 0;
  m_yrOutFile_ncid.fileid = -1; m_num_yearly_vars = 0;
  m_multiyrOutFile_ncid.fileid = -1; m_num_multiyr_vars = 0;
  m_singleOutFile_ncid.fileid = -1; m_num_single_vars = 0;
  
  // enum MAPSSparameterSetName {defaultParams, ORNLparams, US10kmFAOparams, US10kmSCSparams, WiesYoseParams };
  MpSN[0] = (char *)"defaultParams";
  MpSN[1] = (char *)"ORNLparams";
  MpSN[2] = (char *)"US10kmFAOparams";
  MpSN[3] = (char *)"US10kmSCSparams";
  MpSN[4] = (char *)"WiesYoseParams"; 
    
  m_cmdFileName = cmdFileName;
  m_baseCalibration = interpretCalibrationName(baseCalibrationName);
  if (m_baseCalibration==unknownBaseCalibration) return;  
  getModelParameters(m_baseCalibration);
  
  getRunParameters(cmdFileName);
  
  { // Pass parameter values and CO2 ramp into Century
  extern CenParamStruct treeParams[];
  extern CenParamStruct cropParams[];
  extern CenParamStruct pprdwcParams[];
  int rtnval;

  rtnval = setCenturyParameter(treeParams, "CO2ITR", 0, modelParams.double_CO2_tree_AT_mult); assert(rtnval==1);
  rtnval = setCenturyParameter(treeParams, "CO2IPR", 0, modelParams.double_CO2_tree_NPP_mult); assert(rtnval==1);
  for (int jval = 0; jval<5; jval++)
      {
      rtnval = setCenturyParameter(treeParams, "MAXLAI", jval, modelParams.max_LAI[jval]); assert(rtnval==1);
      rtnval = setCenturyParameter(treeParams, "PRDX(4)", jval, modelParams.max_tree_NPP[jval]); assert(rtnval==1);
      }
  rtnval = setCenturyParameter(cropParams, "CO2ITR", 0, modelParams.double_CO2_grass_AT_mult); assert(rtnval==1);
  rtnval = setCenturyParameter(cropParams, "CO2IPR", 0, modelParams.double_CO2_grass_NPP_mult); assert(rtnval==1);
  for (int jval = 0; jval<3; jval++)
      {
      rtnval = setCenturyParameter(cropParams, "PRDX(1)", jval, modelParams.max_grass_NPP[jval]); assert(rtnval==1);
      rtnval = setCenturyParameter(pprdwcParams, "pprdwc_grass", jval, modelParams.pprdwc_grass[jval]); assert(rtnval==1);    
      rtnval = setCenturyParameter(pprdwcParams, "pprdwc_tree", jval, modelParams.pprdwc_tree[jval]); assert(rtnval==1);    
      }
      
  if (runParams.runMode!=MAPSS_EQ)
      {
      if (strlen(runParams.CO2_file)>0 && strcmp(runParams.CO2_file, "NONE")!=0)
         {
         rtnval = read_co2ramp_C(runParams.CO2_file, runParams.first_calendar_year_of_run, runParams.years_to_run, runParams.years_offset_into_input_data);
         if (rtnval!=1) 
            {
            printf("runParams.CO2_file = %s\n", runParams.CO2_file);
            err_exit("Error when reading CO2 ramp file.");
            }
         }
      else if (runParams.runMode==CENTURY_EQ || runParams.runMode==SPINUP)
            set_constant_CO2(294.842); // preindustrial atmospheric CO2 conc, ppm
      else err_exit("Missing spec for CO2 ramp");
      }
  } // end of block to pass parameter values and CO2 ramp into Century
  
  int num_rows = runParams.row_end - runParams.row_begin + 1;
  int num_cols = runParams.col_end - runParams.col_begin + 1;
  
  // m_pActiveCellArray is free'd in the destructor
  m_pActiveCellArray = (bool *)malloc(num_cols*num_rows*sizeof(bool)); assert(m_pActiveCellArray!=NULL);  
  
   if (runParams.latlonFlag) 
      { float N_lat, S_lat;
      N_lat = north_south_centroid(runParams.row_begin);
      S_lat = north_south_centroid(runParams.row_end);
   
      bool localNorthernHemisphereFlag = N_lat>0.;
      if (localNorthernHemisphereFlag && modelParams.southernHemisphereFlag) printf("*** Warning: Simulation begins in the northern hemisphere, but modelParams.southernHemisphereFlag is set.\n");

      bool localSouthernHemisphereFlag = S_lat<0.;
      if (localSouthernHemisphereFlag && !modelParams.southernHemisphereFlag) printf("*** Warning: Simulation ends in the southern hemisphere, but modelParams.southernHemisphereFlag is not set.\n");
      }
  
  m_num_vpr_tmp_issues = 0;
  m_num_target_cells = m_num_active_cells = m_num_failed_cells = -9999;
  
  } // end of Simulation(baseCalibrationName, cmdFileName, cmdline) constructor


Simulation::~Simulation() 
   { 
   if (m_cmdLine!=NULL) free(m_cmdLine);
   if (m_pActiveCellArray!=NULL) free(m_pActiveCellArray);
   };


void Simulation::run()
   {
   int num_rows = runParams.row_end - runParams.row_begin + 1;
   int num_cols = runParams.col_end - runParams.col_begin + 1;
   bool rtn_flag;
   
   rtn_flag = openInputFiles(); 
   if (!rtn_flag)
   {
      printf("*** run(): One or more necessary input files is missing, inaccessible, or incomplete.\n");
      return;
   }

   m_num_target_cells = initializeActiveCellArray();
   
   // Now that we know the length of the climate data and the specified number of years to run,
   // determine the actual number of years to run and the number of multiyr intervals.
   runParams.actual_years_to_run = runParams.years_to_run;
   switch (runParams.runMode)
      { int first_calendar_year_of_interval;
      case MAPSS_EQ: runParams.actual_years_to_run = 1; break;
      case CENTURY_EQ: break; // runParams.actual_years_to_run = runParams.years_to_run
      case SPINUP:
         if (runParams.yearlyOutputSwitchCount==0) runParams.yearlyOutputFlag = false;
      case TRANSIENT:
         assert(runParams.multiyr_len>0);
         assert(runParams.multiyr_start>=0);

         if (runParams.runMode==SPINUP) 
            {
            runParams.actual_years_to_run = 
               runParams.years_of_climate_data*(runParams.years_to_run/runParams.years_of_climate_data)
               + (runParams.years_to_run%runParams.years_of_climate_data>0 ? 
               runParams.years_of_climate_data : 0);
            assert(runParams.first_calendar_year_of_run==0);
            runParams.multiyr_len = runParams.years_of_climate_data;
            }
         else 
            { // TRANSIENT run mode
            assert(runParams.runMode==TRANSIENT);
            if (runParams.first_calendar_year_of_run==0) 
                  err_exit("Please specify first_calendar_year_of_run.\n");
                  
            int yrs_avail = runParams.years_of_climate_data - runParams.years_offset_into_input_data;
            // Since climate data starts in January, but southern hemisphere years start in July, 
            // when we're in the southern hemisphere we need an extra year of climate data
            if (modelParams.southernHemisphereFlag) yrs_avail--;
            runParams.actual_years_to_run = runParams.years_to_run>yrs_avail ? yrs_avail : runParams.years_to_run;
            }

         runParams.nominal_multiyr_start = runParams.runMode==SPINUP ? 0 : runParams.multiyr_start;
         if (runParams.nominal_multiyr_start<runParams.first_calendar_year_of_run)
            {
            while (runParams.nominal_multiyr_start<=runParams.first_calendar_year_of_run)
                  runParams.nominal_multiyr_start += runParams.multiyr_len;
            runParams.nominal_multiyr_start -= runParams.multiyr_len;
            }
         runParams.last_calendar_year_of_run = runParams.first_calendar_year_of_run + runParams.actual_years_to_run - 1;
               
         // During spinup, allow 1 extra interval for saving initial conditions.
         runParams.num_multiyr_intervals = runParams.runMode==SPINUP ? 1 : 0;
         
         first_calendar_year_of_interval = runParams.nominal_multiyr_start;
         while (first_calendar_year_of_interval<=runParams.last_calendar_year_of_run)
            {
            runParams.num_multiyr_intervals++;
            first_calendar_year_of_interval += runParams.multiyr_len;
            }
         break;
      default: assert(false); break;
      }
      
  initializeOutputFiles();
  
  printf("spaceBeforeTimeFlag, years_to_run, num_rows, num_cols = %s, %d, %d, %d\n", 
      runParams.spaceBeforeTimeFlag ? "true" : "false", runParams.years_to_run, num_rows, num_cols);
  int m_num_failed_cells = 0;
  for (int year_index = 0; year_index<runParams.actual_years_to_run; year_index++) 
    {
    m_num_active_cells = 0;
    for (int row_index = 0; row_index<num_rows; row_index++)
    for (int col_index = 0; col_index<num_cols; col_index++)
      { bool activeCell_flag;
      activeCell_flag = *(m_pActiveCellArray + row_index*num_cols + col_index); 
      if (!activeCell_flag) continue;
      
      m_pProcessModelList = new std::vector<ProcessModel *>;
    
      switch (runParams.runMode)
         {
         case MAPSS_EQ:
            runParams.spaceBeforeTimeFlag = false;
            runParams.monthlyOutputFlag = runParams.yearlyOutputFlag = runParams.multiyrOutputFlag = runParams.timeInvariantOutputFlag = false;
            m_pProcessModelList->push_back(new MAPSS_BiogeogModel(this));
            break;
         case CENTURY_EQ:
            runParams.spaceBeforeTimeFlag = false;
            runParams.monthlyOutputFlag = runParams.yearlyOutputFlag = runParams.multiyrOutputFlag = runParams.timeInvariantOutputFlag = false;
            m_pProcessModelList->push_back(new CENTURY_BiogeochemModel(this));
            break;
         case SPINUP:
            runParams.spaceBeforeTimeFlag = false;
            m_pProcessModelList->push_back(new CENTURY_BiogeochemModel(this));
            break;
         case TRANSIENT:
            runParams.spaceBeforeTimeFlag = false; // true;
            m_pProcessModelList->push_back(new CENTURY_BiogeochemModel(this));
            // if (runParams.fireModelFlag) m_pProcessModelList->push_back(new FireModel(this));
            // m_pProcessModelList->push_back(new BiogeogModel(this));
            break;
         case unknownRunMode:
         default:
            assert(0);
            break;
         }
          
      std::vector<ProcessModel*>::iterator modelIterator;

      for (modelIterator = m_pProcessModelList->begin(); activeCell_flag && modelIterator<m_pProcessModelList->end(); ++modelIterator)
         {
         activeCell_flag = activeCell_flag && (*modelIterator)->runModelAndCollectData(year_index, row_index, col_index);
         delete *modelIterator;
         }

      *(m_pActiveCellArray + row_index*num_cols + col_index) = activeCell_flag;
      if (activeCell_flag) 
         {
         m_num_active_cells++;
         printf("year_index, grid_row, grid_col, num_active_cells = %d, %d, %d, %d\n", 
               year_index, row_index + runParams.row_begin, col_index + runParams.col_begin, m_num_active_cells);
         }
      else 
         {
         m_num_failed_cells++;
         if (mask_ncid.fileid!=-1)
               printf("Simulation::run(): Active cell failed, grid_row, grid_col, num_failed_cells = %d, %d, %d\n",
               row_index + runParams.row_begin, col_index + runParams.col_begin, m_num_failed_cells);
         }
         
      delete m_pProcessModelList;
      
      } // end of loops on row and col
      
    printf("\nAt end of this iteration across space, year_index, num_active_cells = %d, %d\n", year_index, m_num_active_cells);
    if (!runParams.spaceBeforeTimeFlag) break;
    } // end of loop on year_index
  
  if (mask_ncid.fileid==-1) printf("No mask file was found.\n");
  else printf("mask file %s was read to determine which cells are active.\n", mask_ncid.file_name_and_path);
  printf("num_target_cells = %d\n", m_num_target_cells);
  printf("num_active_cells = %d\n", m_num_active_cells);
  printf("num_failed_cells = %d\n", m_num_failed_cells);
  if (m_num_target_cells!=(m_num_active_cells + m_num_failed_cells))
      printf("\n*** warning: num_active_cells + num_failed_cells is not equal to num_target_cells\n");
  saveOutputData();
  
  } // end of Simulation::run()


void Simulation::describeFSvars(OutVarDescription * FSoutputList)
   {
   // int dimids[3];
  
   // dimids[0] = ncidP->time_dimid;    
   // dimids[1] = ncidP->row_dimid;
   // dimids[2] = ncidP->col_dimid;
  
   FSoutputList[FSsum_ann_ppt] = describeVar("sum_ann_ppt", "", "", NC_FLOAT);
   FSoutputList[FS_Pprev] = describeVar("Pprev", "", "", NC_FLOAT);
   FSoutputList[FS_Kprev] = describeVar("Kprev", "", "", NC_FLOAT);
   FSoutputList[FS_Ksum] = describeVar("Ksum", "", "", NC_FLOAT);
   FSoutputList[FSprev_l1hr] = describeVar("FSprev_l1hr", "", "", NC_FLOAT);
   FSoutputList[FSprev_dstand] = describeVar("FSprev_dstand", "", "", NC_FLOAT);
   FSoutputList[FSprev_litter] = describeVar("FSprev_litter", "", "", NC_FLOAT);
   FSoutputList[FSprev_1hr] = describeVar("FSprev_1hr", 
         "1hr dead fuel in previous month", "g C m-2", NC_FLOAT);
   FSoutputList[FSprev_10hr] = describeVar("FSprev_10hr", 
         "10hr dead fuel in previous month", "g C m-2", NC_FLOAT);
   FSoutputList[FSprev_100hr] = describeVar("FSprev_100hr", 
         "100hr dead fuel in previous month", "g C m-2", NC_FLOAT);
   FSoutputList[FSprev_1000hr] = describeVar("FSprev_1000hr", 
         "1000hr dead fuel in previous month", "g C m-2", NC_FLOAT);
   FSoutputList[FSprev_day_mc_100hr] = describeVar("FSprev_day_mc_100hr", 
         "100hr dead fuel moisture content in previous day", "moisture content %", NC_FLOAT);
   FSoutputList[FSprev_day_mc_1000hr] = describeVar("FSprev_day_mc_1000hr", 
         "1000hr dead fuel moisture content in previous day", "moisture content %", NC_FLOAT);
   FSoutputList[FSprev_mxd] = describeVar("FSprev_mxd", "", "", NC_FLOAT);
   FSoutputList[FSprev_mc_grass] = describeVar("FSprev_mc_grass", "", "", NC_FLOAT);
   FSoutputList[FSprev_mc_tree] = describeVar("FSprev_mc_tree", "", "", NC_FLOAT);
   FSoutputList[FSprev_depth] = describeVar("FSprev_depth", "", "", NC_FLOAT);
   FSoutputList[FSclai_in_midseason] = describeVar("clai_in_midseason", "", "", NC_FLOAT);
   FSoutputList[FSprev_snow] = describeVar("prev_snow", "", "", NC_FLOAT);
   FSoutputList[FS_F1] = describeVar("F1", "", "", NC_FLOAT);
   FSoutputList[FS_F2] = describeVar("F2", "", "", NC_FLOAT);
   FSoutputList[FSrand_seed] = describeVar("FSrand_seed", "", "", NC_FLOAT);
   FSoutputList[FSyrs_in_sum_ann_ppt] = describeVar("yrs_in_sum_ann_ppt", "", "", NC_FLOAT);
   FSoutputList[FSyrs_since_fire] = describeVar("yrs_since_fire", "years since last unsuppressed fire", "yrs", NC_FLOAT);
   FSoutputList[FSffmc_prev] = describeVar("FSffmc_prev", "", "", NC_FLOAT);
   FSoutputList[FSdmc_prev] = describeVar("FSdmc_prev", "", "", NC_FLOAT);
   FSoutputList[FSdc_prev] = describeVar("FSdc_prev", "", "", NC_FLOAT);
   
   } // end of Simulation::describeFSvars()


bool Simulation::openInputDataNCfile(ncid_type * ncidP, char * directory, const char * file_name)
   {
   char * var_name;
   char suffix[4];
   bool rtnflag;
   
   assert(strlen(file_name)>3);
   strcpy(suffix, file_name + (strlen(file_name) - 3));
   assert(strcmp(suffix, ".nc")==0);
   
   var_name = (char *)malloc(strlen(file_name) - 2); assert(var_name!=NULL);
   strncpy(var_name, file_name, strlen(file_name) - 3);
   var_name[strlen(file_name) - 3] = 0;
   
   rtnflag = openInputDataNCfile(ncidP, directory, file_name, var_name);
   free(var_name);
   return(rtnflag);
   
   } // end of openInputDataNCfile() where var name is the same as the file name


bool Simulation::openInputDataNCfile(ncid_type * ncidP, char * directory, const char * file_name, const char * var_name)
   {
   char * file_name_and_path;
   char suffix[4];
   int rtnval;
   bool rtnflag;
   
   if (!(ncidP->fileid==-1))
      assert(false); // File should not have already been opened.
   ncidP->row_dimid = ncidP->col_dimid = ncidP->time_dimid = ncidP->ws_dimid = -1;
   ncidP->rowNdxPos = ncidP->colNdxPos = ncidP->timeNdxPos = ncidP->wsNdxPos = -1;
   
   if (file_name==NULL || strlen(file_name)==0)
      {
      printf("*** No file has been identified from which to read variable '%s'.\n", var_name);
      return(false);
      }
      
   strcpy(suffix, file_name + (strlen(file_name) - 3));
   assert(strcmp(suffix, ".nc")==0);
   
   // file_name_and_path is never free'd
   file_name_and_path = (char *)malloc(strlen(directory) + strlen(file_name) + 2); assert(file_name_and_path!=NULL);
   strcpy(file_name_and_path, directory);
   strcat(file_name_and_path, "/");
   strcat(file_name_and_path, file_name);
   ncidP->file_name_and_path = file_name_and_path; 
   
   assert(strlen(var_name)>0);
   
   // ncidP->var_name is never free'd
   ncidP->var_name = (char *)malloc(strlen(var_name) + 1); assert(ncidP->var_name!=NULL);
   strcpy(ncidP->var_name, var_name);
   
   rtnval = nc_open(file_name_and_path, NC_NOWRITE, &(ncidP->fileid)); 
   if (chk_nc(rtnval)) 
      { int varndims = 0;
      rtnflag = chk_nc(nc_inq_varid(ncidP->fileid, var_name, &(ncidP->input_data_varid)));
      rtnflag &= chk_nc(nc_inq_varndims(ncidP->fileid, ncidP->input_data_varid, &varndims));
      if (rtnflag) assert(varndims>=2 && varndims<=3);
      rtnflag &= chk_nc(nc_inq_vardimid(ncidP->fileid, ncidP->input_data_varid, ncidP->vardimids)); 
      rtnflag &= getRowInfo(ncidP);
      rtnflag &= getColInfo(ncidP);
      if (rtnflag) 
         {
         getOffsets(ncidP);
         getTimeInfo(ncidP); // Not every input variable has a time coordinate.
         getVarInfo(ncidP);
         }
      ncidP->rowNdxPos = ncidP->colNdxPos = ncidP->timeNdxPos = -1;
      for (int idim = 0; idim<varndims; idim++)
         {
         if (ncidP->vardimids[idim]==ncidP->row_dimid) ncidP->rowNdxPos = idim;
         else if (ncidP->vardimids[idim]==ncidP->col_dimid) ncidP->colNdxPos = idim;
         else if (ncidP->vardimids[idim]==ncidP->time_dimid) ncidP->timeNdxPos = idim;
         else assert(0);
         }
      }
   else rtnflag = false;
   if (!rtnflag) ncidP->fileid = -1;
   
   if (rtnflag) 
      {
      printf("Opened %s successfully.\n", file_name_and_path);
      // if (ncidP->row_offset==-9999) { printf("*** no row_offset attribute for %s\n", ncidP->var_name); ncidP->row_offset = 0; }
      // if (ncidP->col_offset==-9999) { printf("*** no col_offset attribute for %s\n", ncidP->var_name); ncidP->col_offset = 0; }
      if (ncidP->scale_factor==0.f) { printf("*** no scale_factor attribute for %s\n", ncidP->var_name); ncidP->scale_factor = 1.0f; }
      }
   else printf("Unable to open and interpret %s.\n", file_name_and_path);
   
   return(rtnflag);  
   } // end of openInputDataNCfile(ncid_type * ncidP, char * directory, char * file_name)


bool Simulation::openWarmstartNCfile(ncid_type * ncidP, char * file_name_and_path)
   {
   char suffix[4];
   int rtnval;
   bool rtnflag;
   static char dummy_var_name[] = "dummy_var_name";
   
   assert(file_name_and_path!=NULL);
   ncidP->file_name_and_path = file_name_and_path; 
   
   assert(strlen(file_name_and_path)>3);
   strcpy(suffix, file_name_and_path + (strlen(file_name_and_path) - 3));
   assert(strcmp(suffix, ".nc")==0);
   
   assert(ncidP->fileid==-1);
   rtnval = nc_open(file_name_and_path, NC_NOWRITE, &(ncidP->fileid)); 
   if (chk_nc(rtnval)) 
      { 
      rtnflag = getRowInfo(ncidP);
      rtnflag &= getColInfo(ncidP);
      rtnflag &= getOffsets(ncidP);
      getTimeInfo(ncidP);
      setNdxPos(ncidP);
      ncidP->scale_factor_flag = false;
      ncidP->scale_factor = 1.0f;
      ncidP->var_name = &(dummy_var_name[0]);
      ncidP->missing_value_flag = false;
      ncidP->fill_value = NC_FILL_FLOAT;
      ncidP->fill_value_flag = true;
      ncidP->fill_value = NC_FILL_FLOAT;
      }
   else rtnflag = false;
   if (!rtnflag) ncidP->fileid = -1;
   
   if (rtnflag) 
      {
      printf("Opened %s successfully.\n", file_name_and_path);
      
      switch (runParams.runMode)
         {
         case CENTURY_EQ:
            initializeMAPSSoutputList();
            for (int i = 0; i<NUM_MAPSS_OUTVARS; i++)
               { bool local_rtnflag;
               local_rtnflag = chk_nc(nc_inq_varid(ncidP->fileid, MAPSSoutputList[i].name, &MAPSSoutputList[i].var_id));
               if (!local_rtnflag)
                  {
                  printf("*** openWarmstartNCfile(): variable %i %s is missing from the _mapss.nc file.\n", i, MAPSSoutputList[i].name);
                  rtnflag = false;
                  }
               } // end of loop thru the MAPSS variables
            if (!rtnflag) 
                  printf("*** openWarmstartNCfile(): one or more variables are missing from the _mapss.nc file.\n");
            break;      
         case SPINUP:
            initializeEQoutputList(EQandWSinputList, NUM_EQ_OUTVARS);
            for (int i = 0; i<NUM_EQ_OUTVARS; i++) if (EQandWSinputList[i].show)
               rtnflag &= chk_nc(nc_inq_varid(ncidP->fileid, EQandWSinputList[i].name, &EQandWSinputList[i].var_id));
            if (!rtnflag) printf("*** openWarmstartNCfile(): one or more variables are missing from the _eq.nc file.\n");
            break;      
         case TRANSIENT:
            initializeWSoutputList(EQandWSinputList, NUM_WS_OUTVARS);
            ncidP->ws_varid = -1;
            if (modelParams.code_flags[FAST_WARMSTART_IN_FLAG])
              {
              rtnflag = false;
              rtnval = nc_inq_varid(ncidP->fileid, "ws_data", &(ncidP->ws_varid));
              if (rtnval==NC_NOERR)
                {
                rtnval = nc_inq_dimid(ncidP->fileid, "ws_ndx", &(ncidP->ws_dimid));
                if (rtnval==NC_NOERR) 
                  {
                  ncidP->wsNdxPos = 3;
                  rtnflag = true;
                  }
                else
                  {
                  ncidP->ws_dimid = -1;
                  ncidP->wsNdxPos = -1;
                  ncidP->ws_varid = -1;
                  }
                } 
              }
            else 
              {
              rtnflag = true;
              for (int i = 0; i<NUM_WS_OUTVARS; i++) if (EQandWSinputList[i].show)
                 rtnflag &= chk_nc(nc_inq_varid(ncidP->fileid, EQandWSinputList[i].name, &EQandWSinputList[i].var_id));
              }
            if (!rtnflag) printf("*** openWarmstartNCfile(): one or more variables are missing from the _ws.nc file.\n");
            break;      
         case MAPSS_EQ:
         default:
            assert(0);
            break;
         }
      }
   else printf("Unable to open and interpret %s.\n", file_name_and_path);
   
   return(rtnflag);  
   } // end of openWarmstartNCfile(ncid_type * ncidP, char * file_name_and_path)


void Simulation::getVarInfo(ncid_type * ncidP)
   {
   int rtnval;
   float scale_factor, scaled;
   
   /* Get the variable's type. */  
   rtnval = nc_inq_vartype(ncidP->fileid, ncidP->input_data_varid, &(ncidP->input_data_vartype)); assert(chk_nc(rtnval));
   
   
   /* Get the scale factor for the variable. */
   rtnval = nc_get_att_float(ncidP->fileid, ncidP->input_data_varid, "scale_factor", &scale_factor);
   if (rtnval==NC_NOERR)
      {
      assert(scale_factor!=0.0f);
      ncidP->scale_factor_flag = true;
      }
   else 
      {
      rtnval = nc_get_att_float(ncidP->fileid, ncidP->input_data_varid, "scaled", &scaled);
      if (rtnval==NC_NOERR) 
         {
         assert(scaled!=0.f);
         scale_factor = 1.f/scaled;
         ncidP->scale_factor_flag = true;
         }
      else ncidP->scale_factor_flag = false;
      }
   
   ncidP->scale_factor = ncidP->scale_factor_flag ? scale_factor : 1.0f;
      
   /* Get the fill value and/or missing data value, if they exist. */
   getAttInfo("missing_value", ncidP, &(ncidP->missing_value_flag), &(ncidP->missing_value));
   getAttInfo("_FillValue", ncidP, &(ncidP->fill_value_flag), &(ncidP->fill_value));
   
   } // end of getVarInfo()


void Simulation::getAttInfo(const char * att_name, ncid_type * pNcid, bool * pFlag, float * pValue) 
   {
   int rtnval;
   nc_type att_type;
   int intval;
   short shortval;
   signed char byteVal;
   double doubleVal;
   
   *pValue = 0.0;
   *pFlag = false;
   
   rtnval = nc_inq_atttype(pNcid->fileid, pNcid->input_data_varid, att_name, &att_type);
   if (rtnval!=NC_NOERR) return;
   
   switch (att_type)
      {
      case NC_FLOAT:
         rtnval = nc_get_att_float(pNcid->fileid, pNcid->input_data_varid, att_name, pValue);
         assert(chk_nc(rtnval));
         break;
      case NC_DOUBLE:
         rtnval = nc_get_att_double(pNcid->fileid, pNcid->input_data_varid, att_name, &doubleVal);
         assert(chk_nc(rtnval));
         *pValue = (float)doubleVal;
         break;
      case NC_BYTE:
         rtnval = nc_get_att_schar(pNcid->fileid, pNcid->input_data_varid, att_name, &byteVal);
         assert(chk_nc(rtnval));
         *pValue = (float)byteVal;
         break;
      case NC_SHORT:
         rtnval = nc_get_att_short(pNcid->fileid, pNcid->input_data_varid, att_name, &shortval);
         assert(chk_nc(rtnval));
         *pValue = (float)shortval;
         break;
      case NC_INT:
         rtnval = nc_get_att_int(pNcid->fileid, pNcid->input_data_varid, att_name, &intval);
         assert(chk_nc(rtnval));
         *pValue = (float)intval;
         break;
      default: assert(false); break;
      }
   *pFlag = true;
   
   return;
   } // end of Simulation::getAttInfo()
  

bool Simulation::getRowInfo(ncid_type * ncidP)
   {
   int rtnval;
   
	rtnval = nc_inq_dimid(ncidP->fileid, "lat", &(ncidP->row_dimid));
	if (rtnval!=NC_NOERR) 
      {
      rtnval = nc_inq_dimid(ncidP->fileid, "row", &(ncidP->row_dimid));
      if (rtnval!=NC_NOERR) 
         {
         ncidP->row_dimid = -1;
         ncidP->rowNdxPos = -1;
         ncidP->row_dimlen = 0;
         ncidP->row_offset = -9999;
         ncidP->row_varid = -1;
         return(false);
         }
      else 
         { // coordinate is called "row"
         rtnval = nc_inq_varid(ncidP->fileid, "row", &(ncidP->row_varid));
         if (rtnval!=NC_NOERR) ncidP->row_varid = -1;
         }
      }
   else 
      { // coordinate is called "lat"
      rtnval = nc_inq_varid(ncidP->fileid, "lat", &(ncidP->row_varid));
      if (rtnval!=NC_NOERR) ncidP->row_varid = -1;
      }
   
   rtnval = nc_inq_dimlen(ncidP->fileid, ncidP->row_dimid, &(ncidP->row_dimlen)); assert(chk_nc(rtnval));
   if (ncidP->row_varid!=-1)
      {
      rtnval = nc_inq_vartype(ncidP->fileid, ncidP->row_varid, &(ncidP->row_vartype)); assert(chk_nc(rtnval));
      }
/*
   rtnval = nc_inq_att(ncidP->fileid, NC_GLOBAL, "row_offset", &offset_type, &offset_length);
   if (rtnval==NC_NOERR)
      { int intval; float floatval;
      assert(offset_length==1);
      switch (offset_type)
         {
         case NC_SHORT: 
            rtnval = nc_get_att_short(ncidP->fileid, NC_GLOBAL, "row_offset", &(ncidP->row_offset)); assert(chk_nc(rtnval));
            break;
         case NC_INT: 
            rtnval = nc_get_att_int(ncidP->fileid, NC_GLOBAL, "row_offset", &intval); assert(chk_nc(rtnval));
            ncidP->row_offset = (short)intval;
            break;
         case NC_FLOAT:
            rtnval = nc_get_att_float(ncidP->fileid, NC_GLOBAL, "row_offset", &floatval); assert(chk_nc(rtnval));
            ncidP->row_offset = (short)floatval;
            break;
         default: assert(0); break;
         }
      }
   else ncidP->row_offset = -9999;
*/   
   return(true);
   } // end of getRowInfo()
   

bool Simulation::getOffsets(ncid_type * ncidP)
   {
   char row_offset_str[] = "row_offset";
   char col_offset_str[] = "col_offset";
   bool rtnFlag;
   
   rtnFlag = getAttVal(ncidP, row_offset_str, &(ncidP->row_offset));
   rtnFlag &= getAttVal(ncidP, col_offset_str, &(ncidP->col_offset));
   
   if (!rtnFlag)
      {
      ncidP->row_offset = ncidP->col_offset = 0;
      printf("*** getOffsets(): For file %s, one or both of the row and column offset attributes are missing or in an unexpected format.\n"
            "row and column offsets have been set to 0.\n", ncidP->file_name_and_path);
      }
      
   return(rtnFlag);
   
   } // end of Simulation::getOffsets()


bool Simulation::getAttVal(ncid_type * ncidP, const char * att_name, short * short_attValP)
   {
   nc_type att_type;
   size_t att_length;
   int rtnval;
   
   rtnval = nc_inq_att(ncidP->fileid, NC_GLOBAL, att_name, &att_type, &att_length);
   if (rtnval==NC_NOERR)
      { int intval; float floatval;
      assert(att_length==1);
      switch (att_type)
         {
         case NC_SHORT: 
            rtnval = nc_get_att_short(ncidP->fileid, NC_GLOBAL, att_name, short_attValP); assert(chk_nc(rtnval));
            break;
         case NC_INT: 
            rtnval = nc_get_att_int(ncidP->fileid, NC_GLOBAL, att_name, &intval); assert(chk_nc(rtnval));
            *short_attValP = (short)intval;
            break;
         case NC_FLOAT:
            rtnval = nc_get_att_float(ncidP->fileid, NC_GLOBAL, att_name, &floatval); assert(chk_nc(rtnval));
            *short_attValP = (short)floatval;
            break;
         default: assert(0); break;
         }
      }
      
   return(rtnval==NC_NOERR);
   
   } // end of Simulation::getAttVal()
   

bool Simulation::getColInfo(ncid_type * ncidP)
{
  int rtnval;
  
	rtnval = nc_inq_dimid(ncidP->fileid, "lon", &(ncidP->col_dimid));
	if (rtnval!=NC_NOERR) 
	{
	  rtnval = nc_inq_dimid(ncidP->fileid, "col", &(ncidP->col_dimid));
	  if (rtnval!=NC_NOERR) 
    {
      ncidP->col_dimid = -1;
      ncidP->colNdxPos = -1;
      ncidP->col_dimlen = 0;
      ncidP->col_offset = -9999;
      ncidP->col_varid = -1;
      return(false);
    }
    else 
    { // coordinate is called "col"
      rtnval = nc_inq_varid(ncidP->fileid, "col", &(ncidP->col_varid));
      if (rtnval!=NC_NOERR) ncidP->col_varid = -1;
    }
	}
  else 
  { // coordinate is called "lon"
    rtnval = nc_inq_varid(ncidP->fileid, "lon", &(ncidP->col_varid));
    if (rtnval!=NC_NOERR) ncidP->col_varid = -1;
  }
  
  rtnval = nc_inq_dimlen(ncidP->fileid, ncidP->col_dimid, &(ncidP->col_dimlen)); assert(chk_nc(rtnval));
  if (ncidP->col_varid!=-1)
  {
    rtnval = nc_inq_vartype(ncidP->fileid, ncidP->col_varid, &(ncidP->col_vartype)); assert(chk_nc(rtnval));
  }
/*
  rtnval = nc_inq_att(ncidP->fileid, NC_GLOBAL, "col_offset", &offset_type, &offset_length);
  if (rtnval==NC_NOERR)
  { int intval;
    assert(offset_length==1);
    if (offset_type==NC_SHORT) {rtnval = nc_get_att_short(ncidP->fileid, NC_GLOBAL, "col_offset", &(ncidP->col_offset)); assert(chk_nc(rtnval));}
    else if (offset_type==NC_INT)
    {
      rtnval = nc_get_att(ncidP->fileid, NC_GLOBAL, "col_offset", &intval); assert(chk_nc(rtnval));
      ncidP->col_offset = (short)intval;
    }
    else assert(0);
  }
  else ncidP->col_offset = -9999;
*/  
  return(true);
} // end of getColInfo()


bool Simulation::getTimeInfo(ncid_type * ncidP)
   {
   int rtnval;
  
   rtnval = nc_inq_dimid(ncidP->fileid, "month", &(ncidP->time_dimid));
	if (rtnval==NC_NOERR) 
      { // coordinate is called "month"
      rtnval = nc_inq_varid(ncidP->fileid, "month", &(ncidP->time_varid));
      if (rtnval!=NC_NOERR) ncidP->time_varid = -1;
      }
   else 
      {
      rtnval = nc_inq_dimid(ncidP->fileid, "year", &(ncidP->time_dimid));
      if (rtnval==NC_NOERR)
         { // coordinate is called "year"
         rtnval = nc_inq_varid(ncidP->fileid, "year", &(ncidP->time_varid));
         if (rtnval!=NC_NOERR) ncidP->time_varid = -1;
         }
      else 
         {
         rtnval = nc_inq_dimid(ncidP->fileid, "band", &(ncidP->time_dimid));
         if (rtnval==NC_NOERR)
            { // coordinate is called "band"
            rtnval = nc_inq_varid(ncidP->fileid, "band", &(ncidP->time_varid));
            if (rtnval!=NC_NOERR) ncidP->time_varid = -1;
            }
         else 
            {
            rtnval = nc_inq_dimid(ncidP->fileid, "time", &(ncidP->time_dimid));
            if (rtnval==NC_NOERR)
               { // coordinate is called "time"
               rtnval = nc_inq_varid(ncidP->fileid, "time", &(ncidP->time_varid));
               if (rtnval!=NC_NOERR) ncidP->time_varid = -1;
               }
            else  
               {
               ncidP->time_dimid = -1;
               ncidP->timeNdxPos = -1;
               ncidP->time_dimlen = 0;
               ncidP->time_varid = -1;
               return(false);
               }
            }
         }
      }
  
   rtnval = nc_inq_dimlen(ncidP->fileid, ncidP->time_dimid, &(ncidP->time_dimlen)); assert(chk_nc(rtnval));
   if (ncidP->time_varid!=-1)
      {
      rtnval = nc_inq_vartype(ncidP->fileid, ncidP->time_varid, &(ncidP->time_vartype)); assert(chk_nc(rtnval));
      }
  
   return(true);
   } // end of getTimeInfo()


bool Simulation::openInputFiles()
   {
   bool rtnflag, file_open_rtn_flag;
   unsigned int num_months;
   bool time_len_flag;
   
   assert(runParams.earth_data_directory!=NULL && runParams.soil_data_file!=NULL 
        && runParams.climate_data_directory!=NULL && runParams.soil_bulk_density_file!=NULL);
   
   file_open_rtn_flag = openInputDataNCfile(&mask_ncid, runParams.earth_data_directory, runParams.mask_file, "mask");
   // Next line is for backward compatibility with MC1...
   if (!file_open_rtn_flag) openInputDataNCfile(&mask_ncid, runParams.climate_data_directory, runParams.mask_file, "mask");
   
   rtnflag = openInputDataNCfile(&soilData_ncid, runParams.earth_data_directory, runParams.soil_data_file, "soils");   
   if (!rtnflag) rtnflag = openInputDataNCfile(&soilData_ncid, runParams.earth_data_directory, runParams.soil_data_file, "soils_scs"); // for backward compatibility with old soil data files   
   rtnflag &= runParams.runMode==MAPSS_EQ || openInputDataNCfile(&bd_ncid, runParams.earth_data_directory, runParams.soil_bulk_density_file, "bd");
   
   rtnflag &= openInputDataNCfile(&elev_ncid, runParams.earth_data_directory, "elev.nc");   
   rtnflag &= openInputDataNCfile(&ppt_ncid, runParams.climate_data_directory, "ppt.nc");
   num_months = ppt_ncid.time_dimlen; time_len_flag = num_months>=12;
   
   file_open_rtn_flag = openInputDataNCfile(&tmin_ncid, runParams.climate_data_directory, "tmin.nc");
   if (file_open_rtn_flag) time_len_flag &= tmin_ncid.time_dimlen==num_months;
   rtnflag &= file_open_rtn_flag || runParams.runMode==MAPSS_EQ; 
   
   file_open_rtn_flag = openInputDataNCfile(&tmax_ncid, runParams.climate_data_directory, "tmax.nc");
   if (file_open_rtn_flag) time_len_flag &= tmax_ncid.time_dimlen==num_months;
   rtnflag &= file_open_rtn_flag || runParams.runMode==MAPSS_EQ; 
   
   file_open_rtn_flag = openInputDataNCfile(&tmp_ncid, runParams.climate_data_directory, "tmp.nc");
   if (file_open_rtn_flag) time_len_flag &= tmp_ncid.time_dimlen==num_months;
   rtnflag &= file_open_rtn_flag || (tmax_ncid.fileid!=-1 && tmin_ncid.fileid!=-1); 
   
   file_open_rtn_flag = openInputDataNCfile(&vpr_ncid, runParams.climate_data_directory, "vpr.nc");
   if (file_open_rtn_flag) time_len_flag &= vpr_ncid.time_dimlen==num_months;
   else 
      {
      file_open_rtn_flag = openInputDataNCfile(&tdmean_ncid, runParams.climate_data_directory, "tdmean.nc");
      if (file_open_rtn_flag) time_len_flag &= tdmean_ncid.time_dimlen==num_months;
      }  
   rtnflag &= file_open_rtn_flag; 
   
   if (!runParams.dummyWindFlag)
      {
      // Look first in the climate data directory, then in the earth data directory.   
      file_open_rtn_flag = openInputDataNCfile(&wnd_ncid, runParams.climate_data_directory, "wnd.nc");
      if (!file_open_rtn_flag) file_open_rtn_flag = openInputDataNCfile(&wnd_ncid, runParams.earth_data_directory, "wnd.nc");
      if (file_open_rtn_flag) time_len_flag &= wnd_ncid.time_dimlen==12 || wnd_ncid.time_dimlen==num_months;
      else printf("*** openInputFiles(): Unable to open wnd.nc.  If you wish to use dummy "
            "wind data, put \'dummy_wind_switch = ON\' in the command file.\n");
      rtnflag &= file_open_rtn_flag;
      }
      
   if (runParams.baseCalibration==mc2WCR)
      { // Open files for fog factor, etc.
      file_open_rtn_flag = openInputDataNCfile(&fog_ncid, runParams.earth_data_directory, "fog.nc");
      rtnflag &= file_open_rtn_flag;
      assert(fog_ncid.timeNdxPos<0);
      
      file_open_rtn_flag = openInputDataNCfile(&topomoist_ncid, runParams.earth_data_directory, "topomoist.nc");
      rtnflag &= file_open_rtn_flag;
      assert(topomoist_ncid.timeNdxPos<0);
      
      file_open_rtn_flag = openInputDataNCfile(&deltaTsl_ncid, runParams.earth_data_directory, "deltaTsl.nc");
      rtnflag &= file_open_rtn_flag;
      assert(deltaTsl_ncid.timeNdxPos<0);
      
      file_open_rtn_flag = openInputDataNCfile(&aspect_ncid, runParams.earth_data_directory, "aspect.nc");
      rtnflag &= file_open_rtn_flag;
      assert(aspect_ncid.timeNdxPos<0);
      
      file_open_rtn_flag = openInputDataNCfile(&sw_ncid, runParams.earth_data_directory, "shortwave.nc");
      rtnflag &= file_open_rtn_flag;
      assert(sw_ncid.timeNdxPos<0);
      
      file_open_rtn_flag = openInputDataNCfile(&cad_ncid, runParams.earth_data_directory, "cad.nc");
      rtnflag &= file_open_rtn_flag;
      assert(cad_ncid.timeNdxPos<0);
      }
   
   if (time_len_flag) runParams.years_of_climate_data = num_months/12;
   else 
      {
      rtnflag = false;
      runParams.years_of_climate_data = 0;
      printf("*** openInputFiles(): the number of months of climate data is less than 12 or is not the same in all "
             "the climate data input files.\n");
      }
   // assert(runParams.years_of_climate_data>=1);
   if (!(runParams.years_of_climate_data>=1))
      assert(0);
   
   // Open the warmstart file, if appropriate.
   switch (runParams.runMode)
      {
      case MAPSS_EQ: // no warmstart file for MAPSS phase
         break; 
      
      case CENTURY_EQ: // open the file created by the MAPSS phase
         // the path and file name are in runParams.warmstart_file
         file_open_rtn_flag = openWarmstartNCfile(&MAPSSdata_ncid, runParams.warmstart_file); 
         rtnflag &= file_open_rtn_flag;
         break;
      
      case SPINUP: // open the file containing the MC2centuryEQdataOutputs object created by the CENTURY_EQ phase
         // the path and file name are in runParams.warmstart_file
         file_open_rtn_flag = openWarmstartNCfile(&EQdata_ncid, runParams.warmstart_file); 
         rtnflag &= file_open_rtn_flag;
         break;
      
      case TRANSIENT:
         // read a file containing the StateVariables object created by a spinup or transient phase
         file_open_rtn_flag = openWarmstartNCfile(&WSdata_ncid, runParams.warmstart_file); 
         rtnflag &= file_open_rtn_flag;
         break;
      
      case unknownRunMode:
      default: assert(0); break;
      }
   
   return(rtnflag);
   } // end of openInputFiles()


void instructions()
{
  printf("*** " COMMAND_NAME ": Invoke as '" COMMAND_NAME " <calibration> <command file>\n"
      "  where 'calibration' is one of: GLOBAL, CONUS, CONUS_LC, WCR\n"
	  "  Example: " COMMAND_NAME " GLOBAL runGlobal_Spinup.txt\n\n");
  exit(-1);
} // end of instructions()


void Simulation::initializeOutputFiles()
   {
   char * cmd_line = m_cmdLine;
   char * cmd_file_name = m_cmdFileName;
   
   switch (runParams.runMode)
      {
      case MAPSS_EQ: initializeMAPSSoutputFile(cmd_line, cmd_file_name); break;
      case CENTURY_EQ: initializeEQoutputFile(cmd_line, cmd_file_name); break;
      case SPINUP: 
         runParams.warmstartOutputFlag = true;
      case TRANSIENT:
         if (runParams.warmstartOutputFlag) initializeWSoutputFile(cmd_line, cmd_file_name);
         initializeDiscretionaryOutputFiles(cmd_line, cmd_file_name);
         break;
      default: assert(0); break;
      }
   
   } // end of Simulation::initializeOutputFiles()


void Simulation::initializeMAPSSoutputFile(char * cmd_line, char * cmd_file_name)
   {
   int rtnval;
   bool rtnFlag;
   char * file_name;
   int dimids[3];
   
   file_name = makePath(runParams.output_file_prefix, "_mapss.nc");
   rtnval = nc_create(file_name, NC_CLOBBER, &(MAPSSoutFile_ncid.fileid)); assert(chk_nc(rtnval));
   free(file_name);
   writeStdAttributes(&MAPSSoutFile_ncid);
   writeCommands(&MAPSSoutFile_ncid, cmd_line, cmd_file_name);
   writeRunParameters(&MAPSSoutFile_ncid);
   writeModelParameters(&MAPSSoutFile_ncid);
   rtnFlag = getRowInfo(&MAPSSoutFile_ncid); assert(rtnFlag);
   rtnFlag = getColInfo(&MAPSSoutFile_ncid); assert(rtnFlag);
   rtnFlag = getOffsets(&MAPSSoutFile_ncid); assert(rtnFlag);
   rtnFlag = getTimeInfo(&MAPSSoutFile_ncid); assert(rtnFlag); 
   setNdxPos(&MAPSSoutFile_ncid);

   dimids[0] = MAPSSoutFile_ncid.time_dimid;    
   dimids[1] = MAPSSoutFile_ncid.row_dimid;
   dimids[2] = MAPSSoutFile_ncid.col_dimid;
   
   for (int i = 0; i<NUM_MAPSS_OUTVARS; i++) MAPSSoutputList[i].show = false;
   MAPSSoutputList[MAPSSmclass] = describeVar("mclass", "MAPSS vegetation class", "", NC_SHORT);
   MAPSSoutputList[MAPSSvclass] = describeVar("vclass", "VEMAP vegetation class", "", NC_SHORT);
   MAPSSoutputList[MAPSScanopy] = describeVar("canopy", "1=C3Dominance, 2=C3C4Mixed, 3=C4Dominance", "", NC_FLOAT);
   MAPSSoutputList[MAPSSzone] = describeVar("MAPSS_climate_zone", "0-3: boreal, temperate, subtropical, tropical", "zone", NC_SHORT);
   
   defineNCvarsInList(&MAPSSoutFile_ncid, MAPSSoutputList, NUM_MAPSS_OUTVARS);
   
   rtnval = nc_enddef(MAPSSoutFile_ncid.fileid); assert(chk_nc(rtnval));

   writeOutputVarList("MAPSS_output_vars.txt", MAPSSoutputList, NUM_MAPSS_OUTVARS);
   
   } // end of Simulation::initializeMAPSSoutputFile()


void Simulation::makeOutputListEntry(OutVarDescription * f, const char * name, const char * description, const char * units, nc_type type)
   {
   f->name     = name;
   f->descrip  = description;
   f->units    = units;
   f->type     = type;
   f->scale    = 1.000000 ;
   f->interval = 0;
   } // end of Simulation::makeOutputListEntry()


void Simulation::initializeMAPSSoutputList()
   {
   for (int i = 0; i<NUM_MAPSS_OUTVARS; i++) MAPSSoutputList[i].show = false;
   makeOutputListEntry(&MAPSSoutputList[MAPSSmclass], "mclass", "MAPSS vegetation class", "", NC_SHORT);
   makeOutputListEntry(&MAPSSoutputList[MAPSSvclass], "vclass", "VEMAP vegetation class", "", NC_SHORT);
   makeOutputListEntry(&MAPSSoutputList[MAPSScanopy], "canopy","1=C3Dominance, 2=C3C4Mixed, 3=C4Dominance", "%", NC_FLOAT);
   makeOutputListEntry(&MAPSSoutputList[MAPSSzone], "MAPSS_climate_zone","0-3: boreal, temperate, subtropical, tropical", "zone", NC_SHORT);
   } // end of Simulation::initializeMAPSSoutputList()


void Simulation::initializeEQoutputFile(char * cmd_line, char * cmd_file_name)
   {
   int rtnval;
   bool rtnFlag;
   char * file_name;
   
   file_name = makePath(runParams.output_file_prefix, "_eq.nc");
   rtnval = nc_create(file_name, NC_CLOBBER, &(EQoutFile_ncid.fileid)); assert(chk_nc(rtnval));
   free(file_name);
   writeStdAttributes(&EQoutFile_ncid);
   writeCommands(&EQoutFile_ncid, cmd_line, cmd_file_name);
   writeRunParameters(&EQoutFile_ncid);
   writeModelParameters(&EQoutFile_ncid);
   rtnFlag = getRowInfo(&EQoutFile_ncid); assert(rtnFlag);
   rtnFlag = getColInfo(&EQoutFile_ncid); assert(rtnFlag);
   rtnFlag = getOffsets(&EQoutFile_ncid); assert(rtnFlag);
   rtnFlag = getTimeInfo(&EQoutFile_ncid); assert(rtnFlag); 
   setNdxPos(&EQoutFile_ncid);
   
   initializeEQoutputList(EQandWSoutputList, sizeof(EQandWSoutputList)/sizeof(OutVarDescription));   
   defineNCvarsInList(&EQoutFile_ncid, EQandWSoutputList, sizeof(EQandWSoutputList)/sizeof(OutVarDescription));
   rtnval = nc_enddef(EQoutFile_ncid.fileid); assert(chk_nc(rtnval));

   writeOutputVarList("EQ_output_vars.txt", EQandWSoutputList, sizeof(EQandWSoutputList)/sizeof(OutVarDescription));

   } // end of Simulation::initializeEQoutputFile()  
  

void Simulation::initializeWSoutputFile(char * cmd_line, char * cmd_file_name)
{
   int rtnval;
   bool rtnFlag;
   char * file_name;
   
   file_name = makePath(runParams.output_file_prefix, "_ws.nc");
   rtnval = nc_create(file_name, NC_CLOBBER, &(WSoutFile_ncid.fileid)); assert(chk_nc(rtnval));
   free(file_name);
   writeStdAttributes(&WSoutFile_ncid);
   writeCommands(&WSoutFile_ncid, cmd_line, cmd_file_name);
   writeRunParameters(&WSoutFile_ncid);
   writeModelParameters(&WSoutFile_ncid);
   rtnFlag = getRowInfo(&WSoutFile_ncid); assert(rtnFlag);
   rtnFlag = getColInfo(&WSoutFile_ncid); assert(rtnFlag);
   rtnFlag = getOffsets(&WSoutFile_ncid); assert(rtnFlag);
   rtnFlag = getTimeInfo(&WSoutFile_ncid); assert(rtnFlag); 
   setNdxPos(&WSoutFile_ncid);
   
   initializeWSoutputList(EQandWSoutputList, sizeof(EQandWSoutputList)/sizeof(OutVarDescription));
   
   if (modelParams.code_flags[FAST_WARMSTART_OUT_FLAG])
     {
     int dimids[4]; 
      
     WSoutFile_ncid.vardimids[3] = WSoutFile_ncid.ws_dimid = ncdimdef(WSoutFile_ncid.fileid, "ws_ndx", NUM_WS_OUTVARS);
     assert(WSoutFile_ncid.ws_dimid==3); WSoutFile_ncid.wsNdxPos = 3; 
     dimids[0] = WSoutFile_ncid.time_dimid;
     dimids[1] = WSoutFile_ncid.row_dimid;
     dimids[2] = WSoutFile_ncid.col_dimid;
     dimids[3] = WSoutFile_ncid.ws_dimid;
     WSoutFile_ncid.ws_varid = ncvardef(WSoutFile_ncid.fileid, "ws_data", NC_FLOAT, 4, dimids); 
     ncattput(WSoutFile_ncid.fileid, WSoutFile_ncid.ws_varid, "long_name", NC_CHAR, strlen("warmstart data array"), "warmstart data array");
     }
   else defineNCvarsInList(&WSoutFile_ncid, EQandWSoutputList, sizeof(EQandWSoutputList)/sizeof(OutVarDescription));

   rtnval = nc_enddef(WSoutFile_ncid.fileid); assert(chk_nc(rtnval));
   
   writeOutputVarList("WS_output_vars.txt", EQandWSoutputList, sizeof(EQandWSoutputList)/sizeof(OutVarDescription));

} // end of Simulation::initializeWSoutputFile()  


void Simulation::initializeDiscretionaryOutputFiles(char * cmd_line, char * cmd_file_name)
{
	FILE * fp;
#define NAME_LEN 100
	char buffer[NAME_LEN];
	char name[NAME_LEN];
  
   if (strlen(runParams.output_variable_file)==0 || strcmp(runParams.output_variable_file, "NONE")==0) 
         printf("No output variable file was specified.\n");
   else
      { // Process the output variable file.
      printf("Loading output variable file: %s\n", runParams.output_variable_file);
      fp = fopen(runParams.output_variable_file, "r"); 
      if (fp==NULL)
         {
         fprintf(stderr,"Can't find output filter file\n");
         assert(false);
         }

      while (!feof(fp)) 
         {
         fgets(buffer, NAME_LEN, fp);
         if (feof(fp)) break;
         memset((char *)name, '\0', sizeof(name));
         sscanf(buffer,"%s",name);

         if (strlen(name) == 0 || name[0] == '#') continue; // ignore lines w/ # and empty lines
         
         addVarToList(name);
         } // end of loop to read the output variable file
         
      fclose(fp);
      } // end of block to process the output variable file
  
   // Deal with monthly output data.
   if (m_num_monthly_vars>0)
      { 
      int moBufLen; 

      moBufLen = 12*runParams.actual_years_to_run*m_num_monthly_vars;
      
      // m_moBufP is never free'd
      m_moBufP = (float *)malloc(sizeof(float)*moBufLen); assert(m_moBufP!=NULL);
      
      if (runParams.monthlyOutputFlag)
         {
         int rtnval;
         size_t coords[1];
      
         initializeDiscretionaryOutputFile(&m_moOutFile_ncid, "_month.nc", m_moOutList, m_num_monthly_vars, cmd_line, cmd_file_name);

         // Write out values of time coordinate variable.
         rtnval = nc_enddef(m_moOutFile_ncid.fileid); assert(chk_nc(rtnval));
         for (int moNdx = 0; moNdx<12*runParams.actual_years_to_run; moNdx++)
            {
            float year = runParams.first_calendar_year_of_run + moNdx/12 + moNdx%12*(1./12.) + 1./24.;
            coords[0] = moNdx;
            rtnval = nc_put_var1_float(m_moOutFile_ncid.fileid, m_moOutFile_ncid.time_varid, coords, &year); assert(chk_nc(rtnval));
            }
         } // end of if (runParams.monthlyOutputFlag)
      } // end of if (m_num_monthly_vars>0)
   else m_moBufP = NULL;

   // Deal with yearly output data.
   if (m_num_yearly_vars>0)
      {
      int yrBufLen;
      
      yrBufLen = runParams.actual_years_to_run*m_num_yearly_vars;
      
      // m_yrBufP is never free'd
      m_yrBufP = (float *)malloc(sizeof(float)*yrBufLen); assert(m_yrBufP!=NULL);

      if (runParams.yearlyOutputFlag)
         { 
         int rtnval;
         size_t coords[1];

         initializeDiscretionaryOutputFile(&m_yrOutFile_ncid, "_year.nc", m_yrOutList, m_num_yearly_vars, cmd_line, cmd_file_name);
      
         // Write out values of time coordinate variable.
         rtnval = nc_enddef(m_yrOutFile_ncid.fileid); assert(chk_nc(rtnval));
         for (int yrNdx = 0; yrNdx<runParams.actual_years_to_run; yrNdx++)
            {
            float year = (float)(runParams.first_calendar_year_of_run + yrNdx);
            coords[0] = yrNdx;
            rtnval = nc_put_var1_float(m_yrOutFile_ncid.fileid, m_yrOutFile_ncid.time_varid, coords, &year); assert(chk_nc(rtnval));
            } // end of for (int yrNdx = 0; yrNdx<runParams.actual_years_to_run; yrNdx++)
         } // end of if (runParams.yearlyOutputFlag)
      } // end of if (m_num_yearly_vars>0)
   else m_yrBufP = NULL;
   
   // Deal with multiyr output data.
   if (m_num_multiyr_vars>0)
      { 
      int multiyrBufLen, rtnval, ending_year_var_id, interval_length_var_id;
      int dimids[1];
      size_t coords[1];
      
      multiyrBufLen = runParams.num_multiyr_intervals*m_num_multiyr_vars;

      // m_multiyrBufP and m_multiyrIntervalBufP are never free'd
      m_multiyrBufP = (float *)malloc(sizeof(float)*multiyrBufLen); assert(m_multiyrBufP!=NULL);
      m_multiyrIntervalBufP = (short *)malloc(sizeof(short)*2*runParams.num_multiyr_intervals); assert(m_multiyrBufP!=NULL);
      
      // Calculate the ending year and length of each interval
      int offset = 0;
      for (int interval_index = 0; interval_index<runParams.num_multiyr_intervals; interval_index++)
         { 
         int first_yr_of_interval, last_yr_of_interval, current_interval_len;
         
         if (interval_index==0 && runParams.runMode==SPINUP)
            {
            first_yr_of_interval = last_yr_of_interval = 0;
            current_interval_len = 1;
            offset = 1;
            }
         else
            {
            first_yr_of_interval = runParams.nominal_multiyr_start + (interval_index - offset)*runParams.multiyr_len;
            last_yr_of_interval = first_yr_of_interval + runParams.multiyr_len - 1;
            }
         if (first_yr_of_interval<runParams.first_calendar_year_of_run) 
               first_yr_of_interval = runParams.first_calendar_year_of_run;
         if (last_yr_of_interval>runParams.last_calendar_year_of_run) 
               last_yr_of_interval = runParams.last_calendar_year_of_run;
         current_interval_len = last_yr_of_interval - first_yr_of_interval + 1;
         assert(current_interval_len>=1);
         *(m_multiyrIntervalBufP + 2*interval_index) = (short)last_yr_of_interval;
         *(m_multiyrIntervalBufP + 2*interval_index + 1) = (short)current_interval_len;
         } // end of loop to calculate interval ending years and lengths
         
      if (runParams.multiyrOutputFlag)
         {
         initializeDiscretionaryOutputFile(&m_multiyrOutFile_ncid, "_multiyr.nc", m_multiyrOutList, m_num_multiyr_vars, cmd_line, cmd_file_name);
      
         // Define the "ending_year" and "interval_length" variables in the _multiyr.nc file.
         dimids[0] = m_multiyrOutFile_ncid.time_dimid;
         rtnval = nc_def_var(m_multiyrOutFile_ncid.fileid, "ending_year", NC_SHORT, 1, dimids, &(ending_year_var_id)); assert(chk_nc(rtnval));
         rtnval = nc_def_var(m_multiyrOutFile_ncid.fileid, "interval_length", NC_SHORT, 1, dimids, &(interval_length_var_id)); assert(chk_nc(rtnval));
         rtnval = nc_enddef(m_multiyrOutFile_ncid.fileid); assert(chk_nc(rtnval)); 
         
         // Write the data to the ending_year and interval_length variables.
         for (int interval_index = 0; interval_index<runParams.num_multiyr_intervals; interval_index++)
            {
            coords[0] = interval_index;
            rtnval = nc_put_var1_short(m_multiyrOutFile_ncid.fileid, ending_year_var_id, coords, m_multiyrIntervalBufP + 2*interval_index); assert(chk_nc(rtnval));
            float year = (float)(*(m_multiyrIntervalBufP + 2*interval_index));
            rtnval = nc_put_var1_float(m_multiyrOutFile_ncid.fileid, m_multiyrOutFile_ncid.time_varid, coords, &year); assert(chk_nc(rtnval));
            rtnval = nc_put_var1_short(m_multiyrOutFile_ncid.fileid, interval_length_var_id, coords, m_multiyrIntervalBufP + 2*interval_index + 1); assert(chk_nc(rtnval));
            } 
         } // end of block for if (runParams.multiyrOutputFlag)
      } // end of block for if (m_num_multiyr_vars>0)
   else 
      {
      m_multiyrBufP = NULL;
      m_multiyrIntervalBufP = NULL;
      }
   
   // Deal with time invariant output data.
   if (m_num_single_vars>0 && runParams.timeInvariantOutputFlag)
      { 
      int rtnval;
      
      initializeDiscretionaryOutputFile(&m_singleOutFile_ncid, "_single.nc", m_singleOutList, m_num_single_vars, cmd_line, cmd_file_name);      
      rtnval = nc_enddef(m_singleOutFile_ncid.fileid); assert(chk_nc(rtnval));
      }

} // end of Simulation::initializeDiscretionaryOutputFiles()  
   
   
void Simulation::initializeDiscretionaryOutputFile(ncid_type * ncidP, const char * suffix, OutVarDescription outList[], int num_vars,
      char * cmd_line, char * cmd_file_name)
   { 
   char * file_name;
   int rtnval;
   bool rtnFlag;
   
   assert(num_vars>0);
   
   file_name = makePath(runParams.output_file_prefix, suffix);
   rtnval = nc_create(file_name, NC_CLOBBER, &(ncidP->fileid)); assert(chk_nc(rtnval));
   free(file_name);
   writeStdAttributes(ncidP);
   writeCommands(ncidP, cmd_line, cmd_file_name);
   writeRunParameters(ncidP);
   writeModelParameters(ncidP);
   rtnFlag = getRowInfo(ncidP); assert(rtnFlag);
   rtnFlag = getColInfo(ncidP); assert(rtnFlag);
   rtnFlag = getOffsets(ncidP); assert(rtnFlag);
   rtnFlag = getTimeInfo(ncidP); assert(rtnFlag); 
   setNdxPos(ncidP);
   
   defineNCvarsInList(ncidP, outList, num_vars);
   
   } // end of Simulation::initializeDiscretionaryOutputFile()


void Simulation::setNdxPos(ncid_type * ncidP)
   { 
   int unlimdimid; 
   int rtnval;
   
   assert(ncidP->time_dimid==0); ncidP->timeNdxPos = 0; ncidP->vardimids[0] = 0;
   assert(ncidP->row_dimid==1); ncidP->rowNdxPos = 1; ncidP->vardimids[1] = 1;
   assert(ncidP->col_dimid==2); ncidP->colNdxPos = 2; ncidP->vardimids[2] = 2;
   rtnval = nc_inq_unlimdim(ncidP->fileid, &unlimdimid); assert(chk_nc(rtnval));
   assert(ncidP->time_dimid==unlimdimid);
   
   } // end of Simulation::setNdxPos()



void Simulation::defineNCvarsInList(ncid_type * ncidP, OutVarDescription * list, int filter_length)
   {
   int dimids[3];
  
   dimids[0] = ncidP->time_dimid;    
   dimids[1] = ncidP->row_dimid;
   dimids[2] = ncidP->col_dimid;
  
   for (int i = 0; i<filter_length; i++) if (list[i].show) list[i].var_id = 
         defineNCvar(ncidP->fileid, dimids, list[i].name, list[i].descrip, list[i].units, list[i].type);
         
   } // end of Simulation::defineNCvarsInList()

 
 void Simulation::writeStdAttributes(ncid_type * ncidP)
   {
   short nc_short;
   int rtnval;
   int dimids[3];
   float lat, lon;
   size_t coords[3];
   
   nc_short = runParams.row_begin;
   ncattput(ncidP->fileid, NC_GLOBAL, "row_offset", NC_SHORT, 1, &(nc_short) );
   nc_short = runParams.col_begin;
   ncattput(ncidP->fileid, NC_GLOBAL, "col_offset", NC_SHORT, 1, &(nc_short) );
   
   rtnval = nc_def_dim(ncidP->fileid, "year", NC_UNLIMITED, &(ncidP->time_dimid)); assert(chk_nc(rtnval));
   ncidP->row_dimid = ncdimdef(ncidP->fileid, "lat", runParams.row_end - runParams.row_begin + 1);
   ncidP->col_dimid = ncdimdef(ncidP->fileid, "lon", runParams.col_end - runParams.col_begin + 1);
   
	/* row coordinate variable */
	dimids[0] = ncidP->row_dimid;
	ncidP->row_varid = ncvardef(ncidP->fileid, "lat", NC_FLOAT, 1, dimids); 
   ncattput(ncidP->fileid, ncidP->row_varid, "long_name", NC_CHAR, strlen("latitude"), "latitude");
   ncattput(ncidP->fileid, ncidP->row_varid, "standard_name", NC_CHAR, strlen("latitude"), "latitude");
   ncattput(ncidP->fileid, ncidP->row_varid, "units", NC_CHAR, strlen("degrees_north"), "degrees_north");
   
	/* column coordinate variable */
	dimids[0] = ncidP->col_dimid;
	ncidP->col_varid = ncvardef(ncidP->fileid, "lon", NC_FLOAT, 1, dimids); 
   ncattput(ncidP->fileid, ncidP->col_varid, "long_name", NC_CHAR, strlen("longitude"), "longitude");
   ncattput(ncidP->fileid, ncidP->col_varid, "standard_name", NC_CHAR, strlen("longitude"), "longitude");
   ncattput(ncidP->fileid, ncidP->col_varid, "units", NC_CHAR, strlen("degrees_east"), "degrees_east");
   
	/* time coordinate variable */
   dimids[0] = ncidP->time_dimid;
   rtnval = nc_def_var(ncidP->fileid, "year", NC_FLOAT, 1, dimids, &(ncidP->time_varid)); assert(chk_nc(rtnval));
   ncattput(ncidP->fileid, ncidP->time_varid, "standard_name", NC_CHAR, strlen("year"), "year");
   ncattput(ncidP->fileid, ncidP->time_varid, "units", NC_CHAR, strlen("year"), "year");
   
   // write out the row and column coordinate variables
   rtnval = nc_enddef(ncidP->fileid); assert(chk_nc(rtnval));
   for (int row = runParams.row_begin; row<=runParams.row_end; row++)
   for (int col = runParams.col_begin; col<=runParams.col_end; col++)
      {
      lat = (float)(runParams.latitude0 - runParams.cell_spacing*row - runParams.cell_spacing/2.);
      lon = (float)(runParams.longitude0 + runParams.cell_spacing*col + runParams.cell_spacing/2.); 
      coords[0] = row - runParams.row_begin;
      rtnval = nc_put_var1_float(ncidP->fileid, ncidP->row_varid, coords, &lat); assert(chk_nc(rtnval));
      coords[0] = col - runParams.col_begin;
      rtnval = nc_put_var1_float(ncidP->fileid, ncidP->col_varid, coords, &lon); assert(chk_nc(rtnval));
      }
   rtnval = nc_redef(ncidP->fileid); assert(chk_nc(rtnval));
  
   } // end of Simulation::writeStdAttributes()


OutVarDescription Simulation::describeVar(const char * var_name, const char * var_descrip, const char * var_units, nc_type var_type)
   {
   OutVarDescription outVarRecord;
   
   outVarRecord.descrip = var_descrip;
   outVarRecord.type = var_type;
   outVarRecord.scale = 1.0;
   outVarRecord.name = var_name;
   outVarRecord.interval = MONTH_INTERVAL;
   outVarRecord.show = true; outVarRecord.m_pS = this;
   outVarRecord.units = var_units;
   // outVarRecord.var_id = defineNCvar(fileid, dimids, var_name, var_descrip, var_units, var_type);
   
   return(outVarRecord);
   
   } // end of Simulation::describeVar()


void Simulation::makeVarDictEntry(int varNdx, const char * var_name, const char * var_descrip, const char * var_units, nc_type var_type, int interval)
   {
   OutVarDescription outVarRecord = describeVar(var_name, var_descrip, var_units, var_type);

   outVarRecord.interval = interval;
   outVarRecord.show = false; outVarRecord.m_pS = NULL;
   outVarRecord.dictNdx = varNdx;
   
   printf("%s, %d, %d\n", var_name, varNdx, VarDict[varNdx].dictNdx);
   if (VarDict[varNdx].dictNdx!=-1 && VarDict[varNdx].dictNdx!=varNdx)
      assert(0);
   VarDict[varNdx] = outVarRecord;
   
   } // end of Simulation::makeVarDictEntry()


void Simulation::makeCategoricalVarDictEntry(int varNdx, const char * var_name, const char * var_descrip, const char * var_units, int max_categories)
   {
   makeVarDictEntry(varNdx, var_name, var_descrip, var_units, NC_SHORT, YEAR_INTERVAL);
   VarDict[varNdx].categoricalFlag = true;
   VarDict[varNdx].max_categories = max_categories;
   } // end of makeCategoricalVarDictEntry()


int Simulation::defineNCvar(int fileid, int dimids[], const char * var_name, const char * var_descrip, const char * var_units, nc_type var_type)
   {
   int var_id;
   int rtnval;
   int ndims = 3;
   short nc_fill_short = NC_FILL_SHORT;
   float nc_fill_float = NC_FILL_FLOAT;
   int nc_fill_int = NC_FILL_INT;
   
   switch (var_type)
      {
      case NC_SHORT:
         rtnval = nc_def_var(fileid, var_name, var_type, ndims, dimids, &(var_id)); assert(chk_nc(rtnval));
         rtnval = nc_put_att_short(fileid, var_id, "_FillValue", NC_SHORT, 1, &nc_fill_short); assert(chk_nc(rtnval));
         break;
      case NC_FLOAT:
         rtnval = nc_def_var(fileid, var_name, var_type, ndims, dimids, &(var_id)); if (!(chk_nc(rtnval)))
            assert(false);
         rtnval = nc_put_att_float(fileid, var_id, "_FillValue", NC_FLOAT, 1, &nc_fill_float); assert(chk_nc(rtnval));
         break;
      case NC_INT:
         rtnval = nc_def_var(fileid, var_name, var_type, ndims, dimids, &(var_id)); assert(chk_nc(rtnval));
         rtnval = nc_put_att_int(fileid, var_id, "_FillValue", NC_INT, 1, &nc_fill_int); assert(chk_nc(rtnval));
         break;
      default:
         printf("var_name = %s, var_type = %d\n", var_name, var_type);
         err_exit("Did not recognize the type of this output variable.");
         break;
      }
   
   ncattput(fileid, var_id, "long_name", NC_CHAR, strlen(var_descrip), var_descrip);
   ncattput(fileid, var_id, "units", NC_CHAR, strlen(var_units), var_units); 
   // ncattput(fileid, var_id, "scale_factor", NC_FLOAT, 1, &(var_scale)); 
   
   return(var_id);
   
   } // end of Simulation::defineVar()


char * Simulation::makePath(char * prefix, const char * suffix)
   {
   char * name;
   
   // name is never free'd
   name = (char *)malloc(strlen(prefix) + strlen(suffix) + 3); assert(name!=NULL);
   strcpy(name, prefix);
   strcat(name, suffix);
   return(name);
   } // end of Simulation::makePath()


void Simulation::saveOutputData()
   {
   if (runParams.runMode==MAPSS_EQ || runParams.runMode==CENTURY_EQ || runParams.warmstartOutputFlag)
   switch (runParams.runMode)
      {
      case MAPSS_EQ: saveOutputFile(&MAPSSoutFile_ncid); break;
      case CENTURY_EQ: saveOutputFile(&EQoutFile_ncid); break;
      case SPINUP: saveOutputFile(&WSoutFile_ncid); break;
      case TRANSIENT:
         if (runParams.warmstartOutputFlag) saveOutputFile(&WSoutFile_ncid); 
         break;
      default: assert(0); break;
      }

   if (m_num_monthly_vars>0 && runParams.monthlyOutputFlag) saveOutputFile(&m_moOutFile_ncid);
   if (m_num_yearly_vars>0 && runParams.yearlyOutputFlag) saveOutputFile(&m_yrOutFile_ncid);
   if (m_num_multiyr_vars>0 && runParams.multiyrOutputFlag) saveOutputFile(&m_multiyrOutFile_ncid);
   if (m_num_single_vars>0 && runParams.timeInvariantOutputFlag) saveOutputFile(&m_singleOutFile_ncid);
      
   if (m_num_vpr_tmp_issues>0) printf("*** saveOutputData: m_num_vpr_tmp_issues = %d\n", m_num_vpr_tmp_issues);

   } // end of Simulation::saveOutputData()


void Simulation::saveOutputFile(ncid_type * ncidP)
   { int rtnval;
   rtnval = nc_redef(ncidP->fileid); assert(chk_nc(rtnval));
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "num_target_cells", NC_INT, 1, &m_num_target_cells); assert(chk_nc(rtnval)); 
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "num_active_cells", NC_INT, 1, &m_num_active_cells); assert(chk_nc(rtnval)); 
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "num_failed_cells", NC_INT, 1, &m_num_failed_cells); assert(chk_nc(rtnval)); 
   rtnval = nc_put_att_int(ncidP->fileid, NC_GLOBAL, "num_vpr_tmp_issues", NC_INT, 1, &m_num_vpr_tmp_issues); assert(chk_nc(rtnval)); 
   rtnval = nc_enddef(ncidP->fileid); assert(chk_nc(rtnval));  
   rtnval = nc_close(ncidP->fileid); assert(chk_nc(rtnval));
   } // end of saveOutputFile()


void Simulation::verifyLatLon(ncid_type * ncidP, int row, int col, double cell_lat, double cell_lon)
   {
   size_t coords[1];
   int rtnval;
   double lat_double, lon_double;
   float lat_float, lon_float;
   
   coords[0] = row - ncidP->row_offset;
   switch (ncidP->row_vartype)
      {
      case NC_DOUBLE: rtnval = nc_get_var1_double(ncidP->fileid, ncidP->row_varid, coords, &lat_double); assert(chk_nc(rtnval)); break;
      case NC_FLOAT: 
         rtnval = nc_get_var1_float(ncidP->fileid, ncidP->row_varid, coords, &lat_float); assert(chk_nc(rtnval)); 
         lat_double = lat_float; 
         break;
      default: return; // assert(0);
      }
   
   coords[0] = col - ncidP->col_offset;
   switch (ncidP->col_vartype)
      {
      case NC_DOUBLE: rtnval = nc_get_var1_double(ncidP->fileid, ncidP->col_varid, coords, &lon_double); assert(chk_nc(rtnval)); break;
      case NC_FLOAT: 
         rtnval = nc_get_var1_float(ncidP->fileid, ncidP->col_varid, coords, &lon_float); assert(chk_nc(rtnval)); 
         lon_double = lon_float; 
         break;
      default: assert(0);
      }
   
   if (!(sciFn.close_enough(lat_double, cell_lat, 0.0003, runParams.cell_spacing/10.) && sciFn.close_enough(lon_double, cell_lon, 0.0003, runParams.cell_spacing/10.)))
      {
      printf("*** cell lat, lon %lf, %lf isn't close enough to %s data lat, lon %lf, %lf for row %d, col %d\n"
             "i.e. row %d, col %d of grid %s\n",
             cell_lat, cell_lon, ncidP->var_name, lat_double, lon_double, row - runParams.row_begin, col - runParams.col_begin, row, col, 
             runParams.grid_name);
      assert(0);
      }

   } // end of verifyLatLon()


void Simulation::getEarthData(InputDataClass * pInputData, SoilDataStruct * soil_dataP, int grid_row, int grid_col, bool ckLatLonFlag, double cell_lat, double cell_lon)
   { 
      size_t coords[3];
      size_t count[3];
      int rtnval;
      
      if (ckLatLonFlag) verifyLatLon(&elev_ncid, grid_row, grid_col, cell_lat, cell_lon);  
      coords[elev_ncid.rowNdxPos] = grid_row - elev_ncid.row_offset;
      coords[elev_ncid.colNdxPos] = grid_col - elev_ncid.col_offset;
      if (elev_ncid.timeNdxPos>=0) coords[elev_ncid.timeNdxPos] = 0;
      float * elevP = &(pInputData->elev);
      rtnval = nc_get_var1_float(elev_ncid.fileid, elev_ncid.input_data_varid, coords, elevP); 
      if (!chk_nc(rtnval))
         assert(false);
      *elevP *= elev_ncid.scale_factor;
      
      if (bd_ncid.fileid!=-1)
         {
         if (ckLatLonFlag) verifyLatLon(&bd_ncid, grid_row, grid_col, cell_lat, cell_lon);
         coords[bd_ncid.rowNdxPos] = grid_row - bd_ncid.row_offset;
         coords[bd_ncid.colNdxPos] = grid_col - bd_ncid.col_offset;
         if (bd_ncid.timeNdxPos>=0) coords[bd_ncid.timeNdxPos] = 0;
         switch (bd_ncid.input_data_vartype)
            { short shortval; int intval;
            case NC_FLOAT:
               rtnval = nc_get_var1_float(bd_ncid.fileid, bd_ncid.input_data_varid, coords, &(soil_dataP->bd)); assert(chk_nc(rtnval));
               break;
            case NC_SHORT:
               rtnval = nc_get_var1_short(bd_ncid.fileid, bd_ncid.input_data_varid, coords, &shortval); assert(chk_nc(rtnval));
               soil_dataP->bd = (float)shortval;
               break;
            case NC_INT:
               rtnval = nc_get_var1_int(bd_ncid.fileid, bd_ncid.input_data_varid, coords, &intval); assert(chk_nc(rtnval));
               soil_dataP->bd = (float)intval;
               break;
            default:
               assert(0);
               break;
            }
         soil_dataP->bd *= bd_ncid.scale_factor;
         }
      else soil_dataP->bd = NC_FILL_FLOAT;
      
      if (ckLatLonFlag) verifyLatLon(&soilData_ncid, grid_row, grid_col, cell_lat, cell_lon); 
      coords[soilData_ncid.rowNdxPos] = grid_row - soilData_ncid.row_offset; count[soilData_ncid.rowNdxPos] = 1;
      coords[soilData_ncid.colNdxPos] = grid_col - soilData_ncid.col_offset; count[soilData_ncid.colNdxPos] = 1;
      coords[soilData_ncid.timeNdxPos] = 0; count[soilData_ncid.timeNdxPos] = 10;
      switch (soilData_ncid.input_data_vartype)
         { short shortvals[10]; int intvals[10];
         case NC_FLOAT:
            rtnval = nc_get_vara_float(soilData_ncid.fileid, soilData_ncid.input_data_varid, coords, count, &(soil_dataP->soils_record[0])); assert(chk_nc(rtnval));
            break;
         case NC_SHORT:
            rtnval = nc_get_vara_short(soilData_ncid.fileid, soilData_ncid.input_data_varid, coords, count, shortvals); 
            assert(chk_nc(rtnval));
            for (int k = 0; k<10; k++) soil_dataP->soils_record[k] = (float)shortvals[k];
            break;
         case NC_INT:
            rtnval = nc_get_vara_int(soilData_ncid.fileid, soilData_ncid.input_data_varid, coords, count, intvals); assert(chk_nc(rtnval));
            for (int k = 0; k<10; k++) soil_dataP->soils_record[k] = (float)intvals[k];
            break;
         default:
            assert(0);
            break;
         }
      if (soilData_ncid.scale_factor!=1.0f) for (int k = 0; k<10; k++) soil_dataP->soils_record[k] *= soilData_ncid.scale_factor;
      
      if (runParams.baseCalibration==mc2WCR)
         { 
         float * fogP = &(pInputData->fog);
         if (fog_ncid.fileid!=-1)
            { // read the fog factor
            if (ckLatLonFlag) verifyLatLon(&fog_ncid, grid_row, grid_col, cell_lat, cell_lon);  
            coords[fog_ncid.rowNdxPos] = grid_row - fog_ncid.row_offset;
            coords[fog_ncid.colNdxPos] = grid_col - fog_ncid.col_offset;
            rtnval = nc_get_var1_float(fog_ncid.fileid, fog_ncid.input_data_varid, coords, fogP); assert(chk_nc(rtnval));
            *fogP *= fog_ncid.scale_factor;
            }
         else *fogP = NC_FILL_FLOAT;

         float * topomoistP = &(pInputData->topomoist);
         if (topomoist_ncid.fileid!=-1)
            { // read the topographic moisture factor
            if (ckLatLonFlag) verifyLatLon(&topomoist_ncid, grid_row, grid_col, cell_lat, cell_lon);  
            coords[topomoist_ncid.rowNdxPos] = grid_row - topomoist_ncid.row_offset;
            coords[topomoist_ncid.colNdxPos] = grid_col - topomoist_ncid.col_offset;
            rtnval = nc_get_var1_float(topomoist_ncid.fileid, topomoist_ncid.input_data_varid, coords, topomoistP); assert(chk_nc(rtnval));
            *topomoistP *= topomoist_ncid.scale_factor;
            }
         else *topomoistP = NC_FILL_FLOAT;
         
         float * deltaTslP = &(pInputData->deltaTsl);
         if (deltaTsl_ncid.fileid!=-1)
            { // read the deltaTsl layer 
            if (ckLatLonFlag) verifyLatLon(&deltaTsl_ncid, grid_row, grid_col, cell_lat, cell_lon);  
            coords[deltaTsl_ncid.rowNdxPos] = grid_row - deltaTsl_ncid.row_offset;
            coords[deltaTsl_ncid.colNdxPos] = grid_col - deltaTsl_ncid.col_offset;
            rtnval = nc_get_var1_float(deltaTsl_ncid.fileid, deltaTsl_ncid.input_data_varid, coords, deltaTslP); assert(chk_nc(rtnval));
            *deltaTslP *=deltaTsl_ncid.scale_factor;
            }
         else *deltaTslP = NC_FILL_FLOAT;
         
         float * aspectP = &(pInputData->aspect);
         if (aspect_ncid.fileid!=-1)
            { // read the aspect layer
            if (ckLatLonFlag) verifyLatLon(&aspect_ncid, grid_row, grid_col, cell_lat, cell_lon);  
            coords[aspect_ncid.rowNdxPos] = grid_row - aspect_ncid.row_offset;
            coords[aspect_ncid.colNdxPos] = grid_col - aspect_ncid.col_offset;
            rtnval = nc_get_var1_float(aspect_ncid.fileid, aspect_ncid.input_data_varid, coords, aspectP); assert(chk_nc(rtnval));
            *aspectP *= aspect_ncid.scale_factor;
            }
         else *aspectP = NC_FILL_FLOAT;
         
         float * swP = &(pInputData->sw);
         if (sw_ncid.fileid!=-1)
            { // read the shortwave layer
            if (ckLatLonFlag) verifyLatLon(&sw_ncid, grid_row, grid_col, cell_lat, cell_lon);  
            coords[sw_ncid.rowNdxPos] = grid_row - sw_ncid.row_offset;
            coords[sw_ncid.colNdxPos] = grid_col - sw_ncid.col_offset;
            rtnval = nc_get_var1_float(sw_ncid.fileid, sw_ncid.input_data_varid, coords, swP); assert(chk_nc(rtnval));
            *swP *= sw_ncid.scale_factor;
            }
         else *swP = NC_FILL_FLOAT;
         
         float * cadP = &(pInputData->cad);
         if (cad_ncid.fileid!=-1)
            { // read the cold air drainage factor
            if (ckLatLonFlag) verifyLatLon(&cad_ncid, grid_row, grid_col, cell_lat, cell_lon);  
            coords[cad_ncid.rowNdxPos] = grid_row - cad_ncid.row_offset;
            coords[cad_ncid.colNdxPos] = grid_col - cad_ncid.col_offset;
            rtnval = nc_get_var1_float(cad_ncid.fileid, cad_ncid.input_data_varid, coords, cadP); assert(chk_nc(rtnval));
            *cadP *= cad_ncid.scale_factor;
            }
         else *cadP = NC_FILL_FLOAT;
         
         } // end of if (runParams.baseCalibration==mc2WCR)
      
   } // end of getEarthData()


bool Simulation::getClimateVar(ncid_type * ncidP, float * varP, int grid_row, int grid_col, int num_months, bool ckLatLonFlag, double cell_lat, double cell_lon)
   {
   int months_offset = 12*runParams.years_offset_into_input_data;
   bool rtn_flag;
   
   rtn_flag = getClimateVar(ncidP, varP, grid_row, grid_col, num_months, ckLatLonFlag, cell_lat, cell_lon,
         months_offset);
   return(rtn_flag);
   
   } // end of getClimateVar(ncidP, varP, grid_row, grid_col, num_months, ckLatLonFlag, cell_lat, cell_lon)
   
   
bool Simulation::getClimateVar(ncid_type * ncidP, float * varP, int grid_row, int grid_col, int num_months, bool ckLatLonFlag, double cell_lat, double cell_lon, int months_offset)
   {
   size_t coords[3], count[3];
   int rtnval;
   short * shortvalP;
   int * intvalP;
   double * doublevalP;
   
   assert(ncidP->fileid!=-1 && ncidP->row_offset!=-9999 && ncidP->col_offset!=-9999);
   if (ckLatLonFlag) verifyLatLon(ncidP, grid_row, grid_col, cell_lat, cell_lon);  
   coords[ncidP->rowNdxPos] = grid_row - ncidP->row_offset; count[ncidP->rowNdxPos] = 1;
   coords[ncidP->colNdxPos] = grid_col - ncidP->col_offset; count[ncidP->colNdxPos] = 1;

   // coords[ncidP->timeNdxPos] = 0; count[ncidP->timeNdxPos] = num_months;
   coords[ncidP->timeNdxPos] = months_offset; count[ncidP->timeNdxPos] = num_months;
   
   // logic to deal with shorts, ints, etc. here
   switch (ncidP->input_data_vartype)
      {
      case NC_FLOAT:
         rtnval = nc_get_vara_float(ncidP->fileid, ncidP->input_data_varid, coords, count, varP); assert(chk_nc(rtnval));
         break;
      case NC_DOUBLE:
         doublevalP = (double *)malloc(num_months*sizeof(double)); assert(doublevalP!=NULL);
         rtnval = nc_get_vara_double(ncidP->fileid, ncidP->input_data_varid, coords, count, doublevalP); assert(chk_nc(rtnval));
         for (int k = 0; k<num_months; k++) *(varP + k) = (float)(*(doublevalP + k));         
         break;
      case NC_SHORT:
         shortvalP = (short *)malloc(num_months*sizeof(short)); assert(shortvalP!=NULL);
         rtnval = nc_get_vara_short(ncidP->fileid, ncidP->input_data_varid, coords, count, shortvalP); 
         if (!chk_nc(rtnval))
               assert(false);
         for (int k = 0; k<num_months; k++) *(varP + k) = *(shortvalP + k);
         free(shortvalP);
         break;
      case NC_INT:
         intvalP = (int *)malloc(num_months*sizeof(int)); assert(intvalP!=NULL);
         rtnval = nc_get_vara_int(ncidP->fileid, ncidP->input_data_varid, coords, count, intvalP); assert(chk_nc(rtnval));
         for (int k = 0; k<num_months; k++) *(varP + k) = *(intvalP + k);
         free(intvalP);
         break;
      default:
         assert(0);
         break;
      }
   
   // check for missing data or fill values, but only for the first value in the time series
   if (ncidP->fill_value_flag && ncidP->fill_value==(*varP)) return(false);
   if (ncidP->missing_value_flag && ncidP->missing_value==(*varP)) return(false);
   
   // Scale the data appropriately.
   if (ncidP->scale_factor!=1.0) for (int k = 0; k<num_months; k++) *(varP + k) *= ncidP->scale_factor;
   
   return(true);
   
   } // end of getClimateVar()


/*
bool Simulation::maskedOut(int row_index, int col_index)
  { 
  int grid_row, grid_col;
  size_t coords[3];
  short shortVal;
  int rtnval;
  int maskval;
  bool maskedOutFlag;
  signed char scharVal;
  float floatVal;

  if (mask_ncid.fileid==-1) return(false);
  
  grid_row = row_index + runParams.row_begin;
  grid_col = col_index + runParams.col_begin;
  coords[mask_ncid.rowNdxPos] = grid_row - mask_ncid.row_offset;
  coords[mask_ncid.colNdxPos] = grid_col - mask_ncid.col_offset;
  if (mask_ncid.timeNdxPos>=0) coords[mask_ncid.timeNdxPos] = 0;
  
  switch (mask_ncid.input_data_vartype)
    {
    case NC_INT:
      rtnval = nc_get_var1_int(mask_ncid.fileid, mask_ncid.input_data_varid, coords, &maskval); 
      break;
    case NC_SHORT:
      rtnval = nc_get_var1_short(mask_ncid.fileid, mask_ncid.input_data_varid, coords, &shortVal); 
      maskval = shortVal;
      break;
    case NC_BYTE:
      rtnval = nc_get_var1_schar(mask_ncid.fileid, mask_ncid.input_data_varid, coords, &scharVal);
      maskval = scharVal;
      break;
    case NC_FLOAT:
      rtnval = nc_get_var1_float(mask_ncid.fileid, mask_ncid.input_data_varid, coords, &floatVal);
      maskval = (int)floatVal;
      break;
    
    default:
      assert(false);
      break;
    }
  
  maskedOutFlag = rtnval!=NC_NOERR || maskval<=0;
  return(maskedOutFlag); 
  } // end of Simulation::maskedOut()
*/
  
  
bool Simulation::getInputData(int month, int grid_row, int grid_col, bool ckLatLonFlag, InputDataClass * pInputData)
   {
   size_t coords[3];
   int rtnval;
   bool rtn_flag = true;
   
   pInputData->lat = pInputData->northing = north_south_centroid(grid_row);
   pInputData->lon = pInputData->easting = east_west_centroid(grid_col);
   
   getEarthData(pInputData, &pInputData->soilData, grid_row, grid_col, ckLatLonFlag, pInputData->lat, pInputData->lon);

   switch (runParams.runMode)
      {
      case MAPSS_EQ:
         // MAPSS doesn't use bulk density.
         // If there is no bulk density input data file, put in a nominal bulk density value so that 
         // the Soil_MAPSS constructor doesn't flag the soil data as invalid.
         if (bd_ncid.fileid==-1) pInputData->soilData.bd = 1.5; 
         pInputData->num_months = 12.;
         break;
         
      case CENTURY_EQ:
         pInputData->num_months = 12;
         
         // read the file created by the MAPSS phase
         assert(MAPSSdata_ncid.fileid!=-1 && MAPSSdata_ncid.rowNdxPos>=0 && MAPSSdata_ncid.colNdxPos>=0);
         coords[MAPSSdata_ncid.rowNdxPos] = grid_row - MAPSSdata_ncid.row_offset;
         coords[MAPSSdata_ncid.colNdxPos] = grid_col - MAPSSdata_ncid.col_offset;
         if (MAPSSdata_ncid.timeNdxPos>=0) coords[MAPSSdata_ncid.timeNdxPos] = 0;
         
         assert(pInputData->inFromMAPSS==NULL);      
         pInputData->inFromMAPSS = (float *)malloc(NUM_MAPSS_OUTVARS*sizeof(float)); 
         // Gets free'd in the destructor for the InputDataClass object
         assert(pInputData->inFromMAPSS!=NULL);
         
         for (int i = 0; i<NUM_MAPSS_OUTVARS; i++) 
            { float floatval; short shortval;
            switch (MAPSSoutputList[i].type)
               {
               case NC_SHORT:
                  rtnval = nc_get_var1_short(MAPSSdata_ncid.fileid, MAPSSoutputList[i].var_id, coords, &shortval);
                  floatval = (float)shortval;
                  break;
               case NC_FLOAT:
                  rtnval = nc_get_var1_float(MAPSSdata_ncid.fileid, MAPSSoutputList[i].var_id, coords, &floatval);
                  break;
               default:
                  assert(0);
                  break;
               } // end of switch on type        
            if (chk_nc(rtnval)) *(pInputData->inFromMAPSS + i) = floatval;
            else 
               {
               if (rtnval==NC_EINVALCOORDS) printf("Data for row %d, col %d is missing from the .mapss file.\n", grid_row, grid_col);
               return(false);
               }
            }
         
         break; // end of case CENTURY_EQ
         
      case SPINUP: // read the file containing the MC2centuryEQdataOutputs object created by the CENTURY_EQ phase
         assert(ppt_ncid.fileid!=-1);
         pInputData->num_months = ppt_ncid.time_dimlen;
         assert(EQdata_ncid.fileid!=-1);
         coords[EQdata_ncid.rowNdxPos] = grid_row - EQdata_ncid.row_offset;
         coords[EQdata_ncid.colNdxPos] = grid_col - EQdata_ncid.col_offset;
         if (EQdata_ncid.timeNdxPos>=0) coords[EQdata_ncid.timeNdxPos] = 0;
         
         assert(pInputData->inFromEQ==NULL);
         pInputData->inFromEQ = (float *)malloc(NUM_EQ_OUTVARS*sizeof(float)); 
         // Gets free'd in the destructor for the InputDataClass object
         assert(pInputData->inFromEQ!=NULL);
         for (int i = 0; i<NUM_EQ_OUTVARS; i++) if (EQandWSinputList[i].show)
            { float floatval;
            rtnval = nc_get_var1_float(EQdata_ncid.fileid, EQandWSinputList[i].var_id, coords, &floatval);
            if (chk_nc(rtnval)) *((pInputData->inFromEQ) + i) = floatval;
            else 
                 {
                 printf("Data for var %d, row %d, col %d is missing from the .eq file.\n", i, grid_row, grid_col);
                 printf("EQandWSinputList[i]. var_id, name, show = %d, %s, %d\n", 
                     EQandWSinputList[i].var_id, EQandWSinputList[i].name, EQandWSinputList[i].show);
                 printf("rtnval = %d\n", rtnval);
                 return(false);
                 }
            }
         break; // end of case SPINUP
                  
      case TRANSIENT:
         pInputData->num_months = runParams.actual_years_to_run*12;
         assert(ppt_ncid.fileid!=-1 
               && ppt_ncid.time_dimlen>=
               ((unsigned int)pInputData->num_months + 12*runParams.years_offset_into_input_data));
         assert(WSdata_ncid.fileid!=-1);
         coords[WSdata_ncid.rowNdxPos] = grid_row - WSdata_ncid.row_offset;
         coords[WSdata_ncid.colNdxPos] = grid_col - WSdata_ncid.col_offset;
         if (WSdata_ncid.timeNdxPos>=0) coords[WSdata_ncid.timeNdxPos] = 0;
         
         assert(pInputData->inFromWS==NULL);
         pInputData->inFromWS = (float *)malloc(NUM_WS_OUTVARS*sizeof(float)); 
         // Gets free'd in the destructor for the InputDataClass object
         assert(pInputData->inFromWS!=NULL);
         
         if (modelParams.code_flags[FAST_WARMSTART_IN_FLAG])
            { 
            assert(WSdata_ncid.ws_varid!=-1);
            size_t ws_coords[4], ws_counts[4];
            ws_coords[WSdata_ncid.timeNdxPos] = 0;
            ws_coords[WSdata_ncid.rowNdxPos] = coords[WSdata_ncid.rowNdxPos];
            ws_coords[WSdata_ncid.colNdxPos] = coords[WSdata_ncid.colNdxPos];
            ws_counts[WSdata_ncid.timeNdxPos] = ws_counts[WSdata_ncid.rowNdxPos] = ws_counts[WSdata_ncid.colNdxPos] = 1;
            ws_coords[WSdata_ncid.wsNdxPos] = 0;
            ws_counts[WSdata_ncid.wsNdxPos] = NUM_WS_OUTVARS;
            rtnval = nc_get_vara_float(WSdata_ncid.fileid, WSdata_ncid.ws_varid, ws_coords, ws_counts, pInputData->inFromWS);
            if (!chk_nc(rtnval))
              {
              printf("*** getInputData(): ws_varid, ws_coords[], ws_counts[] = %d\n %ld, %ld, %ld, %ld\n %ld, %ld, %ld, %ld\n",
                  WSdata_ncid.ws_varid,
                  ws_coords[0], ws_coords[1], ws_coords[2], ws_coords[3], 
                  ws_counts[0], ws_counts[1], ws_counts[2], ws_counts[3]);
              assert(chk_nc(rtnval));
              }
            rtn_flag = true;
            }
         else 
            {
            for (int i = 0; i<NUM_WS_OUTVARS; i++) if (EQandWSinputList[i].show)
              { float floatval;
              rtnval = nc_get_var1_float(WSdata_ncid.fileid, EQandWSinputList[i].var_id, coords, &floatval);
              if (chk_nc(rtnval)) *(pInputData->inFromWS + i) = floatval;
              else 
                 {
                 printf("Data for var %d, row %d, col %d is missing from the .ws file.\n", i, grid_row, grid_col);
                 printf("EQandWSinputList[i]. var_id, name, show = %d, %s, %d\n", 
                     EQandWSinputList[i].var_id, EQandWSinputList[i].name, EQandWSinputList[i].show);
                 printf("rtnval = %d\n", rtnval);
                 return(false);
                 }
              }
            }
         break; // end of case TRANSIENT
         
      default:
         assert(0);
         break;
      }
   pInputData->pptP = (float *)malloc(pInputData->num_months*sizeof(float)); assert(pInputData->pptP!=NULL);
   rtn_flag &= getClimateVar(&ppt_ncid, pInputData->pptP, grid_row, grid_col, pInputData->num_months, ckLatLonFlag, pInputData->lat, pInputData->lon);  
   
   if (tmin_ncid.fileid!=-1)
      {
      pInputData->tminP = (float *)malloc(pInputData->num_months*sizeof(float)); assert(pInputData->tminP!=NULL);
      rtn_flag &= getClimateVar(&tmin_ncid, pInputData->tminP, grid_row, grid_col, pInputData->num_months, ckLatLonFlag, pInputData->lat, pInputData->lon);  
      }
   else pInputData->tminP = NULL;
   
   if (tmax_ncid.fileid!=-1)
      {
      pInputData->tmaxP = (float *)malloc(pInputData->num_months*sizeof(float)); assert(pInputData->tmaxP!=NULL);
      rtn_flag &= getClimateVar(&tmax_ncid, pInputData->tmaxP, grid_row, grid_col, pInputData->num_months, ckLatLonFlag, pInputData->lat, pInputData->lon);  
      }
   else pInputData->tmaxP = NULL;
   
   pInputData->tmpP = (float *)malloc(pInputData->num_months*sizeof(float)); assert(pInputData->tmpP!=NULL);
   if (tmp_ncid.fileid!=-1) rtn_flag &= getClimateVar(&tmp_ncid, pInputData->tmpP, grid_row, grid_col, pInputData->num_months, ckLatLonFlag, pInputData->lat, pInputData->lon);  
   else 
      { float tmin, tmax, tmp;
      assert(tmin_ncid.fileid!=-1 && tmax_ncid.fileid!=-1);
      for (int k = 0; k<pInputData->num_months; k++) // *(pInputData->tmpP + k) = (*(pInputData->tminP + k) + *(pInputData->tmaxP + k))/2.f;
         { 
         tmin = *(pInputData->tminP + k);
         tmax = *(pInputData->tmaxP + k);
         tmp = (tmin + tmax)/2.f;
         *(pInputData->tmpP + k) = tmp;
         }
      }
   
   if (tdmean_ncid.fileid!=-1) 
      {
      pInputData->tdmeanP = (float *)malloc(pInputData->num_months*sizeof(float)); assert(pInputData->tdmeanP!=NULL);
      rtn_flag &= getClimateVar(&tdmean_ncid, pInputData->tdmeanP, grid_row, grid_col, pInputData->num_months, ckLatLonFlag, pInputData->lat, pInputData->lon);  
      }
   else pInputData->tdmeanP = NULL;
   
   pInputData->vprP = (float *)malloc(pInputData->num_months*sizeof(float)); assert(pInputData->vprP!=NULL);
   if (vpr_ncid.fileid!=-1) rtn_flag &= getClimateVar(&vpr_ncid, pInputData->vprP, grid_row, grid_col, pInputData->num_months, ckLatLonFlag, pInputData->lat, pInputData->lon); 
   else 
      {
      assert(tdmean_ncid.fileid!=-1);
      for (int k = 0; k<pInputData->num_months; k++)
            { 
            float tdmeanVal = *(pInputData->tdmeanP + k);
            assert(tdmeanVal>-FREEZE);
            *(pInputData->vprP + k) = sciFn.satvp(tdmeanVal);  
            }  
      }
   
   if (wnd_ncid.fileid!=-1) 
      {
      pInputData->wndP = (float *)malloc(12*sizeof(float)); assert(pInputData->wndP!=NULL);
      rtn_flag &= getClimateVar(&wnd_ncid, pInputData->wndP, grid_row, grid_col, 12, ckLatLonFlag, pInputData->lat, pInputData->lon, 0);
      }
   else pInputData->wndP = NULL;
   
   return(rtn_flag);
  
   } // end of Simulation::getInputData()


bool Simulation::interpretSwitch(char * switchStr)
   {
   bool rtnFlag = false;
   
   if (switchStr==NULL) rtnFlag = FALSE;
   else if (strcmp(switchStr, "ON")==0) rtnFlag = true;
   else if (strcmp(switchStr, "OFF")!=0) err_exit("switch value was not recognized as ON or OFF");
   
   return(rtnFlag);
   } // end of Simulation::interpretSwitch()


RunModeEnum RunParamsClass::interpretRunMode(char * run_mode)
{
  RunModeEnum rm;
  
  if (run_mode==NULL) rm = unknownRunMode;
  else if (strcmp(run_mode, "MAPSS_EQ")==0) rm = MAPSS_EQ;
  else if (strcmp(run_mode, "CENTURY_EQ")==0) rm = CENTURY_EQ;
  else if (strcmp(run_mode, "SPINUP")==0) rm = SPINUP;
  else if (strcmp(run_mode, "TRANSIENT")==0) rm = TRANSIENT;
  else rm = unknownRunMode;
  
  return(rm);
} // end of interpretRunMode()


GridNameEnum RunParamsClass::interpretGridName(char * grid_name)
{
  GridNameEnum gn;
  
  if (grid_name==NULL) gn = unknownGridName;
  else if (strcmp(grid_name, "Global")==0) gn = Global;
  else if (strcmp(grid_name, "VEMAP")==0) gn = VEMAP;
  else if (strcmp(grid_name, "US12km")==0) gn = US12km;
  else if (strcmp(grid_name, "NA8km")==0) gn = NA8km;
  else if (strcmp(grid_name, "US4km")==0) gn = US4km;
  else if (strcmp(grid_name, "USeast4km")==0) gn = USeast4km;
  else if (strcmp(grid_name, "US800m")==0) gn = US800m;
  else if (strcmp(grid_name, "USwest800m")==0) gn = USwest800m;
  else if (strcmp(grid_name, "BearSoil800m")==0) gn = BearSoil800m;
  else if (strcmp(grid_name, "ATtest800m")==0) gn = ATtest800m;
  else if (strcmp(grid_name, "CA800m")==0) gn = CA800m;
  else if (strcmp(grid_name, "YNP800m")==0) gn = YNP800m;
  else if (strcmp(grid_name, "PNW800m")==0) gn = PNW800m;
  else if (strcmp(grid_name, "Yellowstone800m")==0) gn = Yellowstone800m;
  else if (strcmp(grid_name, "CA12km")==0) gn = CA12km;
  else if (strcmp(grid_name, "CA10kmAlbers")==0) gn = CA10kmAlbers;
  else if (strcmp(grid_name, "US10kmAlbers")==0) gn = US10kmAlbers;
  else if (strcmp(grid_name, "BLM_PSW4kmAlbers")==0) gn = BLM_PSW4kmAlbers;
  else gn = unknownGridName;
  
  int nrows, ncols;
  getGrid(grid_name, &latitude0, &longitude0, &cell_spacing, &nrows, &ncols);
  
  return(gn);
} // end of interpretGridName()


bool Simulation::getModelParameters(BaseCalibrationEnum bc)
   { 
      // Initialize fire thresholds by climate zone and tree type
      // typedef enum {UNKNOWNzone = 0, ARCTICzone, BOREALzone, TEMPERATEzone, SUBTROPICALzone, TROPICALzone} ClimateZone;
      // typedef enum {UNKNOWNtree_typ=0, EN_TREES, EN_DB_TREES, DB_TREES, DB_EB_TREES, EN_EB_TREES, EB_TREES, 
      // DN_TREES, DN_EN_TREES} TreeType;
      const float ffmcThreshold[6][9] = { // eight tree types for each of 5 climate zones
      { -9999., -9999., -9999., -9999., -9999., -9999., -9999., -9999., -9999.}, // zone 0 UNKNOWNzone
      //           EN    EN_DB     DB    DB_EB   EN_EB     EB      DN    DN_EN
      { -9999.,    94.,    94.,    94.,    94.,    94.,    94.,    94.,    94.}, // 1 ARCTICzone
      { -9999.,  83.95,    92.,    92.,    87.,  91.42,  91.42,  83.95,  83.95}, // 2 BOREALzone,
      { -9999., 89.13744,  92.,    92.,    87.,  91.42,  91.42, 89.13744, 89.13744}, // 3 TEMPERATEzone
      { -9999., 89.13744,  92.,    92.,    87.,  91.42,  91.42, 89.13744, 89.13744}, // 4 SUBTROPICALzone
      { -9999., 89.13744,  92.,    92.,    87.,  91.42,  91.42, 89.13744, 89.13744}}; // 5 TROPICALzone 
      const float buiThreshold[6][9] = {  // eight tree types for each of 5 climate zones
      { -9999., -9999., -9999., -9999., -9999., -9999., -9999., -9999., -9999.}, // zone 0 UNKNOWNzone
      { -9999.,   147.,   147.,   147.,   147.,   147.,   147.,   147.,   147.}, // 1 ARCTICzone
      { -9999.,  38.20,   150.,   150.,   110., 223.11, 223.11,  38.20,  38.20}, // 2 BOREALzone,
      { -9999.,   245.,   150.,   150.,   110., 223.11, 223.11,   245.,   245.}, // 3 TEMPERATEzone
      { -9999.,   245.,   150.,   150.,   110., 223.11, 223.11,   245.,   245.}, // 4 SUBTROPICALzone
      { -9999.,   245.,   150.,   150.,   110., 223.11, 223.11,   245.,   245.}}; // 5 TROPICALzone 
      for (int zone = 0; zone<6; zone++) for (int treeType = 0; treeType<9; treeType++) 
         {
         modelParams.ffmc_threshold[zone][treeType] = ffmcThreshold[zone][treeType];
         modelParams.bui_threshold[zone][treeType] = buiThreshold[zone][treeType];
         }
/*
      TreeType treeTypeFromVtype[MAX_VTYPE+1];
      treeTypeFromVtype[UNKNOWNveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[COLD_BARRENveg] = UNKNOWNtree_typ; 
      treeTypeFromVtype[TUNDRAveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[TAIGA_TUNDRAveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[BOREAL_NEEDLELEAF_FORESTveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[BOREAL_WOODLANDveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[SUBALPINE_FORESTveg] = EN_TREES;
      treeTypeFromVtype[MARITIME_EN_FORESTveg] = EN_TREES;
      treeTypeFromVtype[MESIC_TEMPERATE_NEEDLELEAF_FORESTveg] = EN_TREES;
      treeTypeFromVtype[TEMPERATE_DB_FORESTveg] = DB_TREES;
      treeTypeFromVtype[COOL_MIXED_FORESTveg] = EN_DB_TREES;
      treeTypeFromVtype[TEMPERATE_WARM_MIXED_FORESTveg] = EN_DB_TREES;
      treeTypeFromVtype[TEMPERATE_EN_WOODLANDveg] = EN_TREES;
      treeTypeFromVtype[TEMPERATE_DB_WOODLANDveg] = DB_TREES;
      treeTypeFromVtype[TEMPERATE_COOL_MIXED_WOODLANDveg] = EN_DB_TREES;
      treeTypeFromVtype[TEMPERATE_WARM_MIXED_WOODLANDveg] = EN_DB_TREES;
      treeTypeFromVtype[C3SHRUBveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[C3GRASSveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[TEMPERATE_DESERTveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[SUBTROPICAL_EN_FORESTveg] = EN_TREES;
      treeTypeFromVtype[SUBTROPICAL_DB_FORESTveg] = DB_TREES;
      treeTypeFromVtype[WARM_EB_FORESTveg] = EB_TREES;
      treeTypeFromVtype[SUBTROPICAL_MIXED_FORESTveg] = EN_DB_TREES;
      treeTypeFromVtype[SUBTROPICAL_EN_WOODLANDveg] = EN_TREES; 
      treeTypeFromVtype[SUBTROPICAL_DB_WOODLANDveg] = DB_TREES;
      treeTypeFromVtype[SUBTROPICAL_EB_WOODLANDveg] = EB_TREES;
      treeTypeFromVtype[SUBTROPICAL_MIXED_WOODLANDveg] = EN_DB_TREES;
      treeTypeFromVtype[C4SHRUBveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[C4GRASSveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[SUBTROPICAL_DESERTveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[TROPICAL_EB_FORESTveg] = EB_TREES;
      treeTypeFromVtype[TROPICAL_DECIDUOUS_WOODLANDveg] = DB_TREES;
      treeTypeFromVtype[TROPICAL_SAVANNAveg] = EN_EB_TREES;
      treeTypeFromVtype[VTYPE33veg] = UNKNOWNtree_typ;
      treeTypeFromVtype[VTYPE34veg] = UNKNOWNtree_typ;
      treeTypeFromVtype[TROPICAL_DESERTveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[MOIST_TEMPERATE_NEEDLELEAF_FORESTveg] = EN_TREES;
      treeTypeFromVtype[VTYPE37veg] = UNKNOWNtree_typ;
      treeTypeFromVtype[SUBALPINE_MEADOWveg] = EN_TREES;
      treeTypeFromVtype[WATERveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[NATURAL_BARRENveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[DEVELOPEDveg] = UNKNOWNtree_typ;
      treeTypeFromVtype[LARCH_FORESTveg] = DN_TREES;
      treeTypeFromVtype[SSZveg] = EN_TREES;
      treeTypeFromVtype[WHZveg] = EN_TREES;
      treeTypeFromVtype[PSFZveg] = EN_TREES;
      treeTypeFromVtype[MHZveg] = EN_TREES;
      treeTypeFromVtype[SAFZveg] = EN_TREES;
      treeTypeFromVtype[PKLZveg] = EN_TREES;
      treeTypeFromVtype[DRY_TEMPERATE_NEEDLELEAF_FORESTveg] = EN_TREES;
      // treeTypeFromVtype[XERIC_NEEDLELEAF_WOODLANDveg] = EN_TREES; 
*/
      for (int vt = 0; vt<=MAX_VTYPE; vt++)
         {
         modelParams.ffmc_threshold_by_vtype[vt] = -9999.f; 
         modelParams.bui_threshold_by_vtype[vt] = -9999.f; 
         }
          
      
      // Initialize fire return parameters.
      const float vveg2mfriGlobal[][2] = {
      {1000., 1000.}, // 0 unused
      {9000., 9000.}, // 1 ice aka barren COLD_BARRENveg
      {1000., 1000.}, // 2 tundra aka alpine TUNDRAveg
      {1000., 1000.}, // 3 taiga-tundra TAIGA_TUNDRAveg
      {50., 50.}, // 4 boreal needleleaf forest BOREAL_NEEDLELEAF_FORESTveg
      {50., 50.}, // 5 boreal mixed woodland BOREAL_WOODLANDveg
      {300., 300.}, // 6 subalpine forest SUBALPINE_FORESTveg
      {150., 150.}, // 7 maritime needleleaf forest MARITIME_EN_FORESTveg
      {50., 50.}, // 8 temperate needleleaf forest MESIC_TEMPERATE_NEEDLELEAF_FORESTveg
      {155., 237.}, // 9 temperate deciduous broadleaf forest TEMPERATE_DB_FORESTveg
      {100., 100.}, // 10 temperate cool mixed forest COOL_MIXED_FORESTveg
      {75., 75.}, // 11 temperate warm mixed forest TEMPERATE_WARM_MIXED_FORESTveg
      {15., 15.}, // 12 temperate needleleaf woodland TEMPERATE_EN_WOODLANDveg
      {15., 15.}, // 13 temperate deciduous broadleaf woodland TEMPERATE_DB_WOODLANDveg
      {15., 15.}, // 14 temperate cool mixed woodland TEMPERATE_COOL_MIXED_WOODLANDveg
      {39., 75.}, // 15 temperate warm mixed woodland TEMPERATE_WARM_MIXED_WOODLANDveg
      {30., 30.}, // 16 C3 shrubland C3SHRUBveg
      {15., 32.}, // 17 C3 grassland C3GRASSveg
      {200., 200.}, // 18 temperate desert TEMPERATE_DESERTveg
      {50., 50.}, // 19 subtropical needleleaf forest SUBTROPICAL_EN_FORESTveg
      {100., 100.}, // 20 subtropical deciduous broadleaf forest SUBTROPICAL_DB_FORESTveg
      {40., 40.}, // 21 subtropical evergreen broadleaf forest WARM_EB_FORESTveg
      {40., 40.}, // 22 subtropical mixed forest SUBTROPICAL_MIXED_FORESTveg
      {15., 15.}, // 23 subtropical needleleaf woodland SUBTROPICAL_EN_WOODLANDveg
      {15., 15.}, // 24 subtropical deciduous broadleaf woodland SUBTROPICAL_DB_WOODLANDveg
      {15., 15.}, // 25 subtropical evergreen broadleaf woodland SUBTROPICAL_EB_WOODLANDveg
      {39., 75.}, // 26 subtropical mixed woodland SUBTROPICAL_MIXED_WOODLANDveg
      {30., 30.}, // 27 C4 shrubland C4SHRUBveg
      {12., 24.}, // 28 C4 grassland C4GRASSveg
      {200., 200.}, // 29 subtropical desert SUBTROPICAL_DESERTveg
      {150., 150.}, // 30 tropical evergreen broadleaf forest TROPICAL_EB_FORESTveg
      {15., 15.}, // 31 tropical deciduous woodland TROPICAL_DECIDUOUS_WOODLANDveg
      {15., 15.}, // 32 tropical savanna TROPICAL_SAVANNAveg
      {0., 0.,}, // 33 unused VTYPE33veg
      {0., 0.}, // 34 unused VTYPE34veg
      {200., 200.}, // 35 tropical desert TROPICAL_DESERTveg
      {50., 150.}, // 36 cool moist needleleaf forest MOIST_TEMPERATE_NEEDLELEAF_FORESTveg
      {0., 0.}, // 37 unused VTYPE37veg
      {30., 30.}, // 38 subalpine meadow SUBALPINE_MEADOWveg
      {1000., 1000.}, // 39 water and wetlands WATERveg
      {1000., 1000.}, // 40 natural barren NATURAL_BARRENveg
      {1000., 1000.}, // 41 developed DEVELOPEDveg
      {50., 50.}, // 42 larch forest 
      // Data for SSZ, WHZ, PSFZ, MHZ, and SAFZ is taken from Henderson et al., 
      // Forested Plant Associations of the Olympic National Forest, R6-ECOL-TP 001-88, USFS 1989, p. 19
      // and from an email from Josh Halofsky to Dave Conklin on 11/27/12
      {900., 1000.}, // 43 SSZ
      {234., 234.}, // 44 WHZ
      {625., 629.}, // 45 PSFZ
      {844., 844.}, // 46 MHZ
      {208., 222.}, // 47 SAFZ
      {50., 50.}, // 48 PKLZ, from 5 boreal mixed woodland BOREAL_WOODLANDveg
      {80., 190.}, // 49 DRY_TEMPERATE_NEEDLELEAF_FORESTveg from Rapid Assessment Reference Condition Model R#SPFI 8/11/08
      {130., 130.}, // 50 XERIC_NEEDLELEAF_WOODLANDveg from Rapid Assessment Reference Condition Model R#PIPOxe 8/11/08
      }; // end of vveg2mfriGlobal[][2]
      int num_records = sizeof(vveg2mfriGlobal)/(2*sizeof(float));
      assert(num_records>=(MAX_VTYPE+1));
      for (int vtypeNdx = 0; vtypeNdx<num_records; vtypeNdx++)
         {
         modelParams.vveg2mfri[vtypeNdx][0] = vtypeNdx;
         modelParams.vveg2mfri[vtypeNdx][1] = vveg2mfriGlobal[vtypeNdx][0];
         modelParams.vveg2mfri[vtypeNdx][2] = vveg2mfriGlobal[vtypeNdx][1];
         }
         
      modelParams.MAPSS_parameter_set = MpSN[defaultParams]; 
         
      // Set values of any parameters which are different from those in ModelParameters_B41 
      switch (bc)
      {
      case unknownBaseCalibration: assert(0); break;
/*
      case mc2VEMAP: 
         modelParams.MAPSS_parameter_set = MpSN[ORNLparams]; 
         break;
      case mc2NA8km:
         modelParams.MAPSS_parameter_set = MpSN[US10kmSCSparams]; 
         modelParams.subalpine_threshold = 1600.f;
         modelParams.bui_EN_threshold = 122.85f;
         modelParams.bui_EN_DB_threshold = 42.39f;
         modelParams.ffmc_EN_DB_threshold = 84.743f;
         err_exit("some parameter values have not been entered for NA8km (aka LYNX)"); 
         break;
      case mc2CA08:
         modelParams.MAPSS_parameter_set = MpSN[defaultParams]; 
         modelParams.maritime_threshold = 14.f;
         modelParams.subalpine_threshold = 2000.f;
         modelParams.bui_EN_threshold = 290.f;
         modelParams.bui_EN_DB_threshold = 310.f;
         modelParams.ffmc_EN_threshold = 87.f;
         modelParams.ffmc_EN_DB_threshold = 87.f;
         break;
      case mc2YOSE:
         modelParams.MAPSS_parameter_set = MpSN[defaultParams]; 
         modelParams.maritime_threshold = 14.f;
         modelParams.bui_EN_DB_threshold = 265.f;
         modelParams.ffmc_EN_threshold = 90.25f;
         modelParams.ffmc_EN_DB_threshold = 84.74f;
         err_exit("some parameter values have not been entered for YOSE");
         break;
      case mc2WIES_YOSE: 
         modelParams.MAPSS_parameter_set = MpSN[WiesYoseParams];     
         err_exit("MC2 parameter values have not been entered for WIES_YOSE"); 
         break;
      case mc2VINCERA:
         modelParams.MAPSS_parameter_set = MpSN[US10kmSCSparams]; 
         runParams.fire_suppression_switch = (char *)ONstr;
         runParams.fire_suppression_first_year = 1951;
         modelParams.maritime_threshold = 16.f;
         modelParams.fire_min = 0.0f;
         err_exit("some parameter values have not been entered for VINCERA"); 
         break;
      case mc2US10kmAlbers: 
         modelParams.MAPSS_parameter_set = MpSN[US10kmFAOparams];
         runParams.soil_data_file = (char *)FAOsoilsStr;
         err_exit("some parameter values have not been entered for US10kmAlbers"); 
         break;
*/
      case mc2GLOBAL: 
         break;
         
      case mc2WCR:
      case mc2ConUS:
         modelParams.maritime_threshold = 16.;
         modelParams.tmmin_threshold = 1.5;  
                
         // modelParams.bui_EN_threshold = 122.85;
         modelParams.bui_threshold[TEMPERATEzone][EN_TREES] = 122.85;         
         modelParams.bui_threshold[TEMPERATEzone][DN_TREES] = 122.85;         
         modelParams.bui_threshold[TEMPERATEzone][DN_EN_TREES] = 122.85;         
         modelParams.bui_threshold[SUBTROPICALzone][EN_TREES] = 122.85;         
         modelParams.bui_threshold[SUBTROPICALzone][DN_TREES] = 122.85;         
         modelParams.bui_threshold[SUBTROPICALzone][DN_EN_TREES] = 122.85;         
         modelParams.bui_threshold[TROPICALzone][EN_TREES] = 122.85;         
         modelParams.bui_threshold[TROPICALzone][DN_TREES] = 122.85;         
         modelParams.bui_threshold[TROPICALzone][DN_EN_TREES] = 122.85;         
         // modelParams.ffmc_EN_DB_threshold = 87.;
         modelParams.ffmc_threshold[BOREALzone][EN_DB_TREES] = 87.;
         modelParams.ffmc_threshold[BOREALzone][DB_TREES] = 87.;
         modelParams.ffmc_threshold[TEMPERATEzone][EN_DB_TREES] = 87.;
         modelParams.ffmc_threshold[TEMPERATEzone][DB_TREES] = 87.;
         modelParams.ffmc_threshold[SUBTROPICALzone][EN_DB_TREES] = 87.;
         modelParams.ffmc_threshold[SUBTROPICALzone][DB_TREES] = 87.;
         modelParams.ffmc_threshold[TROPICALzone][EN_DB_TREES] = 87.;
         modelParams.ffmc_threshold[TROPICALzone][DB_TREES] = 87.;
         // modelParams.bui_EN_DB_threshold = 110.;
         modelParams.bui_threshold[BOREALzone][EN_DB_TREES] = 110.;
         modelParams.bui_threshold[BOREALzone][DB_TREES] = 110.;
         modelParams.bui_threshold[TEMPERATEzone][EN_DB_TREES] = 110.;
         modelParams.bui_threshold[TEMPERATEzone][DB_TREES] = 110.;
         modelParams.bui_threshold[SUBTROPICALzone][EN_DB_TREES] = 110.;
         modelParams.bui_threshold[SUBTROPICALzone][DB_TREES] = 110.;
         modelParams.bui_threshold[TROPICALzone][EN_DB_TREES] = 110.;
         modelParams.bui_threshold[TROPICALzone][DB_TREES] = 110.;
         break;  
             
      case mc2ConUS_LC:
         modelParams.maritime_threshold = 16.;
         modelParams.tmmin_threshold = 1.5; 
                 
         modelParams.ffmc_threshold[BOREALzone][EN_TREES] = 86.0;
         modelParams.bui_threshold[BOREALzone][EN_TREES] =60.0; ;
         modelParams.ffmc_threshold[BOREALzone][DN_TREES] = 86.0;
         modelParams.bui_threshold[BOREALzone][DN_TREES] = 60.0;
         modelParams.ffmc_threshold[BOREALzone][DN_EN_TREES] = 86.0;
         modelParams.bui_threshold[BOREALzone][DN_EN_TREES] = 60.0;
         modelParams.ffmc_threshold[BOREALzone][EN_DB_TREES] = 87.0;
         modelParams.bui_threshold[BOREALzone][EN_DB_TREES] = 110.0;
         modelParams.ffmc_threshold[BOREALzone][DB_TREES] = 87.0;
         modelParams.bui_threshold[BOREALzone][DB_TREES] = 110.0;
         
         modelParams.ffmc_threshold[TEMPERATEzone][EN_TREES] = 89.5;
         modelParams.bui_threshold[TEMPERATEzone][EN_TREES] = 130.0; ;
         modelParams.bui_threshold[TEMPERATEzone][DN_TREES] = 122.85;
         modelParams.bui_threshold[TEMPERATEzone][DN_EN_TREES] = 122.85;         
         modelParams.ffmc_threshold[TEMPERATEzone][EN_DB_TREES] = 86.75;
         modelParams.bui_threshold[TEMPERATEzone][EN_DB_TREES] = 70.0;
         modelParams.ffmc_threshold[TEMPERATEzone][DB_TREES] = 87.5;
         modelParams.bui_threshold[TEMPERATEzone][DB_TREES] = 90.0;
         modelParams.ffmc_threshold[TEMPERATEzone][DB_EB_TREES] = 87.5;
         modelParams.bui_threshold[TEMPERATEzone][DB_EB_TREES] = 80.0;         
         modelParams.ffmc_threshold[TEMPERATEzone][EN_EB_TREES] = 88.;
         modelParams.bui_threshold[TEMPERATEzone][EN_EB_TREES] = 180.0;
         
         modelParams.bui_threshold[SUBTROPICALzone][EN_TREES] = 122.85; 
         modelParams.bui_threshold[SUBTROPICALzone][DN_TREES] = 122.85;
         modelParams.bui_threshold[SUBTROPICALzone][DN_EN_TREES] = 122.85;         
         modelParams.ffmc_threshold[SUBTROPICALzone][EN_DB_TREES] = 87.0;
         modelParams.bui_threshold[SUBTROPICALzone][EN_DB_TREES] = 110.;         
         modelParams.ffmc_threshold[SUBTROPICALzone][DB_TREES] = 86.5;
         modelParams.bui_threshold[SUBTROPICALzone][DB_TREES] = 60.;         
         modelParams.ffmc_threshold[SUBTROPICALzone][DB_EB_TREES] = 85.5;
         modelParams.bui_threshold[SUBTROPICALzone][DB_EB_TREES] = 80.;
         modelParams.ffmc_threshold[SUBTROPICALzone][EN_EB_TREES] = 90.;
         modelParams.bui_threshold[SUBTROPICALzone][EN_EB_TREES] = 310.;
         modelParams.ffmc_threshold[SUBTROPICALzone][EB_TREES] = 86.;
         modelParams.bui_threshold[SUBTROPICALzone][EB_TREES] = 70.;
         
         modelParams.bui_threshold[TROPICALzone][EN_TREES] = 122.85; 
         modelParams.bui_threshold[TROPICALzone][DN_TREES] = 122.85;
         modelParams.bui_threshold[TROPICALzone][DN_EN_TREES] = 122.85;         
         modelParams.ffmc_threshold[TROPICALzone][EN_DB_TREES] = 87.0;
         modelParams.bui_threshold[TROPICALzone][EN_DB_TREES] = 110.;         
         modelParams.ffmc_threshold[TROPICALzone][DB_TREES] = 87.0;
         modelParams.bui_threshold[TROPICALzone][DB_TREES] = 110.;

         changeFRI(Lynx_Boreal_Mixed_Woodland, 25., 25.);
         changeFRI(Lynx_Temperate_Shrubland, 15., 15.);
         changeFRI(Lynx_Subtropical_Shrubland, 15., 15.);
         changeFRI(Lynx_Tropical_Shrubland, 15., 15.);
         changeFRI(Lynx_Tropical_Grassland, 15., 15.);
         changeFRI(Lynx_AgricultureGrazing, 30., 30.);
         break;
         
      case mc2California: // Originally from John Kim's calibration for California at 30 arc-sec resolution
         assert(strcmp(modelParams.century_path, "Input/ModelParameters_MC2_California")==0);

         modelParams.maritime_threshold = 16.;
         modelParams.tmmin_threshold = 1.5; 
         
         modelParams.ffmc_threshold[BOREALzone][EN_TREES] = 86.0;
         modelParams.bui_threshold[BOREALzone][EN_TREES] =60.0; ;
         modelParams.ffmc_threshold[BOREALzone][DN_TREES] = 86.0;
         modelParams.bui_threshold[BOREALzone][DN_TREES] = 60.0;
         modelParams.ffmc_threshold[BOREALzone][DN_EN_TREES] = 86.0;
         modelParams.bui_threshold[BOREALzone][DN_EN_TREES] = 60.0;
         modelParams.ffmc_threshold[BOREALzone][EN_DB_TREES] = 87.0;
         modelParams.bui_threshold[BOREALzone][EN_DB_TREES] = 110.0;
         modelParams.ffmc_threshold[BOREALzone][DB_TREES] = 87.0;
         modelParams.bui_threshold[BOREALzone][DB_TREES] = 110.0;
         
         modelParams.ffmc_threshold[TEMPERATEzone][EN_TREES] = 89.5;
         modelParams.bui_threshold[TEMPERATEzone][EN_TREES] = 130.0; ;
         modelParams.bui_threshold[TEMPERATEzone][DN_TREES] = 122.85;
         modelParams.bui_threshold[TEMPERATEzone][DN_EN_TREES] = 122.85;         
         modelParams.ffmc_threshold[TEMPERATEzone][EN_DB_TREES] = 86.75;
         modelParams.bui_threshold[TEMPERATEzone][EN_DB_TREES] = 70.0;
         modelParams.ffmc_threshold[TEMPERATEzone][DB_TREES] = 87.5;
         modelParams.bui_threshold[TEMPERATEzone][DB_TREES] = 90.0;
         modelParams.ffmc_threshold[TEMPERATEzone][DB_EB_TREES] = 87.5;
         modelParams.bui_threshold[TEMPERATEzone][DB_EB_TREES] = 80.0;         
         modelParams.ffmc_threshold[TEMPERATEzone][EN_EB_TREES] = 88.;
         modelParams.bui_threshold[TEMPERATEzone][EN_EB_TREES] = 180.0;
         
         modelParams.bui_threshold[SUBTROPICALzone][EN_TREES] = 122.85; 
         modelParams.bui_threshold[SUBTROPICALzone][DN_TREES] = 122.85;
         modelParams.bui_threshold[SUBTROPICALzone][DN_EN_TREES] = 122.85;         
         modelParams.ffmc_threshold[SUBTROPICALzone][EN_DB_TREES] = 87.0;
         modelParams.bui_threshold[SUBTROPICALzone][EN_DB_TREES] = 110.;         
         modelParams.ffmc_threshold[SUBTROPICALzone][DB_TREES] = 86.5;
         modelParams.bui_threshold[SUBTROPICALzone][DB_TREES] = 60.;         
         modelParams.ffmc_threshold[SUBTROPICALzone][DB_EB_TREES] = 85.5;
         modelParams.bui_threshold[SUBTROPICALzone][DB_EB_TREES] = 80.;
         modelParams.ffmc_threshold[SUBTROPICALzone][EN_EB_TREES] = 90.;
         modelParams.bui_threshold[SUBTROPICALzone][EN_EB_TREES] = 310.;
         modelParams.ffmc_threshold[SUBTROPICALzone][EB_TREES] = 86.;
         modelParams.bui_threshold[SUBTROPICALzone][EB_TREES] = 70.;
         
         modelParams.bui_threshold[TROPICALzone][EN_TREES] = 122.85; 
         modelParams.bui_threshold[TROPICALzone][DN_TREES] = 122.85;
         modelParams.bui_threshold[TROPICALzone][DN_EN_TREES] = 122.85;         
         modelParams.ffmc_threshold[TROPICALzone][EN_DB_TREES] = 87.0;
         modelParams.bui_threshold[TROPICALzone][EN_DB_TREES] = 110.;         
         modelParams.ffmc_threshold[TROPICALzone][DB_TREES] = 87.0;
         modelParams.bui_threshold[TROPICALzone][DB_TREES] = 110.;
         
         changeFRI(BOREAL_WOODLANDveg, 50., 50.); 
         changeFRI(SUBALPINE_FORESTveg, 250., 250.); // Leenhouts 1998
         changeFRI(MARITIME_EN_FORESTveg, 240., 240.); // Leenhouts 1998
         changeFRI(MESIC_TEMPERATE_NEEDLELEAF_FORESTveg, 76., 76.); // Leenhouts 1998
         changeFRI(TEMPERATE_WARM_MIXED_FORESTveg, 34., 34.); // Leenhouts 1998
         changeFRI(TEMPERATE_EN_WOODLANDveg, 34., 34.); // Leenhouts 1998
         changeFRI(TEMPERATE_DB_WOODLANDveg, 150., 150.); // Krawchuk 2012
         changeFRI(TEMPERATE_COOL_MIXED_WOODLANDveg, 40., 40.); // Leenhouts 1998
         changeFRI(TEMPERATE_WARM_MIXED_WOODLANDveg, 40., 40.); // Leenhouts 1998
         changeFRI(C3SHRUBveg, 120., 120.); // Leenhouts 1998
         changeFRI(C3GRASSveg, 30., 30.); // Leenhouts 1998
         changeFRI(TEMPERATE_DESERTveg, 500., 500.); // Leenhouts 1998
         changeFRI(SUBTROPICAL_EN_FORESTveg, 76., 76.); // same as temperate
         changeFRI(SUBTROPICAL_EN_WOODLANDveg, 50., 50.); // Leenhouts 1998
         changeFRI(SUBTROPICAL_DB_WOODLANDveg, 50., 50.); // Leenhouts 1998
         changeFRI(SUBTROPICAL_EB_WOODLANDveg, 50., 50.); // Leenhouts 1998
         changeFRI(SUBTROPICAL_MIXED_WOODLANDveg, 50., 50.); // Leenhouts 1998
         changeFRI(C4SHRUBveg, 120., 120.); // Leenhouts 1998
         changeFRI(C4GRASSveg, 30., 30.); // Leenhouts 1998
         changeFRI(SUBTROPICAL_DESERTveg, 500., 500.); // Leenhouts 1998
         break;
         
     default: assert(0); break;
      }
      
   modelParams.MAPSSparameterSet = modelParams.interpretMAPSSparameterSet(modelParams.MAPSS_parameter_set);
   modelParams.MAPSS_parameter_set = NULL; // rpStringFcn assumes that if MAPSS_parameter_set is not NULL, then it points to an allocation in the heap
   modelParams.altFuelLoadFlag = interpretSwitch(modelParams.alt_fuel_load_switch);
   modelParams.altTreeAllometryFlag = interpretSwitch(modelParams.alt_tree_allometry_switch);
   modelParams.unlimitedNflag = interpretSwitch(modelParams.unlimited_N_switch);
   modelParams.southernHemisphereFlag = interpretSwitch(modelParams.southern_hemisphere_switch);
   
   return(true);
   } // end of getModelParameters()
   
   
void Simulation::changeFRI(int vtype, int minFRI, int maxFRI)
  {
  assert(modelParams.vveg2mfri[vtype][0]==vtype); 
  modelParams.vveg2mfri[vtype][1] = minFRI;
  modelParams.vveg2mfri[vtype][2] = maxFRI;
  } // end of ChangeFRI()
  

MAPSSparameterSetName ModelParamsClass::interpretMAPSSparameterSet(char * name)
{
  MAPSSparameterSetName mpn;
  
  if (name==NULL) mpn = defaultParams;
  else if (strcmp(name, "ORNLparams")==0) mpn = ORNLparams;
  else if (strcmp(name, "US10kmFAOparams")==0) mpn = US10kmFAOparams;
  else if (strcmp(name, "US10kmSCSparams")==0) mpn = US10kmSCSparams;
  else if (strcmp(name, "WiesYoseParams")==0) mpn = WiesYoseParams;
  else mpn = defaultParams;
  
  return(mpn);
} // end of interpretMAPSSparameterSet()


RunParamsClass::RunParamsClass()
{
  const static char workingDir[] = "./";
  
  first_calendar_year_of_run = 0;
  years_to_run = 1;
  years_offset_into_input_data = 0;
  col_begin = 0;
  col_end = 0;
  row_begin = 0;
  row_end = 0;
  grid_name = NULL; gridName = unknownGridName;
  run_mode = NULL; runMode = unknownRunMode;
  multiyr_start = 1901;
  multiyr_len = 10;
  fire_model_switch = NULL; fireModelFlag = true;
  fire_suppression_switch = NULL; fireSuppressionFlag = false;
  fire_suppression_first_year = 1941;
  suppressed_fire_cell_fraction = 0.0006;
  fire_suppression_fli_threshold = 900.;
  fire_suppression_ros_threshold = 100.;
  fire_suppression_erc_threshold = 60.;
  fire_set_jday = 274;
  fire_set_interval = -1;
  climate_data_directory = (char *)malloc(strlen(workingDir) + 1); assert(climate_data_directory!=NULL); strcpy(climate_data_directory, workingDir);
  earth_data_directory = (char *)malloc(strlen(workingDir) + 1); assert(earth_data_directory!=NULL); strcpy(earth_data_directory, workingDir);
  soil_data_file = (char *)malloc(1); assert(soil_data_file!=NULL); *soil_data_file = 0;
  CO2_file = (char *)malloc(1); assert(CO2_file!=NULL); *CO2_file = 0;
  mask_file =(char *)malloc(1); assert(mask_file!=NULL); *mask_file = 0;
  warmstart_file = (char *)malloc(1); assert(warmstart_file!=NULL); *warmstart_file = 0;
  output_variable_file = (char *)malloc(1); assert(output_variable_file!=NULL); *output_variable_file = 0;
  output_file_prefix = (char *)malloc(1); assert(output_file_prefix!=NULL); *output_file_prefix = 0;
  // model_parameter_set = NULL; 
  baseCalibration = unknownBaseCalibration;
  diags = true;
  spaceBeforeTimeFlag = false;
  monthlyOutputSwitchCount = yearlyOutputSwitchCount = 
         multiyrOutputSwitchCount = timeInvariantOutputSwitchCount =
         warmstartOutputSwitchCount = 0;
  monthlyOutputFlag = yearlyOutputFlag = false;
  multiyrOutputFlag = true;
  timeInvariantOutputFlag = warmstartOutputFlag = false;
  
} // end of default constructor for class RunParamsClass


RunParamsClass::~RunParamsClass()
{
  if (grid_name!=NULL) free(grid_name);
  if (run_mode!=NULL) free(run_mode);
  if (climate_data_directory!=NULL) free(climate_data_directory);
  if (earth_data_directory!=NULL) free(earth_data_directory);
  if (soil_bulk_density_file!=NULL) free(soil_bulk_density_file);
  if (soil_data_file!=NULL) free(soil_data_file);
  if (CO2_file!=NULL) free(CO2_file);
  if (mask_file!=NULL) free(mask_file);
  if (warmstart_file!=NULL) free(warmstart_file);
  if (output_variable_file!=NULL) free(output_variable_file);
  if (output_file_prefix!=NULL) free(output_file_prefix);
  // if (model_parameter_set!=NULL) free(model_parameter_set);
} // end of default destructor for classRunParamsClass
  

ModelParamsClass::ModelParamsClass()
   {
   const static char B41str[] = "B41";
   static char ONstr[] = "ON";
   // const static char modelParametersStr[] = "Input/ModelParameters_B41/";
   
   // Set values of parameters as in MC1 version B41 with WWETAC biogeography
   MAPSS_parameter_set = (char *)malloc(strlen(B41str) + 1); assert(MAPSS_parameter_set!=NULL);
   strcpy(MAPSS_parameter_set, B41str);
   m_century_runoff_x_intercept = 5.f; // that's 5 cm H2O/mo, not millimeters
   m_century_runoff_slope = 0.55f;
   m_lait_lower_limit = 0.2f;
   m_forest_thres_C = 3000.;
   m_woodl_thres_C = 1150.;
   tmmin_threshold = 0.f;
   maritime_threshold = 18.f;
   subalpine_threshold = 1900.f;
   moist_temperate_threshold = 635.f; // 635 mm = 25"
   dry_temperate_threshold = 432.f; // 432 mm = 17"
   p_hi_mult = 1.0f; 
   fire_min = 0.33f;
   alt_fuel_load_switch = ONstr; altFuelLoadFlag = true;
   alt_tree_allometry_switch = ONstr; altTreeAllometryFlag = true;
   unlimited_N_switch = NULL; unlimitedNflag = false;
   efold_t_max = 10.f;
   southern_hemisphere_switch = NULL; southernHemisphereFlag = false;
   
   for (int zone_index = 0; zone_index<NZONES; zone_index++)
      {
      z0[zone_index][GRASS] = 0.0005f;
      z0[zone_index][TREE] = 0.01f;
      z0[zone_index][SHRUB] = 0.001f;
      }  
   z0[mapssTROPICAL][TREE] = 0.02f;
   
   c3_threshold = 55.;
   
   } // end of default constructor for class ModelParamsClass


ModelParamsClass::~ModelParamsClass()
{
  if (MAPSS_parameter_set!=NULL) free(MAPSS_parameter_set);
} // end of default destructor for class ModelParamsClass


BaseCalibrationEnum Simulation::interpretCalibrationName(char * name)
{
  BaseCalibrationEnum bc;
  
/*
  if (strcmp(name, "VEMAP")==0) bc = mc2VEMAP;
  else if (strcmp(name, "NA8km")==0) bc = mc2NA8km;
  else if (strcmp(name, "CA08")==0) bc = mc2CA08;
  else if (strcmp(name, "YOSE")==0) bc = mc2YOSE;
  else if (strcmp(name, "WIES_YOSE")==0) bc = mc2WIES_YOSE;
  else if (strcmp(name, "WWETAC")==0) bc = mc2WWETAC;
  else if (strcmp(name, "VINCERA")==0) bc = mc2VINCERA;
  else if (strcmp(name, "US10kmAlbers")==0) bc = mc2US10kmAlbers;
  else if (strcmp(name, "US50km")==0) bc = mc2US50km;
  else ...
*/
  if (strcmp(name, "GLOBAL")==0) bc = mc2GLOBAL;
  else if (strcmp(name, "CONUS")==0) bc = mc2ConUS;
  else if (strcmp(name, "CONUS_LC")==0) bc = mc2ConUS_LC;
  else if (strcmp(name, "WCR")==0) bc = mc2WCR;
  else if (strcmp(name, "CALIFORNIA")==0) bc = mc2California;
  else bc = unknownBaseCalibration;
  
  return(bc);
} // end of Simulation::interpretCalibrationName()


void err_exit(const char * msg)
   {
   printf("\n*** " COMMAND_NAME ": %s\n", msg);
   exit(-1);
   } // end of err_exit


double Simulation::north_south_centroid(int row)
   {
   double cell_lat;
   
   cell_lat = runParams.latitude0 - runParams.cell_spacing*(row + 0.5);
   return(cell_lat);
   } // end of north_south_centroid()


double Simulation::east_west_centroid(int col)
   {
   double cell_lon;
   
   cell_lon = runParams.longitude0 + runParams.cell_spacing*(col + 0.5);
   return(cell_lon);
   } // end of east_west_centroid()


bool chk_nc(int io_rtnval)
   {
   if (io_rtnval==NC_NOERR) return(TRUE);
   printf("*** %s\n", nc_strerror(io_rtnval));
   return(FALSE);
   } // end of chk_nc()
   
   
void Simulation::saveOneYearOneCellOutputData(int yrNdx)
   {
   if (m_num_yearly_vars<=0) return;
   assert(m_yrBufP!=NULL);
   int startPos = yrNdx*m_num_yearly_vars;
   for (int outNdx = 0; outNdx<m_num_yearly_vars; outNdx++) 
         *(m_yrBufP + startPos + outNdx) = m_yrOutVars[outNdx];
         //*(m_yrBufP + startPos + outNdx) = *(m_yrOutList[outNdx].valP);
   }; // end of Simulation::saveOneYearOneCellOutputData()


void Simulation::saveOneMonthOneCellOutputData(int moNdx)
   {
   if (m_num_monthly_vars<=0) return;
   assert(m_moBufP!=NULL);
   void * destination = (void *)(m_moBufP + moNdx*m_num_monthly_vars);
   memcpy(destination, (void *)m_moOutVars, m_num_monthly_vars*sizeof(float));

   }; // end of Simulation::saveOneMonthOneCellOutputData()


void OutVarDescription::save(float inVal, int interval)
   { 
   if (!show) return;
   switch (interval)
      {
      case MONTH_INTERVAL: 
         assert(moOutNdx>=0);
         m_pS->m_moOutVars[moOutNdx] = inVal; 
         break;
      case YEAR_INTERVAL: 
         assert(yrOutNdx>=0);
         m_pS->m_yrOutVars[yrOutNdx] = inVal; 
         break;
      case MULTIYR_INTERVAL: 
         assert(multiyrOutNdx>=0);
         m_pS->m_multiyrOutVars[multiyrOutNdx] = inVal; 
         break;
      case SINGLE_INTERVAL: 
         assert(singleOutNdx>=0);
         m_pS->m_singleOutVars[singleOutNdx] = inVal; 
         break;
      default: assert(0); break;
      }
   } // end of OutVarDescription::save(inVal, interval) method
   

void OutVarDescription::save(float inVal)
   {
   assert(interval==SINGLE_INTERVAL);
   save(inVal, SINGLE_INTERVAL);
   } // end of OutVarDescription::save(float inVal)
   

void Simulation::writeMonthlyData(size_t coords[])
   {
   int outNdx;
   int num_months = 12*runParams.actual_years_to_run;

   if (!runParams.monthlyOutputFlag) return;
   
   assert(m_num_monthly_vars>0);

   if (runParams.singleCellFlag) 
      { // just one cell - write out the data in text file format
      printf("monthly output variables");
      for (outNdx = 0; outNdx<m_num_monthly_vars; outNdx++) printf(", %s",
            m_moOutList[outNdx].name);
      printf("\n");
      for (int mo = 0; mo<num_months; mo++)
         { 
         outNdx = 0;
         printf("%d", mo);
         while (outNdx<m_num_monthly_vars) 
            {
            printf(", %f", *(m_moBufP + mo*m_num_monthly_vars + outNdx));
            outNdx++;
            }   
         printf("\n");
         }
      printf("\n");
      }
   else 
      { // more than one cell - write out the data to the _month.nc netCDF file
      outNdx = 0;
      while (outNdx<m_num_monthly_vars) 
         { int rtn_val;
         for (int mo = 0; mo<num_months; mo++)
            { float floatval;
            coords[m_moOutFile_ncid.timeNdxPos] = mo;
            floatval = *(m_moBufP + mo*m_num_monthly_vars + outNdx);
            rtn_val = nc_put_var1_float(m_moOutFile_ncid.fileid, m_moOutList[outNdx].var_id, 
                  coords, &floatval); 
            if (!chk_nc(rtn_val))
                  assert(false);
            }
         outNdx++;
         }
      }
   } // end of Simulation::writeMonthlyData()


void Simulation::writeYearlyData(size_t coords[])
   {
   int outNdx;

    if (!runParams.yearlyOutputFlag) return;
    
    assert(m_num_yearly_vars>0);
   
   if (runParams.singleCellFlag) 
      { // just one cell - write out the data in text file format
      printf("yearly output variables");
      for (outNdx = 0; outNdx<m_num_yearly_vars; outNdx++) printf(", %s",
            m_yrOutList[outNdx].name);
      printf("\n");
      for (int yrNdx = 0; yrNdx<runParams.actual_years_to_run; yrNdx++)
         { 
         outNdx = 0;
         printf("%d", yrNdx);
         while (outNdx<m_num_yearly_vars) 
            {
            printf(", %f", *(m_yrBufP + yrNdx*m_num_yearly_vars + outNdx));
            outNdx++;
            }   
         printf("\n");
         }
      printf("\n");
      }
   else
      { // more than one cell - write out the data to the _year.nc netCDF file
      outNdx = 0;
      while (outNdx<m_num_yearly_vars) 
         { int rtn_val;
         for (int yrNdx = 0; yrNdx<runParams.actual_years_to_run; yrNdx++)
            { float floatval;
            coords[m_yrOutFile_ncid.timeNdxPos] = yrNdx;
            floatval = *(m_yrBufP + yrNdx*m_num_yearly_vars + outNdx);
            rtn_val = nc_put_var1_float(m_yrOutFile_ncid.fileid, m_yrOutList[outNdx].var_id, 
                  coords, &floatval); 
            if (!chk_nc(rtn_val))
                  assert(false);
            }
         outNdx++;
         }
      }
   } // end of Simulation::writeYearlyData()


void Simulation::writeMultiyrData(size_t coords[])
   {
   int outNdx;

   if (!runParams.multiyrOutputFlag) return;
    
   assert(m_num_multiyr_vars>0);
   
   if (runParams.singleCellFlag) 
      { // just one cell - write out the data in text file format
      printf("multiyr output variables");
      for (outNdx = 0; outNdx<m_num_multiyr_vars; outNdx++) printf(", %s",
            m_multiyrOutList[outNdx].name);
      printf("\n");
      for (int multiyrNdx = 0; multiyrNdx<runParams.num_multiyr_intervals; multiyrNdx++)
         { 
         outNdx = 0;
         printf("%d", multiyrNdx);
         while (outNdx<m_num_multiyr_vars) 
            {
            printf(", %f", *(m_multiyrBufP + multiyrNdx*m_num_multiyr_vars + outNdx));
            outNdx++;
            }   
         printf("\n");
         }
      printf("\n");
      }
   else
      { // more than one cell - write out the data to the _multiyr.nc netCDF file
      outNdx = 0;
      while (outNdx<m_num_multiyr_vars) 
         { int rtn_val;
         for (int multiyrNdx = 0; multiyrNdx<runParams.num_multiyr_intervals; multiyrNdx++)
            { float floatval;
            coords[m_multiyrOutFile_ncid.timeNdxPos] = multiyrNdx;
            floatval = *(m_multiyrBufP + multiyrNdx*m_num_multiyr_vars + outNdx);
            rtn_val = nc_put_var1_float(m_multiyrOutFile_ncid.fileid, m_multiyrOutList[outNdx].var_id, 
                  coords, &floatval); 
            if (!chk_nc(rtn_val))
                  assert(false);
            }
         outNdx++;
         }
      }
   } // end of Simulation::writeMultiyrData()


void Simulation::writeSingleData(size_t coords[])
   {
   int outNdx;

   if (!runParams.timeInvariantOutputFlag) return;
   
   assert(m_num_single_vars>0);
   
   if (runParams.singleCellFlag) 
      { // just one cell - write out the data in text file format
      printf("time-invariant (..._single.nc) output variables");
      for (outNdx = 0; outNdx<m_num_single_vars; outNdx++) printf(", %s %f\n",
            m_singleOutList[outNdx].name, m_singleOutVars[outNdx]);
      printf("\n");
      }
   else
      { // more than one cell - write out the data to the _single.nc netCDF file
      if (m_singleOutFile_ncid.timeNdxPos>=0) coords[m_singleOutFile_ncid.timeNdxPos] = 0;
      outNdx = 0;
      while (outNdx<m_num_single_vars) 
         { int rtn_val;
         rtn_val = nc_put_var1_float(m_singleOutFile_ncid.fileid, m_singleOutList[outNdx].var_id, 
               coords, m_singleOutVars + outNdx); 
         if (!chk_nc(rtn_val))
               assert(false);
         outNdx++;
         }
      }
   } // end of Simulation::writeSingleData()


void Simulation::saveMultiyrOutputData()
// Called once at the end of the spinup and transient phases.
   {
   int offset = 0;
   for ( int interval_index = 0; interval_index<runParams.num_multiyr_intervals; interval_index++)
      { 
      int first_yr_of_interval, last_yr_of_interval, current_interval_len;
      float * recordP;
      
      if (runParams.runMode==SPINUP && interval_index==0)
         {
         first_yr_of_interval = last_yr_of_interval = 0;
         current_interval_len = 1;
         offset = 1;
         }
      else
         {
         first_yr_of_interval = runParams.nominal_multiyr_start + (interval_index - offset)*runParams.multiyr_len;
         last_yr_of_interval = first_yr_of_interval + runParams.multiyr_len - 1;
         if (first_yr_of_interval<runParams.first_calendar_year_of_run) 
               first_yr_of_interval = runParams.first_calendar_year_of_run;
         if (last_yr_of_interval>runParams.last_calendar_year_of_run) 
               last_yr_of_interval = runParams.last_calendar_year_of_run;
         current_interval_len = last_yr_of_interval - first_yr_of_interval + 1;
         }
      assert(current_interval_len>=1);
      *(m_multiyrIntervalBufP + 2*interval_index) = last_yr_of_interval;
      *(m_multiyrIntervalBufP + 2*interval_index + 1) = current_interval_len;
      recordP = m_multiyrBufP + interval_index*m_num_multiyr_vars;
      
      for (int outNdx = 0; outNdx<m_num_multiyr_vars; outNdx++)
         { 
         int index_in_yrOutList = VarDict[m_multiyrOutList[outNdx].dictNdx].yrOutNdx;
         assert(index_in_yrOutList>=0 && index_in_yrOutList<m_num_yearly_vars);
         int first_yrNdx = first_yr_of_interval - runParams.first_calendar_year_of_run;
         int last_yrNdx = last_yr_of_interval - runParams.first_calendar_year_of_run;
         if (m_multiyrOutList[outNdx].categoricalFlag)
            { // Calculate modal value over the interval.
            int * val_countP;
            
            int max_cat = m_multiyrOutList[outNdx].max_categories + 1; // allow for 0,1,...,max_categories
            assert(max_cat>0);
            val_countP = (int *)malloc(max_cat*sizeof(int)); // free'd just before the end of the block
            for (int catNdx = 0; catNdx<max_cat; catNdx++) *(val_countP + catNdx) = 0;
            for (int yrNdx = first_yrNdx; yrNdx<=last_yrNdx; yrNdx++)
               { int cat;
               cat = (int)(*(m_yrBufP + yrNdx*m_num_yearly_vars + outNdx));
               if (cat>=0 && cat<max_cat) *(val_countP + cat) += 1;
               }
            int modal_val = -9999;
            int modal_count = 0;
            for (int catNdx = 0; catNdx<max_cat; catNdx++) if (*(val_countP + catNdx)>modal_count)
               {
               modal_val = catNdx;
               modal_count = *(val_countP + catNdx);
               }
            *(recordP + outNdx) = modal_val;
            free(val_countP);
            } // end of block to calculate modal value
         else
            { // Calculate mean value over the interval.
            float sum = 0.;
            for (int yrNdx = first_yrNdx; yrNdx<=last_yrNdx; yrNdx++)
                  sum += *(m_yrBufP + yrNdx*m_num_yearly_vars + outNdx);
            *(recordP + outNdx) = sum/current_interval_len;
            }
         } // end of loop over all the multiyr variables         
      } // end of loop over all the intervals
   } // end of Simulation::saveMultiyrOutputData()


int Simulation::initializeActiveCellArray()
   { 
   int i, j;
   int nRows = runParams.row_end - runParams.row_begin + 1;
   int nCols = runParams.col_end - runParams.col_begin + 1;
   
   if (mask_ncid.fileid==-1) 
      {
      for (i = 0; i<nRows; i++) for (j = 0; j<nCols; j++)
            *(m_pActiveCellArray + i*nCols + j) = true;
      return(nRows*nCols);
      }
      
   size_t coords[3], count[3];
   short * shortP;
   signed char * scharP;
   float * floatP;
   int * intP;
   int rtnval;
   
   coords[mask_ncid.rowNdxPos] = runParams.row_begin - mask_ncid.row_offset;
   coords[mask_ncid.colNdxPos] = runParams.col_begin - mask_ncid.col_offset;
   if (mask_ncid.timeNdxPos>=0) coords[mask_ncid.timeNdxPos] = 0;
   
   count[mask_ncid.rowNdxPos] = nRows;
   count[mask_ncid.colNdxPos] = nCols;
   if (mask_ncid.timeNdxPos>=0) count[mask_ncid.timeNdxPos] = 1;
  
   floatP = (float *)malloc(nRows*nCols*sizeof(float)); assert(floatP!=NULL);
   switch (mask_ncid.input_data_vartype)
      {
      case NC_INT:
         intP = (int *)malloc(nRows*nCols*sizeof(int)); assert(intP!=NULL);
         rtnval = nc_get_vara_int(mask_ncid.fileid, mask_ncid.input_data_varid, coords, count, intP);
         assert(chk_nc(rtnval));
         for (i = 0; i<nRows; i++) for (j = 0; j<nCols; j++) 
               *(floatP + i*nCols + j) = *(intP + i*nCols + j);
         free(intP);
         break;
      case NC_SHORT:
         shortP = (short *)malloc(nRows*nCols*sizeof(short)); assert(shortP!=NULL);
         rtnval = nc_get_vara_short(mask_ncid.fileid, mask_ncid.input_data_varid, coords, count, shortP);
         assert(chk_nc(rtnval));
         for (i = 0; i<nRows; i++) for (j = 0; j<nCols; j++) 
               *(floatP + i*nCols + j) = *(shortP + i*nCols + j);
         free(shortP);
         break;
      case NC_BYTE:
         scharP = (signed char *)malloc(nRows*nCols*sizeof(signed char)); assert(scharP!=NULL);
         rtnval = nc_get_vara_schar(mask_ncid.fileid, mask_ncid.input_data_varid, coords, count, scharP);
         assert(chk_nc(rtnval));
         for (i = 0; i<nRows; i++) for (j = 0; j<nCols; j++) 
               *(floatP + i*nCols + j) = *(scharP + i*nCols + j);
         free(scharP);
         break;
      case NC_FLOAT:
         rtnval = nc_get_vara_float(mask_ncid.fileid, mask_ncid.input_data_varid, coords, count, floatP);
         assert(chk_nc(rtnval));
         break;
      default:
         assert(false);
         break;
   }
   
  int num_target_cells = 0;
  for (i = 0; i<nRows; i++) for (j = 0; j<nCols; j++) 
    { bool activeFlag;
    activeFlag = *(floatP + i*nCols + j)>0;
    if (modelParams.code_flags[MASK_GE_997000_CODE_FLAG] && *(floatP + i*nCols + j)>=997000.) activeFlag = false;
    *(m_pActiveCellArray + i*nCols + j) = activeFlag;
    if (activeFlag) num_target_cells++;
    }
  free(floatP);  
  return(num_target_cells); 
  } // end of Simulation::initializeActiveCellArray()
  
  


