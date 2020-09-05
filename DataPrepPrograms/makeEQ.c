/*
 *  makeEQ.c
 *  DataPrepPrograms
 *
 *  Created by Dave Conklin on 7/10/10.
 *  Copyright 2010 Conservation Biology Institute. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <assert.h>
#include <math.h>
#include "netcdf.h"
#include "gridspecs.h"

#define Boolean int
#define TRUE 1
#define FALSE 0
#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)

#define STRING_LEN 200

#define COMMAND_NAME "makeEQ"
#define NAME_ARG 1
#define INPUT_PATH_ARG NAME_ARG+1
#define OUTPUT_PATH_ARG INPUT_PATH_ARG+1
#define TIME_NAME_ARG OUTPUT_PATH_ARG+1
#define FIRST_INPUT_DATA_YEAR_ARG TIME_NAME_ARG+1
#define NOTE_ARG FIRST_INPUT_DATA_YEAR_ARG+1
#define NARGS NOTE_ARG+1

Boolean chk_nc(int io_rtnval);
void err_exit(char * msg);
void instructions();
float * floatP(float * base, int time, short row, short col, short nrows, short ncols);
short * shortP(short * base, int time, short row, short col, short nrows, short ncols);


int main(int argc, char * argv[])
{
  int iarg;
  char varname[STRING_LEN+1], infilename[STRING_LEN+1], outfilename[STRING_LEN+1];
  int infileid, invarid;
  int latdimid, londimid, timedimid;
  int ndims;
  int indimids[NC_MAX_VAR_DIMS];
  size_t NumRows, NumCols;
  nc_type invartype;
  int first_year, last_year, first_input_data_year, NumYrsAv, first_input_layer;
  Boolean timeless_flag, month_flag;
  int NumInputLayers, valCount;
  short * InputShortDataValsP;
  float * InputFloatDataValsP;
  size_t incoords[3], incounts[3], outcoords[3], outcounts[3];
  int row, col, ThisMo, ThisYr;
  long int Sum;
  float Mean;
  short shortVal;
  short * EQValsShortP;
  float floatVal;
  float * EQValsFloatP;
  int outfileid, outvarid, outbanddimid, outlatdimid, outlondimid, outdimids[3];
  int outlatvarid, outlonvarid, inlatvarid, inlonvarid;
  float * lat_arrayInP;
  float * lon_arrayInP;
  char info[STRING_LEN + 1];
  int rtnval;
  Boolean twoD_flag;
  
  if (argc!=NARGS) 
  {
    printf("*** makeEQ: wrong number of arguments.\n");
    instructions();   
  }
  for (iarg = 0; iarg<argc; iarg++) printf("%s ", argv[iarg]);
  printf("\n");
  
  // Open the historical data netCDF file.
  strncpy(varname, argv[NAME_ARG], STRING_LEN); varname[STRING_LEN] = 0;
  strncpy(infilename, argv[INPUT_PATH_ARG], STRING_LEN); infilename[STRING_LEN] = 0;
  strncat(infilename, varname, sizeof(infilename) - strlen(infilename) - 1); infilename[STRING_LEN] = 0;
  strncat(infilename, ".nc", sizeof(infilename) - strlen(infilename) - 1); infilename[STRING_LEN] = 0;
  rtnval = nc_open(infilename, NC_NOWRITE, &infileid); 
  if (!chk_nc(rtnval)) 
  {
    printf("infilename = %s\n", infilename);
    err_exit("error while opening input data file.");
  }


  rtnval = nc_inq_dimid(infileid, "time", &timedimid); 
  if (rtnval!=NC_NOERR) 
  {
    rtnval = nc_inq_dimid(infileid, "year", &timedimid); 
    if (rtnval!=NC_NOERR) 
    {
      rtnval = nc_inq_dimid(infileid, "month", &timedimid); assert(chk_nc(rtnval));
    }
  }
  rtnval = nc_inq_dimid(infileid, "lat", &latdimid); assert(chk_nc(rtnval));
  rtnval = nc_inq_dimid(infileid, "lon", &londimid); assert(chk_nc(rtnval));
  
  rtnval = nc_inq_dimlen(infileid, latdimid, &NumRows); assert(chk_nc(rtnval));  
  rtnval = nc_inq_dimlen(infileid, londimid, &NumCols); assert(chk_nc(rtnval));
  rtnval = nc_inq_varid(infileid, varname, &invarid); assert(chk_nc(rtnval));
  rtnval = nc_inq_vartype(infileid, invarid, &invartype); assert(chk_nc(rtnval));
  assert(invartype==NC_SHORT || invartype==NC_FLOAT);

  rtnval = nc_inq_varndims(infileid, invarid, &ndims); assert(chk_nc(rtnval));
  assert(ndims==3);
  rtnval = nc_inq_vardimid(infileid, invarid, indimids); assert(chk_nc(rtnval));
  assert(indimids[0]==timedimid);
  assert(indimids[1]==latdimid);
  assert(indimids[2]==londimid);

  getTime(argv[TIME_NAME_ARG], &first_year, &last_year, &timeless_flag, &month_flag, &twoD_flag);
  assert(!timeless_flag && !month_flag && !twoD_flag);
  NumYrsAv = last_year - first_year + 1;
  
  if (sscanf(argv[FIRST_INPUT_DATA_YEAR_ARG], "%d", &first_input_data_year)!=1) 
      err_exit("error when reading <first year of input data>.");
      
  first_input_layer = 12*(first_year - first_input_data_year); assert(first_input_layer>=0);
  NumInputLayers = 12*NumYrsAv;

  printf("Getting values from input netCDF file to make EQ...\n");
  valCount = NumCols*NumRows*NumInputLayers;
  incoords[0] = first_input_layer; incoords[1] = 0; incoords[2] = 0;
  incounts[0] = NumInputLayers; incounts[1] = NumRows; incounts[2] = NumCols;
  switch (invartype)
  {
    case NC_SHORT:
      InputShortDataValsP = malloc(sizeof(short)*valCount); assert(InputShortDataValsP!=NULL);
      rtnval = nc_get_vara_short(infileid, invarid, incoords, incounts, InputShortDataValsP);
      assert(chk_nc(rtnval));
      break;
    case NC_FLOAT:
      InputFloatDataValsP = malloc(sizeof(float)*valCount); assert(InputFloatDataValsP!=NULL);
      rtnval = nc_get_vara_float(infileid, invarid, incoords, incounts, InputFloatDataValsP);
      assert(chk_nc(rtnval));
      break;
    default:
      assert(0);
      break;
  }
      

  // Calculate the EQ layer
  valCount = NumCols*NumRows*12;
  switch (invartype)
  {
    case NC_SHORT:
      EQValsShortP = (short *)malloc(sizeof(short)*valCount); assert(EQValsShortP!=NULL);
      printf("Starting main loop now...\n");
      for (row = 0; row<NumRows; row++) 
      {
        printf("Starting row %d now...\n", row);
        for (col = 0; col<NumCols; col++)
        {
          if (*shortP(InputShortDataValsP, 0, row, col, NumRows, NumCols)==-9999)
              for (ThisMo = 0; ThisMo<12; ThisMo++) 
              *shortP(EQValsShortP, ThisMo, row, col, NumRows, NumCols) = -9999;
          else 
          {
            for (ThisMo = 0; ThisMo<12; ThisMo++) 
            {
              Sum = 0;
              for (ThisYr = 0; ThisYr<NumYrsAv; ThisYr++)
              {
                shortVal = *shortP(InputShortDataValsP, ThisYr*12 + ThisMo, row, col, NumRows, NumCols);
                assert(shortVal!=-9999);
                Sum += (long int)shortVal;
              } // end of for ThisYr
              Mean = ((float)Sum)/NumYrsAv;
              shortVal = Mean<0. ? (short)(Mean - 0.5) : (short)(Mean + 0.5);
              *shortP(EQValsShortP, ThisMo, row, col, NumRows, NumCols) = shortVal;
            } // end of for ThisMo  
          }
        } // end of for col
      } // end of for row 
      break;
    case NC_FLOAT:
      EQValsFloatP = (float *)malloc(sizeof(float)*valCount); assert(EQValsFloatP!=NULL);
      printf("Starting main loop now...\n");
      for (row = 0; row<NumRows; row++) 
      {
        printf("Starting row %d now...\n", row);
        for (col = 0; col<NumCols; col++)
        {
          if (*floatP(InputFloatDataValsP, 0, row, col, NumRows, NumCols)==1.E20)
              for (ThisMo = 0; ThisMo<12; ThisMo++) 
              *floatP(EQValsFloatP, ThisMo, row, col, NumRows, NumCols) = 1.E20;
          else 
          {
            for (ThisMo = 0; ThisMo<12; ThisMo++) 
            {
              Sum = 0;
              for (ThisYr = 0; ThisYr<NumYrsAv; ThisYr++)
              {
                floatVal = *floatP(InputFloatDataValsP, ThisYr*12 + ThisMo, row, col, NumRows, NumCols);
                assert(floatVal!=1.E20);
                Sum += (double)floatVal;
              } // end of for ThisYr
              Mean = ((double)Sum)/NumYrsAv;
              *floatP(EQValsFloatP, ThisMo, row, col, NumRows, NumCols) = Mean;
            } // end of for ThisMo  
          }
        } // end of for col
      } // end of for row 
      break;
    default:
      assert(0);
      break;
  } // end of switch(invartype)
      
  // Create the output file.  
  strncpy(outfilename, argv[OUTPUT_PATH_ARG], STRING_LEN); outfilename[STRING_LEN] = 0;
  strncat(outfilename, varname, sizeof(outfilename) - strlen(outfilename) - 1); outfilename[STRING_LEN] = 0;
  strncat(outfilename, ".nc", sizeof(outfilename) - strlen(outfilename) - 1); outfilename[STRING_LEN] = 0;
  assert(strcmp(infilename, outfilename)!=0);
  printf("Creating output netCDF file: %s\n", outfilename);	
  rtnval = nc_create(outfilename, NC_CLOBBER | NC_64BIT_OFFSET, &outfileid); assert(chk_nc(rtnval));
  rtnval = nc_def_dim(outfileid, "month", NC_UNLIMITED, &outbanddimid); assert(chk_nc(rtnval));
  rtnval = nc_def_dim(outfileid, "lat", NumRows, &outlatdimid); assert(chk_nc(rtnval));
  rtnval = nc_def_dim(outfileid, "lon", NumCols, &outlondimid); assert(chk_nc(rtnval));
  rtnval = nc_copy_att(infileid, NC_GLOBAL, "row_offset", outfileid, NC_GLOBAL); 
  rtnval = nc_copy_att(infileid, NC_GLOBAL, "col_offset", outfileid, NC_GLOBAL);
  rtnval = nc_copy_att(infileid, NC_GLOBAL, "Conventions", outfileid, NC_GLOBAL); 
  nc_copy_att(infileid, NC_GLOBAL, "title", outfileid, NC_GLOBAL);
  nc_copy_att(infileid, NC_GLOBAL, "institution", outfileid, NC_GLOBAL);
  nc_copy_att(infileid, NC_GLOBAL, "source", outfileid, NC_GLOBAL);
  nc_copy_att(infileid, NC_GLOBAL, "contact", outfileid, NC_GLOBAL);
  nc_copy_att(infileid, NC_GLOBAL, "project_id", outfileid, NC_GLOBAL);
  nc_copy_att(infileid, NC_GLOBAL, "table_id", outfileid, NC_GLOBAL);
  nc_copy_att(infileid, NC_GLOBAL, "experiment_id", outfileid, NC_GLOBAL);
  nc_copy_att(infileid, NC_GLOBAL, "realization", outfileid, NC_GLOBAL);
  nc_copy_att(infileid, NC_GLOBAL, "cmor_version", outfileid, NC_GLOBAL);
  nc_copy_att(infileid, NC_GLOBAL, "history", outfileid, NC_GLOBAL);
  nc_copy_att(infileid, NC_GLOBAL, "references", outfileid, NC_GLOBAL);
  nc_copy_att(infileid, NC_GLOBAL, "comment", outfileid, NC_GLOBAL);
  info[0] = 0;
  for (iarg=0; iarg<(argc - 1); iarg++)
  {
    if (iarg==NOTE_ARG) continue;
    strncat(info, argv[iarg], STRING_LEN - strlen(info));
    strncat(info, " ", STRING_LEN - strlen(info));
  }  
  rtnval = nc_put_att_text(outfileid, NC_GLOBAL, "cmdline", strlen(info), info); assert(chk_nc(rtnval));
  nc_copy_att(infileid, NC_GLOBAL, "time_series", outfileid, NC_GLOBAL);
  rtnval = nc_put_att_text(outfileid, NC_GLOBAL, "note", strlen(argv[NOTE_ARG]), argv[NOTE_ARG]); assert(chk_nc(rtnval));
  outvarid = -1; // Prevent a compiler warning about an uninitialized variable
  nc_copy_att(infileid, NC_GLOBAL, "note1", outfileid, outvarid);
  nc_copy_att(infileid, NC_GLOBAL, "note2", outfileid, outvarid);
  nc_copy_att(infileid, NC_GLOBAL, "note3", outfileid, outvarid);
  nc_copy_att(infileid, NC_GLOBAL, "note4", outfileid, outvarid);
  nc_copy_att(infileid, NC_GLOBAL, "note5", outfileid, outvarid);
  nc_copy_att(infileid, NC_GLOBAL, "note6", outfileid, outvarid);
  nc_copy_att(infileid, NC_GLOBAL, "note7", outfileid, outvarid);
  nc_copy_att(infileid, NC_GLOBAL, "note8", outfileid, outvarid);  
  
  outdimids[0] = outbanddimid;
  outdimids[1] = outlatdimid;
  outdimids[2] = outlondimid;
  rtnval = nc_def_var(outfileid, varname, invartype, 3, outdimids, &outvarid); assert(chk_nc(rtnval));  
  nc_copy_att(infileid, invarid, "missing_value", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "_FillValue", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "units", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "scale_factor", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "valid_min", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "valid_max", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "long_name", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "scaled", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "band1", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "band2", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "band3", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "band4", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "band5", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "band6", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "band7", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "band8", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "band9", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "band10", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "standard_name", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "original_name", outfileid, outvarid);
  
           
  // Create the coordinate dimension variables.
  outdimids[0] = outlatdimid;
  rtnval = nc_def_var(outfileid, "lat", NC_FLOAT, 1, outdimids, &outlatvarid); assert(chk_nc(rtnval));
  rtnval = nc_put_att_text(outfileid, outlatvarid, "long_name", strlen("latitude"), "latitude"); assert(chk_nc(rtnval));
  rtnval = nc_put_att_text(outfileid, outlatvarid, "standard_name", strlen("latitude"), "latitude"); assert(chk_nc(rtnval));
  rtnval = nc_put_att_text(outfileid, outlatvarid, "units", strlen("degrees_north"), "degrees_north"); assert(chk_nc(rtnval));
  outdimids[0] = outlondimid;
  rtnval = nc_def_var(outfileid, "lon", NC_FLOAT, 1, outdimids, &outlonvarid); assert(chk_nc(rtnval));
  rtnval = nc_put_att_text(outfileid, outlonvarid, "long_name", strlen("longitude"), "longitude"); assert(chk_nc(rtnval));
  rtnval = nc_put_att_text(outfileid, outlonvarid, "standard_name", strlen("longitude"), "longitude"); assert(chk_nc(rtnval));
  rtnval = nc_put_att_text(outfileid, outlonvarid, "units", strlen("degrees_east"), "degrees_east"); assert(chk_nc(rtnval));

  rtnval = nc_enddef(outfileid); assert(chk_nc(rtnval));
   
  lat_arrayInP = (float *)malloc(sizeof(float)*NumRows); assert(lat_arrayInP!=NULL);
  rtnval = nc_inq_varid(infileid, "lat", &inlatvarid); assert(chk_nc(rtnval));
  rtnval = nc_get_var_float(infileid, inlatvarid, lat_arrayInP); assert(chk_nc(rtnval));
  rtnval = nc_put_var_float(outfileid, outlatvarid, lat_arrayInP); assert(chk_nc(rtnval));
  
  lon_arrayInP = (float *)malloc(sizeof(float)*NumCols); assert(lon_arrayInP!=NULL);
  rtnval = nc_inq_varid(infileid, "lon", &inlonvarid); assert(chk_nc(rtnval));
  rtnval = nc_get_var_float(infileid, inlonvarid, lon_arrayInP); assert(chk_nc(rtnval));
  rtnval = nc_put_var_float(outfileid, outlonvarid, lon_arrayInP); assert(chk_nc(rtnval));

  // Now write the EQ data array to the netCDF file.
  outcoords[1] = outcoords[2] = 0;
  outcounts[0] = 1; outcounts[1] = NumRows; outcounts[2] = NumCols;
  for (ThisMo = 0; ThisMo<12; ThisMo++)
  {
    printf("writing out layer for month %d now...\n", ThisMo);
    outcoords[0] = ThisMo; 
    switch (invartype)
    {
      case NC_SHORT:
        rtnval = nc_put_vara_short(outfileid, outvarid, outcoords, outcounts, 
            shortP(EQValsShortP, ThisMo, 0, 0, NumRows, NumCols)); assert(chk_nc(rtnval));
        break;
      case NC_FLOAT:
        rtnval = nc_put_vara_float(outfileid, outvarid, outcoords, outcounts, 
            floatP(EQValsFloatP, ThisMo, 0, 0, NumRows, NumCols)); assert(chk_nc(rtnval));
        break;
      default:
        assert(0);
        break;
    }
  }
  
  rtnval = nc_close(outfileid); assert(chk_nc(rtnval));  
  printf(COMMAND_NAME " has completed.\n");
  return 0;
  
} // end of main()


Boolean chk_nc(int io_rtnval)
{
  if (io_rtnval==NC_NOERR) return(TRUE);
  printf("*** %s\n", nc_strerror(io_rtnval));
  return(FALSE);
} // end of chk_nc()


void err_exit(char * msg)
{
  printf("\n*** " COMMAND_NAME ": %s\n", msg);
  exit(-1);
} // end of err_exit


void instructions()
{
  printf("\n*** " COMMAND_NAME
      ": Invoke as '" COMMAND_NAME " <variable name> <input data path> <EQ output data path> "
          "<time name> <first year of input data> <note>\n"
	    "  Example: " COMMAND_NAME " tmin Input/US800m/HistoricalOregon/ ./ 1895_1950 1895 "
      "'tmin.nc EQ file for Oregon based on the 2009 edition of the PRISM US800m historical climate dataset, "
      "made by David Conklin on 7/10/10 using the " COMMAND_NAME ".c program.'\n"); 
  exit(-1);
} // end of instructions()


short * shortP(short * base, int time, short row, short col, short nrows, short ncols)
{
  short * finalP;
  
  finalP = base + (long)time*(long)nrows*(long)ncols + (long)row*(long)ncols + (long)col;
  return(finalP);
} // end of shortP()


float * floatP(float * base, int time, short row, short col, short nrows, short ncols)
{
  float * finalP;
  
  finalP = base + (long)time*(long)nrows*(long)ncols + (long)row*(long)ncols + (long)col;
  return(finalP);
} // end of floatP()


