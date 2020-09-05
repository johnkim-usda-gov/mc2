/*
 *  makeSpinup.c - makes the spinup data from the historical data
 *  DataPrepPrograms
 *
 *  Created by Dave Conklin on 7/9/10.
 *  Copyright 2010 Conservation Biology Institute. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <assert.h>
#include <math.h>
#include "netcdf.h"

#define Boolean int
#define TRUE 1
#define FALSE 0

#define COMMAND_NAME "makeSpinup"

#define Boolean int
#define TRUE 1
#define FALSE 0
#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)

#define STRING_LEN 200

#define COMMAND_NAME "makeSpinup"
#define NOTE_ARG 1
#define INPUT_PATH_ARG NOTE_ARG+1
#define EQ_PATH_ARG INPUT_PATH_ARG+1
#define OUTPUT_PATH_ARG EQ_PATH_ARG+1
#define MIN_ARG OUTPUT_PATH_ARG+1
#define NAME_ARG MIN_ARG+1
#define NARGS NAME_ARG+1

#define NUM_YRS_WINDOW 30 /* number of years in the moving average window 
note that 30 and 31 have the same effect: the window consists of the current year and 
15 years on either side of the current year */

Boolean chk_nc(int io_rtnval);
void err_exit(char * msg);
void instructions();

static inline float * floatP(float * base, int time, short row, short col, short nrows, short ncols)
{
  float * finalP;
  
  finalP = base + (long)time*(long)nrows*(long)ncols + (long)row*(long)ncols + (long)col;
  return(finalP);
} // end of floatP()


static inline short * shortP(short * base, int time, short row, short col, short nrows, short ncols)
{
  short * finalP;
  
  finalP = base + (long)time*(long)nrows*(long)ncols + (long)row*(long)ncols + (long)col;
  return(finalP);
} // end of shortP()


int main(int argc, char * argv[])
{
  int iarg;
  char varname[STRING_LEN+1], infilename[STRING_LEN+1], eqfilename[STRING_LEN+1], outfilename[STRING_LEN+1];
  int infileid, eqfileid, invarid, eqvarid;
  size_t NumMosHist, NumYrsHist, NumRows, NumCols;
  nc_type invartype;
  int valCount;
  short * InputDataValsShortP;
  short * EQvalsShortP;
  short * SpinupValsShortP;
  float * InputDataValsFloatP;
  float * EQvalsFloatP;
  float * SpinupValsFloatP;
  int row, col, ThisMo, ThisYr;
  int outfileid, outvarid, outbanddimid, outlatdimid, outlondimid, outdimids[3];
  int outlatvarid, outlonvarid, inlatvarid, inlonvarid;
  float * lat_arrayInP;
  float * lon_arrayInP;
  char info[STRING_LEN + 1];
  int rtnval;
  nc_type vartype, attType;
  size_t attlen;
  short missing_value_short, FillValue_short;
  float missing_value_float, FillValue_float;
  int half_window_as_int = NUM_YRS_WINDOW/2;
  int AvStartYr, AvEndYr, NumYrsThisAv, year_inner;
  float windowSum, windowMean, ThisAnomaly, ThisSpinupVal, thisEQval;
  short thisShortVal, shortSpinupVal;
  float thisFloatVal;
  size_t outcoords[3], outcounts[3];
  float min_val;
  int nDims, nVars, nGatts, unlimDimid;
  int dimIndex;
  float scale_factor, scaled;
  
  if (argc!=NARGS) 
  {
    printf(COMMAND_NAME ": wrong number of arguments.\n");
    instructions();   
  }
  for (iarg = 0; iarg<argc; iarg++) printf("%s ", argv[iarg]);
  printf("\n");
  
  // Get the minimum value.
  if (sscanf(argv[MIN_ARG], "%f", &min_val)!=1) err_exit("error when reading <minimum value>.");

  // Open the historical data netCDF file.
  strncpy(varname, argv[NAME_ARG], STRING_LEN); varname[STRING_LEN] = 0;
  strncpy(infilename, argv[INPUT_PATH_ARG], STRING_LEN); infilename[STRING_LEN] = 0;
  strncat(infilename, varname, sizeof(infilename) - strlen(infilename) - 1); infilename[STRING_LEN] = 0;
  strncat(infilename, ".nc", sizeof(infilename) - strlen(infilename) - 1); infilename[STRING_LEN] = 0;
  printf("infilename = %s (file of historical data)\n", infilename);
  rtnval = nc_open(infilename, NC_NOWRITE, &infileid); 
  if (!chk_nc(rtnval)) 
  {
    err_exit("error while opening input data file.");
  }

  rtnval = nc_inq(infileid, &nDims, &nVars, &nGatts, &unlimDimid);
  assert(chk_nc(rtnval));

  for (dimIndex = 0; dimIndex < nDims; ++dimIndex)
    { char dimName[NC_MAX_NAME+1];
      size_t dimSize;

      rtnval = nc_inq_dim(infileid, dimIndex, dimName, &dimSize); assert(chk_nc(rtnval));
      if (strcmp(dimName, "time")==0 ||
	  strcmp(dimName, "month")==0 ||
	  strcmp(dimName, "band")==0 ||
	  strcmp(dimName, "year")==0)  
	{ NumMosHist = dimSize;
          NumYrsHist = NumMosHist/12;
          assert(NumYrsHist*12==NumMosHist);
          printf("NumMosHist = %ld\n", NumMosHist);
	}
      else if (strcmp(dimName, "lat")==0)
        { NumRows = dimSize;
          rtnval = nc_inq_varid(infileid, "lat", &inlatvarid); assert(chk_nc(rtnval));
          rtnval = nc_inq_vartype(infileid, inlatvarid, &vartype); assert(chk_nc(rtnval));
          assert(vartype==NC_FLOAT);
          printf("NumRows = %ld\n", NumRows);
        }
      else if (strcmp(dimName, "lon")==0)
        { NumCols = dimSize;
          rtnval = nc_inq_varid(infileid, "lon", &inlonvarid); assert(chk_nc(rtnval));
          rtnval = nc_inq_vartype(infileid, inlonvarid, &vartype); assert(chk_nc(rtnval));
          assert(vartype==NC_FLOAT);
          printf("NumCols = %ld\n", NumCols);
        }
    }

  rtnval = nc_inq_varid(infileid, varname, &invarid); assert(chk_nc(rtnval));
  rtnval = nc_inq_vartype(infileid, invarid, &invartype); assert(chk_nc(rtnval));
  assert(invartype==NC_SHORT || invartype==NC_FLOAT);
  rtnval = nc_inq_att(infileid, invarid, "missing_value", &attType, &attlen); assert(chk_nc(rtnval));
  assert(attType==invartype && attlen==1);
  rtnval = nc_inq_att(infileid, invarid, "_FillValue", &attType, &attlen); assert(chk_nc(rtnval));
  assert(attType==invartype && attlen==1);
  switch (invartype)
  {
    case NC_SHORT:
      rtnval = nc_get_att_short(infileid, invarid, "missing_value", &missing_value_short); assert(chk_nc(rtnval));
      missing_value_float = (float)missing_value_short;
      rtnval = nc_get_att_short(infileid, invarid, "_FillValue", &FillValue_short); assert(chk_nc(rtnval));
      FillValue_float = (float)FillValue_short;
      break;
    case NC_FLOAT:
      rtnval = nc_get_att_float(infileid, invarid, "missing_value", &missing_value_float); assert(chk_nc(rtnval));
      rtnval = nc_get_att_float(infileid, invarid, "_FillValue", &FillValue_float); assert(chk_nc(rtnval));
      break;
    default: assert(0); 
      break;
  }
  assert(missing_value_float==-9999.);
  assert(FillValue_float==-9999.);
  assert(missing_value_float==FillValue_float);

  // Open the EQ data netCDF file.
  strncpy(eqfilename, argv[EQ_PATH_ARG], STRING_LEN); eqfilename[STRING_LEN] = 0;
  strncat(eqfilename, varname, sizeof(eqfilename) - strlen(eqfilename) - 1); eqfilename[STRING_LEN] = 0;
  strncat(eqfilename, ".nc", sizeof(eqfilename) - strlen(eqfilename) - 1); eqfilename[STRING_LEN] = 0;
  printf("eqfilename = %s\n", eqfilename);
  rtnval = nc_open(eqfilename, NC_NOWRITE, &eqfileid); 
  if (!chk_nc(rtnval)) 
  {
    err_exit("error while opening EQ data file.");
  }
  rtnval = nc_inq_varid(eqfileid, varname, &eqvarid); assert(chk_nc(rtnval));

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
  rtnval = nc_copy_att(infileid, NC_GLOBAL, "row_offset", outfileid, NC_GLOBAL); assert(chk_nc(rtnval));
  rtnval = nc_copy_att(infileid, NC_GLOBAL, "col_offset", outfileid, NC_GLOBAL); assert(chk_nc(rtnval));
  rtnval = nc_copy_att(infileid, NC_GLOBAL, "Conventions", outfileid, NC_GLOBAL); 
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

  // Now create the climate variable.
  outdimids[0] = outbanddimid;
  outdimids[1] = outlatdimid;
  outdimids[2] = outlondimid;
  switch (invartype)
  {
    case NC_SHORT:
      rtnval = nc_def_var(outfileid, varname, NC_SHORT, 3, outdimids, &outvarid); assert(chk_nc(rtnval)); 
      break;
    case NC_FLOAT:
      rtnval = nc_def_var(outfileid, varname, NC_FLOAT, 3, outdimids, &outvarid); assert(chk_nc(rtnval)); 
      break;
    default:
      assert(0);
      break;
  }    
  nc_copy_att(infileid, invarid, "missing_value", outfileid, outvarid); assert(chk_nc(rtnval));
  nc_copy_att(infileid, invarid, "_FillValue", outfileid, outvarid); assert(chk_nc(rtnval));
  nc_copy_att(infileid, invarid, "units", outfileid, outvarid); assert(chk_nc(rtnval));
  rtnval = nc_get_att_float(infileid, invarid, "scale_factor", &scale_factor);
  if (rtnval!=NC_NOERR)
  {
    rtnval = nc_get_att_float(infileid, invarid, "scaled", &scaled);
    if (rtnval==NC_NOERR)
    {
      assert(scaled!=0.);
      scale_factor = 1./scaled;
    }
    else scale_factor = 1.0;
  }
  nc_put_att_float(outfileid, outvarid, "scale_factor", NC_FLOAT, 1, &scale_factor);
  nc_copy_att(infileid, invarid, "valid_min", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "valid_max", outfileid, outvarid);
  nc_copy_att(infileid, invarid, "long_name", outfileid, outvarid);
  rtnval = nc_enddef(outfileid); assert(chk_nc(rtnval));
   
  lat_arrayInP = (float *)malloc(sizeof(float)*NumRows); assert(lat_arrayInP!=NULL);
  rtnval = nc_get_var_float(infileid, inlatvarid, lat_arrayInP); assert(chk_nc(rtnval));
  rtnval = nc_put_var_float(outfileid, outlatvarid, lat_arrayInP); assert(chk_nc(rtnval));
  
  lon_arrayInP = (float *)malloc(sizeof(float)*NumCols); assert(lon_arrayInP!=NULL);
  rtnval = nc_get_var_float(infileid, inlonvarid, lon_arrayInP); assert(chk_nc(rtnval));
  rtnval = nc_put_var_float(outfileid, outlonvarid, lon_arrayInP); assert(chk_nc(rtnval));

  printf("Getting values from netCDF file of EQ data...\n");
  valCount = NumCols*NumRows*12;
  switch (invartype)
  {
    case NC_SHORT:
      EQvalsShortP = (short *)malloc(sizeof(short)*valCount); assert(EQvalsShortP!=NULL);
      rtnval = nc_get_var_short(eqfileid, eqvarid, EQvalsShortP); assert(chk_nc(rtnval));
      break;
    case NC_FLOAT:
      EQvalsFloatP = (float *)malloc(sizeof(float)*valCount); assert(EQvalsFloatP!=NULL);
      rtnval = nc_get_var_float(eqfileid, eqvarid, EQvalsFloatP); assert(chk_nc(rtnval));
      break;
    default: assert(0);
      break;
  }

  printf("Getting values from netCDF file of historical data to make spinup dataset...\n");
  valCount = NumCols*NumRows*NumMosHist;
  switch (invartype)
  {
    case NC_SHORT:
      InputDataValsShortP = (short *)malloc(sizeof(short)*valCount); assert(InputDataValsShortP!=NULL);
      rtnval = nc_get_var_short(infileid, invarid, InputDataValsShortP); assert(chk_nc(rtnval));
      SpinupValsShortP = (short *)malloc(sizeof(short)*valCount); assert(SpinupValsShortP!=NULL);
      break;
    case NC_FLOAT:
      InputDataValsFloatP = (float *)malloc(sizeof(float)*valCount); assert(InputDataValsFloatP!=NULL);
      rtnval = nc_get_var_float(infileid, invarid, InputDataValsFloatP); assert(chk_nc(rtnval));
      SpinupValsFloatP = (float *)malloc(sizeof(float)*valCount); assert(SpinupValsFloatP!=NULL);
      break;
    default: assert(0);
      break;
  }

  printf("Starting main loop now...\n");
  for (row = 0; row<NumRows; row++) 
  {
    printf("Starting row %d now...\n", row);
    for (col = 0; col<NumCols; col++)
    { Boolean missing_value_flag;
      missing_value_flag = FALSE;
      switch (invartype)
      {
        case NC_SHORT:
          thisShortVal = *shortP(InputDataValsShortP, 0, row, col, NumRows, NumCols);
          missing_value_flag = thisShortVal==missing_value_short;
          if (missing_value_flag) for (ThisMo = 0; ThisMo<NumMosHist; ThisMo++) 
              *shortP(SpinupValsShortP, ThisMo, row, col, NumRows, NumCols) = missing_value_short;
          break;
        case NC_FLOAT:
          thisFloatVal = *floatP(InputDataValsFloatP, 0, row, col, NumRows, NumCols);
          missing_value_flag = thisFloatVal==missing_value_float;
          if (missing_value_flag) for (ThisMo = 0; ThisMo<NumMosHist; ThisMo++) 
              *floatP(SpinupValsFloatP, ThisMo, row, col, NumRows, NumCols) = missing_value_float;
          break;
        default: assert(0);
          break;
      }
      if (!missing_value_flag)
      {
        for (ThisMo = 0; ThisMo<12; ThisMo++) for (ThisYr = 0; ThisYr<NumYrsHist; ThisYr++)
        {
          switch (invartype)
          {
            case NC_SHORT:
              thisShortVal = *shortP(InputDataValsShortP, ThisYr*12 + ThisMo, row, col, NumRows, NumCols);
              assert(thisShortVal!=missing_value_short);
              thisEQval = (float)(*shortP(EQvalsShortP, ThisMo, row, col, NumRows, NumCols));
              break;
            case NC_FLOAT:
              thisFloatVal = *floatP(InputDataValsFloatP, ThisYr*12 + ThisMo, row, col, NumRows, NumCols);
              assert(thisFloatVal!=missing_value_float);
              thisEQval = (float)(*floatP(EQvalsFloatP, ThisMo, row, col, NumRows, NumCols));
              break;
            default: assert(0);
              break;
          }              
          if (!(thisEQval*scale_factor>=min_val)) printf("*** makeSpinup: ThisMo, row, col, thisEQval, min_val = %d, %d, %d, %f, %f\n",
              ThisMo, row, col, thisEQval, min_val);
          assert(thisEQval*scale_factor>=min_val);
          AvStartYr = ThisYr - half_window_as_int; if (AvStartYr<0) AvStartYr = 0;
          AvEndYr = ThisYr + half_window_as_int; if (AvEndYr>=NumYrsHist) AvEndYr = NumYrsHist - 1;
          NumYrsThisAv = AvEndYr - AvStartYr + 1;
          windowSum = 0.;
          for (year_inner = 0; year_inner<NumYrsThisAv; year_inner++)
          { 
            switch (invartype)
            {
              case NC_SHORT:
                windowSum += (float)(*shortP(InputDataValsShortP, (year_inner + AvStartYr)*12 + ThisMo, row, col, NumRows, NumCols));
                break;
              case NC_FLOAT:
                windowSum += (float)(*floatP(InputDataValsFloatP, (year_inner + AvStartYr)*12 + ThisMo, row, col, NumRows, NumCols));
                break;
              default: assert(0); break;
            }
          }
          windowMean = windowSum/NumYrsThisAv;
          switch (invartype)
          {
            case NC_SHORT: ThisAnomaly = (float)(thisShortVal) - windowMean; break;
            case NC_FLOAT: ThisAnomaly = thisFloatVal - windowMean; break;
            default: assert(0); break;
          }
          if (ThisAnomaly<0.0 && min_val>=0.0) 
          {
            switch (invartype)
            {
              case NC_SHORT:
                assert(thisShortVal>=min_val/scale_factor && windowMean>min_val/scale_factor);
                ThisAnomaly = ((float)thisShortVal - min_val/scale_factor)/(windowMean - min_val/scale_factor);
                break;
              case NC_FLOAT:
                assert(thisFloatVal>=min_val/scale_factor && windowMean>min_val/scale_factor);
                ThisAnomaly = (thisFloatVal - min_val/scale_factor)/(windowMean - min_val/scale_factor);
                break;
              default: assert(0); break;
            }
            assert(0.<=ThisAnomaly && ThisAnomaly<=1.0);
            ThisSpinupVal = (thisEQval - min_val/scale_factor)*ThisAnomaly + min_val/scale_factor;
            assert(min_val/scale_factor<=ThisSpinupVal && ThisSpinupVal<=thisEQval);
          }
          else 
          {
            ThisSpinupVal = thisEQval + ThisAnomaly;
            if (!(ThisSpinupVal*scale_factor>=min_val)) printf("*** makeSpinup: ThisMo, row, col, ThisSpinupVal, min_val, thisEQval, ThisAnomaly = %d, %d, %d, %f, %f, %f, %f\n",
              ThisMo, row, col, ThisSpinupVal, min_val, thisEQval, ThisAnomaly);
            assert(ThisSpinupVal*scale_factor>=min_val);
          }
          
          // Now write the newly calculated spinup value into the buffer.
          switch (invartype)
          {
            case NC_SHORT:
              shortSpinupVal = ThisSpinupVal<0. ? (short)(ThisSpinupVal - 0.5) : (short)(ThisSpinupVal + 0.5);
              // shortSpinupVal = (short)ThisSpinupVal; // for truncating instead of rounding
              *shortP(SpinupValsShortP, ThisYr*12 + ThisMo, row, col, NumRows, NumCols) = shortSpinupVal;
              break;
            case NC_FLOAT:
              *floatP(SpinupValsFloatP, ThisYr*12 + ThisMo, row, col, NumRows, NumCols) = ThisSpinupVal;
              // if (ThisYr==1 && ThisSpinupVal>0.) printf("*** makeSpinup: ThisYr==1, ThisSpinupVal = %f\n", ThisSpinupVal);
              break;
            default: assert(0); break;
          }
        } // end of for ThisMo for ThisYr  
      }
    } // end of for col
  } // end of for row 

  // Now write the spinup data array to the netCDF file.
  outcoords[1] = outcoords[2] = 0;
  outcounts[0] = 1; outcounts[1] = NumRows; outcounts[2] = NumCols;
  for (ThisMo = 0; ThisMo<NumMosHist; ThisMo++)
  {
    printf("writing out layer for month %d now...\n", ThisMo);
    outcoords[0] = ThisMo; 
    switch (invartype)
    {
      case NC_SHORT:
        rtnval = nc_put_vara_short(outfileid, outvarid, outcoords, outcounts, 
            shortP(SpinupValsShortP, ThisMo, 0, 0, NumRows, NumCols)); assert(chk_nc(rtnval));
        break;
      case NC_FLOAT:
        rtnval = nc_put_vara_float(outfileid, outvarid, outcoords, outcounts, 
            floatP(SpinupValsFloatP, ThisMo, 0, 0, NumRows, NumCols)); assert(chk_nc(rtnval));
        break;
      default: assert(0); break;
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
      ": Invoke as: " COMMAND_NAME " <note> <historical data path> <EQ data path> "
          "<spinup data path> <minimum value> <variable name> \n\n"
	    "  Example for tmin: " COMMAND_NAME " 'spinup data file for Oregon based on the 2009 edition of the PRISM US800m historical climate dataset, "
      "made by David Conklin on 7/14/10 using the " COMMAND_NAME ".c program.' Input/US800m/HistoricalOregon/ Input/US800m/EQoregon/ ./ "
      "-100.0 tmin \n\n" 
	    "  Example for ppt: " COMMAND_NAME " 'spinup data file for Oregon based on the 2009 edition of the PRISM US800m historical climate dataset, "
      "made by David Conklin on 7/14/10 using the " COMMAND_NAME ".c program.' Input/US800m/HistoricalOregon/ Input/US800m/EQoregon/ ./ "
      "0.0 ppt \n\n"); 
  exit(-1);
} // end of instructions()

