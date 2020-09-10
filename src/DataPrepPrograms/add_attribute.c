/*
 *  add_attribute.c
 *  DataPrepPrograms
 *
 *  Created by Dave Conklin on 9/22/10.
 *  Copyright 2010 Conservation Biology Institute. All rights reserved.
 *
 */
// A command line utility to add an attribute to the header of an existing netCDF file

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <assert.h>
#include <math.h>
#include "/opt/local/include/netcdf.h"

#define Boolean int
#define TRUE 1
#define FALSE 0
#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)

#define STRING_LEN 200

#define COMMAND_NAME "add_attribute"
#define FILE_ARG 1 
#define VAR_ARG FILE_ARG+1
#define ATT_NAME_ARG VAR_ARG+1
#define ATT_TYPE_ARG ATT_NAME_ARG+1
#define ATT_VALUE_ARG ATT_TYPE_ARG+1
#define NARGS ATT_VALUE_ARG+1

Boolean chk_nc(int io_rtnval);
void err_exit(char * msg);


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
      ": Invoke as '" COMMAND_NAME " <netCDF file path and name> <variable name | NC_GLOBAL> <attribute name> <attribute type> <attribute value>\n"
	    "  Example: " COMMAND_NAME " soils.nc soil_data missing_value NC_SHORT -9999\n");
  exit(-1);
} // end of instructions()


int main(int argc, char * argv[])
{
  char tgtfilename[STRING_LEN+1];
  char varname[STRING_LEN+1];
  char attname[STRING_LEN+1];
  char att_type[STRING_LEN+1];
  char attvalue_str[STRING_LEN+1];
  int i, tgtfileid, varid, attid, rtnval;
  Boolean global_att_flag;
  signed char attvalue_byte;
  short attvalue_short;
  int attvalue_int;
  float attvalue_float;
  
  if (argc!=NARGS) instructions();   
  for (i=0; i<argc; i++) printf("%s ", argv[i]);
  printf("\n");

  // Open this input file.
  strncpy(tgtfilename, argv[FILE_ARG], STRING_LEN - 1); tgtfilename[STRING_LEN] = 0;
  rtnval = nc_open(tgtfilename, NC_WRITE, &tgtfileid); 
  if (!chk_nc(rtnval)) err_exit("error while opening target data file.");

  // Find the variable, if any.
  global_att_flag = strcmp(argv[VAR_ARG], "NC_GLOBAL")==0;
  if (global_att_flag) varid = NC_GLOBAL;
  else 
  {
    strncpy(varname, argv[VAR_ARG], STRING_LEN - 1); varname[STRING_LEN] = 0;
    rtnval = nc_inq_varid(tgtfileid, varname, &varid); 
    if (!chk_nc(rtnval)) err_exit("couldn't find the variable.");
  }

  // Find out if the attribute already exists.
  strncpy(attname, argv[ATT_NAME_ARG], STRING_LEN - 1); attname[STRING_LEN] = 0;
  rtnval = nc_inq_attid(tgtfileid, varid, attname, &attid);
  if (rtnval!=NC_ENOTATT) err_exit("attribute already exists.");
  
  // Create the attribute.
  strncpy(att_type, argv[ATT_TYPE_ARG], STRING_LEN - 1); att_type[STRING_LEN] = 0;
  strncpy(attvalue_str, argv[ATT_VALUE_ARG], STRING_LEN - 1); attvalue_str[STRING_LEN] = 0;
  rtnval = nc_redef(tgtfileid); assert(chk_nc(rtnval));
  if (strcmp(att_type, "NC_BYTE")==0)
  {
    attvalue_int = atoi(attvalue_str);
    if (attvalue_int<-128 || attvalue_int>127) err_exit("attribute value is out of range for NC_BYTE type.");
    attvalue_byte = (char)attvalue_int;
    rtnval = nc_put_att_schar(tgtfileid, varid, attname, NC_BYTE, 1, &attvalue_byte); assert(chk_nc(rtnval));
  }
  else if (strcmp(att_type, "NC_SHORT")==0)
  {
    attvalue_short = atoi(attvalue_str);
    rtnval = nc_put_att_short(tgtfileid, varid, attname, NC_SHORT, 1, &attvalue_short); assert(chk_nc(rtnval));
  }
  else if (strcmp(att_type, "NC_INT")==0)
  {
    attvalue_int = atoi(attvalue_str);
    rtnval = nc_put_att_int(tgtfileid, varid, attname, NC_INT, 1, &attvalue_int); assert(chk_nc(rtnval));
  }
  else if (strcmp(att_type, "NC_FLOAT")==0)
  {
    attvalue_float = atof(attvalue_str);
    rtnval = nc_put_att_float(tgtfileid, varid, attname, NC_FLOAT, 1, &attvalue_float); assert(chk_nc(rtnval));
  }
  else if (strcmp(att_type, "NC_TEXT")==0)
  {
    rtnval = nc_put_att_text(tgtfileid, varid, attname, strlen(attvalue_str), attvalue_str); assert(chk_nc(rtnval));
  }
  else err_exit("couldn't recognize attribute type; it was not NC_BYTE, NC_SHORT, NC_INT, NC_FLOAT, or NC_TEXT.");
  rtnval = nc_enddef(tgtfileid); assert(chk_nc(rtnval));
  rtnval = nc_close(tgtfileid); assert(chk_nc(rtnval));
  
  printf("\n" COMMAND_NAME " terminated successfully.\n\n");
} // end of main()

