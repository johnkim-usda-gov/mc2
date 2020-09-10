/*
 *  delete_attribute.c
 *  DataPrepPrograms
 *
 *  Created by Dave Conklin on 9/23/10.
 *  Copyright 2010 Conservation Biology Institute. All rights reserved.
 *
 */

// A command line utility to delete an attribute from the header of an existing netCDF file

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

#define COMMAND_NAME "delete_attribute"
#define FILE_ARG 1 
#define VAR_ARG FILE_ARG+1
#define ATT_NAME_ARG VAR_ARG+1
#define NARGS ATT_NAME_ARG+1

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
      ": Invoke as '" COMMAND_NAME " <netCDF file path and name> <variable name | NC_GLOBAL> <attribute name>\n"
	    "  Example: " COMMAND_NAME " soils.nc soil_data units\n");
  exit(-1);
} // end of instructions()


int main(int argc, char * argv[])
{
  char tgtfilename[STRING_LEN+1];
  char varname[STRING_LEN+1];
  char attname[STRING_LEN+1];
  int i, tgtfileid, varid, attid, rtnval;
  Boolean global_att_flag;
  
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

  // Find out if the attribute exists.
  strncpy(attname, argv[ATT_NAME_ARG], STRING_LEN - 1); attname[STRING_LEN] = 0;
  rtnval = nc_inq_attid(tgtfileid, varid, attname, &attid);
  if (rtnval==NC_ENOTATT) err_exit("couldn't find the attribute.");
  
  // Delete the attribute.
  rtnval = nc_redef(tgtfileid); assert(chk_nc(rtnval));
  rtnval = nc_del_att(tgtfileid, varid, attname); assert(chk_nc(rtnval));  
  rtnval = nc_enddef(tgtfileid); assert(chk_nc(rtnval));
  
  rtnval = nc_close(tgtfileid); assert(chk_nc(rtnval));  
  printf("\n" COMMAND_NAME " terminated successfully.\n\n");
} // end of main()


