/*
 *  mergeNCdata.c
 */
// A command line utility to merge data extracted from netCDF files. 

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <assert.h>
#include <math.h>
#include "/opt/local/include/netcdf.h"

#include "../gridspecs.h"

#define Boolean int
#define TRUE 1
#define FALSE 0
#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)

#define STRING_LEN 800

#define COMMAND_NAME "mergeNCdata"
#define CMDFILE_ARG 1
#define NARGS CMDFILE_ARG+1

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
      ": Invoke as '" COMMAND_NAME " <command file path and name>\n"
	  "  Example: " COMMAND_NAME " mergelist.txt\n"
      "  where the first line of mergelist.txt has \n"
      "      <grid name> <output file path and name> <output variable name> <first row> <last row> <first col> <last col> <time_name> <nc_type>\n"
	    "      <grid_name> is 'VEMAP', 'US12km', 'NA8km', 'US800m', 'BearSoil800m', 'ATtest800m', 'CA800m', 'PNW800m', 'Yellowstone800m', 'CA12km', or 'CA10km'.\n"
	    "      <time_name> is either 'TIMELESS', '1900_1999', '1901_2000', '1971_2000', '1895_2006', '1895_2007', '1895_2008', 'MONTHS_1895_2008',\n"
      "            '1998_2007', '2001_2007', '2007_2099', '1895_1896' or '12_MONTHS'\n"
      "      <nc_type> is 'NC_FLOAT', 'NC_SHORT', 'NC_INT', or 'NC_BYTE'.\n"
      "  and succeeding lines have \n"
      "      <input file path and name> <row offset> <col offset> "
      "<input variable name> <index to 1st band> <number of bands>\n"
			"  N.B. index to 1st band starts at 1, not zero!\n"
      "  set <input variable name> to 'multiply', 'additive_normalize', 'difference', 'ratio', 'fire_index',\n"
      "  'neg->0', 'neg->FillValue', '100-last2', '%%change', 'sum' for various special cases.\n");
  exit(-1);
} // end of instructions()


int main(int argc, char * argv[])
{
  char buffer[STRING_LEN+1];
  char cmdfilename[STRING_LEN+1], gridname[STRING_LEN+1];
  FILE * cmdfile;
  Boolean firstline_flag, secondline_flag, inband_flag;
  int i, rtnval;
  short shortval;
  int intval;
  unsigned char byteval;
  float floatval;
  char time_name[STRING_LEN+1];
  char type_name[STRING_LEN+1];
  int first_year, last_year;
  Boolean timeless_flag, month_flag, twoD_flag;
  Boolean stop_flag;
  char info[STRING_LEN + 1];
  int iarg;
	size_t attlen;
  int cell_count;
  size_t latloncoords[1];
  char latlondimname[NC_MAX_NAME+1];
  int ndx_pos;
	        
  // variables pertaining to the output file
  char outfilename[STRING_LEN+1], outvarname[STRING_LEN+1];
  int outfileid, outlatdimid, outlondimid, outbanddimid;
  int outlatvarid, outlonvarid;
  int outvarid;
  nc_type outvartype;
  int outdimids[3];
  size_t outcoords[3];
  int outband, outrow, outcol;
  int nrows_out, ncols_out;
  int outyear, outyearvarid, outmonth;
	int calendar_year;
  float calendar_year_float;
  float outlat, outlon;
  float unity = 1.0;
  float scale_factor;
  
  // variables pertaining to the input files 
  char infilename[STRING_LEN+1], invarname[STRING_LEN+1];
  int infileid, inlatdimid, inlondimid, inbanddimid;
  int in_lat_ndx_pos, in_lon_ndx_pos, in_band_ndx_pos;
  size_t inbanddimlen, nrows_in, ncols_in;
  int invarid, invarndims;
  nc_type invartype;
  int indimids[3];
  size_t incoords[3];
  int inband, inrow, incol;
  int row_offset, col_offset, infirstband, innumbands, inlastband;
  Boolean fire_index_flag, multiply_flag, additive_normalize_flag, difference_flag, ratio_flag, 
      neg2zero_flag, neg2FillValue_flag, p100minuslast2_flag, pctchange_flag, sum_flag;
  int termA_firstband, termB_firstband, termB_band;
  char note_str[] = "notei"; // Will hold "note0", "note1", ...
  char note_str_out[] = "notei_from_input_file";
  int inlatvarid, inlonvarid;
  float inlat, inlon, scale_factor_in, inverse_scale_factor_in; 

  // variables in grid units
  int first_row, last_row, first_col, last_col;
  int gridrow, gridcol;
  int nrows_grid, ncols_grid;  
  
  // miscellaneous variables
  double latitude0, longitude0, cell_spacing;
  float latitude, longitude;
  int cmdfile_line_num, inputfile_num;
  char cmdfile_line_num_str[STRING_LEN+1];
  char inputfile_num_str[STRING_LEN+1];
  float missing_value_float = NC_FILL_FLOAT;
  float fill_value_float = NC_FILL_FLOAT;
  short missing_value_short = NC_FILL_SHORT;
  short fill_value_short = NC_FILL_SHORT;
  int missing_value_int = NC_FILL_INT;
  int fill_value_int = NC_FILL_INT;
  Boolean missing_value_attribute_flag, fill_value_attribute_flag;
  nc_type missing_value_type, fill_value_type, att_type;
  char * rtnptr;
  Boolean continueFlag;
  
  if (argc!=NARGS) instructions();   
  for (i=0; i<argc; i++) printf("%s ", argv[i]);
  printf("\n");

  strncpy(cmdfilename, argv[CMDFILE_ARG], STRING_LEN - 1); cmdfilename[STRING_LEN] = 0;
  if ((cmdfile = fopen(cmdfilename,"r"))==NULL) err_exit("error while opening command file.");
  
  continueFlag = TRUE;    
  while (continueFlag)
  {
    continueFlag = FALSE;
    firstline_flag = TRUE;
    secondline_flag = FALSE;
    outband = 0;
    stop_flag = FALSE;
    cell_count = 0;
    cmdfile_line_num = 0;
    cmdfile_line_num_str[0] = inputfile_num_str[0] = 0;
    while (!feof(cmdfile) && !stop_flag && !continueFlag)
    {
      rtnptr = fgets(buffer, STRING_LEN, cmdfile);
      if (rtnptr==NULL) 
      {
        stop_flag = TRUE;
        continue;
      }
      if (buffer[0]==0 || buffer[0]=='#') continue; // Ignore empty lines and lines that start with '#'.
      
      cmdfile_line_num++;
      printf("%s", buffer);
      assert(cmdfile_line_num>0 && cmdfile_line_num<=99);
      sprintf(cmdfile_line_num_str, "cmdfile_line_%02d", cmdfile_line_num);
      assert(strlen(cmdfile_line_num_str)<=STRING_LEN);
      
      if (firstline_flag)
      {
        // Scan first line of command file to get
        // the grid name, the path and name of the output file,
        // and the name of the output variable.
        memset((char *)gridname, '\0', sizeof(gridname));	
        memset((char *)outfilename, '\0', sizeof(outfilename));	
        memset((char *)outvarname, '\0', sizeof(outvarname));	
        memset((char *)time_name, '\0', sizeof(time_name));	
        memset((char *)type_name, '\0', sizeof(type_name));	
        rtnval = sscanf(buffer, "%s%s%s%d%d%d%d%s%s", gridname, outfilename, outvarname, &first_row, &last_row, &first_col, &last_col, time_name, type_name);
        if (rtnval!=9 || strlen(gridname)==0 || strlen(outfilename)==0 || strlen(outvarname)==0
            || first_row<0 || last_row<first_row || first_col<0 || last_col<first_col || strlen(time_name)==0 || strlen(type_name)==0) 
          {
            printf("rtnval = %d\n", rtnval);
            printf("strlen(gridname) = %d\n", (int)strlen(gridname));
            printf("strlen(outfilename) = %d\n", (int)strlen(outfilename));
            printf("strlen(outvarname) = %d\n", (int)strlen(outvarname));
            printf("first_row, last_row, first_col, last_col = %d, %d, %d, %d\n", first_row, last_row, first_col, last_col);
            printf("strlen(time_name) = %d\n", (int)strlen(time_name));
            printf("strlen(type_name) = %d\n", (int)strlen(type_name));
            err_exit("Error in first line of command file.");
          }
        firstline_flag = FALSE;
        getGrid(gridname, &latitude0, &longitude0, &cell_spacing, &nrows_grid, &ncols_grid);
        if (nrows_grid<=last_row || ncols_grid<=last_col) err_exit("Rows or columns in command line are inconsistent with the grid.");
        nrows_out = last_row - first_row + 1; assert(1<=nrows_out && nrows_out<=nrows_grid);
        ncols_out = last_col - first_col + 1; assert(1<=ncols_out && ncols_out<=ncols_grid);
        getTime(time_name, &first_year, &last_year, &timeless_flag, &month_flag, &twoD_flag); assert(!twoD_flag || timeless_flag);

        // Create the output file.  
        printf("Creating output netCDF file: %s\n", outfilename);	
        rtnval = nc_create(outfilename, NC_NOCLOBBER, &outfileid); chk_nc(rtnval);
        if (rtnval==NC_EEXIST)
        {
          printf("*** Output netCDF file \"%s\" already exists.\n", outfilename);	
          err_exit("Aborting now.");
        }
        assert(chk_nc(rtnval));

        if (!twoD_flag) { rtnval = nc_def_dim(outfileid, timeless_flag ? "band" : "year", NC_UNLIMITED, &outbanddimid); assert(chk_nc(rtnval)); }
        rtnval = nc_def_dim(outfileid, "lat", nrows_out, &outlatdimid); assert(chk_nc(rtnval));
        rtnval = nc_def_dim(outfileid, "lon", ncols_out, &outlondimid); assert(chk_nc(rtnval));
        rtnval = nc_put_att_int(outfileid, NC_GLOBAL, "row_offset", NC_INT, 1, &first_row); assert(chk_nc(rtnval));
        rtnval = nc_put_att_int(outfileid, NC_GLOBAL, "col_offset", NC_INT, 1, &first_col); assert(chk_nc(rtnval));
        info[0] = 0;
        for (iarg=0; iarg<(argc); iarg++)
        {
          strncat(info, argv[iarg], STRING_LEN - strlen(info));
          strncat(info, " ", STRING_LEN - strlen(info));
        }  
        rtnval = nc_put_att_text(outfileid, NC_GLOBAL, "cmdline", strlen(info), info); assert(chk_nc(rtnval));
        rtnval = nc_put_att_text(outfileid, NC_GLOBAL, cmdfile_line_num_str, strlen(buffer) - 1, buffer); assert(chk_nc(rtnval));
        sprintf(info, "Grid name is %s, origin is %f latitude, %f longitude.  Cell spacing is %f.", 
              gridname, latitude0, longitude0, cell_spacing);
        rtnval = nc_put_att(outfileid, NC_GLOBAL, "grid_specs", NC_CHAR, strlen(info), info);
          
        // Create the coordinate dimension variables.
        outdimids[0] = outlatdimid;
        if (latitude0!=NOT_LAT_LON)
        {
          rtnval = nc_def_var(outfileid, "lat", NC_FLOAT, 1, outdimids, &outlatvarid); assert(chk_nc(rtnval));
          rtnval = nc_put_att_text(outfileid, outlatvarid, "long_name", strlen("latitude"), "latitude"); assert(chk_nc(rtnval));
          rtnval = nc_put_att_text(outfileid, outlatvarid, "standard_name", strlen("latitude"), "latitude"); assert(chk_nc(rtnval));
          rtnval = nc_put_att_text(outfileid, outlatvarid, "units", strlen("degrees_north"), "degrees_north"); assert(chk_nc(rtnval));
        }
        outdimids[0] = outlondimid;
        if (longitude0!=NOT_LAT_LON)
        {
          rtnval = nc_def_var(outfileid, "lon", NC_FLOAT, 1, outdimids, &outlonvarid); assert(chk_nc(rtnval));
          rtnval = nc_put_att_text(outfileid, outlonvarid, "long_name", strlen("longitude"), "longitude"); assert(chk_nc(rtnval));
          rtnval = nc_put_att_text(outfileid, outlonvarid, "standard_name", strlen("longitude"), "longitude"); assert(chk_nc(rtnval));
          rtnval = nc_put_att_text(outfileid, outlonvarid, "units", strlen("degrees_east"), "degrees_east"); assert(chk_nc(rtnval));
        }
        if (!timeless_flag)
        {
          outdimids[0] = outbanddimid;
          rtnval = nc_def_var(outfileid, "year", NC_FLOAT, 1, outdimids, &outyearvarid); assert(chk_nc(rtnval));
          rtnval = nc_put_att_text(outfileid, outyearvarid, "long_name", strlen("calendar year"), "calendar year"); assert(chk_nc(rtnval));
          rtnval = nc_put_att_text(outfileid, outyearvarid, "standard_name", strlen("year"), "year"); assert(chk_nc(rtnval));
          rtnval = nc_put_att_text(outfileid, outyearvarid, "units", strlen("year"), "year"); assert(chk_nc(rtnval));
        }

        // Define the output variable.
        printf("output variable type = %s\n", type_name);
        if (strcmp(type_name, "NC_FLOAT")==0) outvartype = NC_FLOAT;
        else if (strcmp(type_name, "NC_SHORT")==0) outvartype = NC_SHORT;
        else if (strcmp(type_name, "NC_INT")==0) outvartype = NC_INT;
        else if (strcmp(type_name, "NC_BYTE")==0) outvartype = NC_BYTE;
        else err_exit("Unable to interpret output variable type.");
        if (twoD_flag)
        {
          outdimids[0] = outlatdimid;
          outdimids[1] = outlondimid;
          rtnval = nc_def_var(outfileid, outvarname, outvartype, 2, outdimids, &outvarid); assert(chk_nc(rtnval));
        }
        else 
        {
          outdimids[0] = outbanddimid;
          outdimids[1] = outlatdimid;
          outdimids[2] = outlondimid;
          rtnval = nc_def_var(outfileid, outvarname, outvartype, 3, outdimids, &outvarid); assert(chk_nc(rtnval));
        }
        outband = 0;
        scale_factor = 1.0;
        if (outvartype==NC_FLOAT) 
        { // Put in a unity scale factor to satisfy older versions of MC1.
          rtnval = nc_put_att_float(outfileid, outvarid, "scale_factor", NC_FLOAT, 1, &unity); assert(chk_nc(rtnval));
        }
    
        // Write out the values of the coordinate dimension variables.      
        rtnval = nc_enddef(outfileid); assert(chk_nc(rtnval));      
        if (latitude0!=NOT_LAT_LON)
        { // Write out lat coordinate variable.
          printf("Writing out latitude coordinate variable now...\n");
          for (outrow = 0; outrow<nrows_out; outrow++) 
          {
          outcoords[0] = outrow;
          gridrow = first_row + outrow;
          latitude = latitude0 - ((float)gridrow + 0.5)*cell_spacing;
          rtnval = nc_put_var1_float(outfileid, outlatvarid, outcoords, &latitude); assert(chk_nc(rtnval));
          } // end of row loop
        }
        if (longitude0!=NOT_LAT_LON)
        { // Write out lon coordinate variable.
          printf("Writing out longitude coordinate variable now...\n");
          for (outcol = 0; outcol<ncols_out; outcol++)
          {
            outcoords[0] = outcol;
            gridcol = first_col + outcol;
            longitude = longitude0 + ((float)gridcol + 0.5)*cell_spacing;
            rtnval = nc_put_var1_float(outfileid, outlonvarid, outcoords, &longitude); assert(chk_nc(rtnval));
          } // end of column loop
        }

        rtnval = nc_sync(outfileid); assert(chk_nc(rtnval));
        
        firstline_flag = FALSE;
        secondline_flag = TRUE;
        continue;
      } // end of logic to process first line of command file

      // Logic for second and subsequent lines of command file
      // Rescan the input buffer to get
      //  <input file path and name> <row offset> <col offset> <input variable name> <index to 1st band> <number of bands>
      rtnval = sscanf(buffer, "%s%d%d%s%d%d", infilename, &row_offset, &col_offset, invarname, &infirstband, &innumbands);
      if (rtnval<1) { stop_flag = TRUE; continue; }
      if (!secondline_flag && rtnval==1 && strcmp(infilename, "continue")==0) { continueFlag = TRUE; continue; }
      if (rtnval!= 6 || nrows_grid<=row_offset || ncols_grid<=col_offset)
      {
        if (rtnval!=6) printf("rtnval = %d.  Should be 6.\n", rtnval);
        else 
        {
          if (nrows_grid<=row_offset) printf("nrows_grid, row_offset = %d, %d. nrows_grid should be > row_offset.\n", nrows_grid, row_offset);
          if (ncols_grid<=col_offset) printf("ncols_grid, col_offset = %d, %d. ncols_grid should be > col_offset.\n", ncols_grid, col_offset);
        }
        err_exit("Input data specs are inconsistent.");
      }
      inband_flag = innumbands>0;
      if (inband_flag && infirstband<1) err_exit("Input band specs are inconsistent."); 
      if (inband_flag) inlastband = infirstband + innumbands - 1;
      else infirstband = inlastband = 1;   
      
      // The next line is for merging data spatially instead of temporally.
      // if (inband_flag && infirstband==1 && innumbands==outband) outband = 0;

      rtnval = nc_redef(outfileid); assert(chk_nc(rtnval));      

      // Copy the command file line to the output file.
      rtnval = nc_put_att_text(outfileid, NC_GLOBAL, cmdfile_line_num_str, strlen(buffer) - 1, buffer); assert(chk_nc(rtnval));
      
      // Open this input file.
      rtnval = nc_open(infilename, NC_NOWRITE, &infileid); 
      if (!chk_nc(rtnval)) err_exit("error while opening input data file.");
      
      // Copy the input file's command line to the output file, if possible.
      attlen = 0;
      rtnval = nc_inq_attlen(infileid, NC_GLOBAL, "cmdline", &attlen); 
      assert(attlen<sizeof(info)); 
      rtnval = nc_get_att_text(infileid, NC_GLOBAL, "cmdline", info); 
      info[attlen] = 0; // Make sure the string is null-terminated.
      if (rtnval==NC_NOERR) 
      {
        // Input file names start on the second line of the command file.
        inputfile_num = cmdfile_line_num - 1;
        assert(inputfile_num>0 && inputfile_num<=99);
        sprintf(inputfile_num_str, "inputfile_%02d_cmdline", inputfile_num);
        assert(strlen(inputfile_num_str)<=STRING_LEN);
        rtnval = nc_put_att_text(outfileid, NC_GLOBAL, inputfile_num_str, strlen(info), info); assert(chk_nc(rtnval));
      }
      
      rtnval = nc_enddef(outfileid); assert(chk_nc(rtnval));      
      
      
      if (inband_flag)
      {
        rtnval = nc_inq_dimid(infileid, "band", &inbanddimid); 
        if (rtnval!=NC_NOERR) rtnval = nc_inq_dimid(infileid, "year", &inbanddimid); 
        if (rtnval!=NC_NOERR) rtnval = nc_inq_dimid(infileid, "month", &inbanddimid); 
        if (rtnval!=NC_NOERR) rtnval = nc_inq_dimid(infileid, "time", &inbanddimid); 
        if (rtnval!=NC_NOERR) err_exit("Didn't recognize the band or time coordinate variable name");
        rtnval = nc_inq_dimlen(infileid, inbanddimid, &inbanddimlen); assert(chk_nc(rtnval));
        if (inlastband>inbanddimlen) err_exit("Not enough bands in the input variable.");
      }  
      else 
      {
        inbanddimid = -1;
        inbanddimlen = 0;
      } // end of if inband_flag else ...
      
      if (!timeless_flag && (secondline_flag || inband_flag))
      { // Write out year coordinate variable.
        printf("Writing out year coordinate variable now...\n");
        if (!month_flag)
        { // yearly timestep
          calendar_year = first_year - 1 + infirstband;
          for (outyear = outband; outyear<(outband + innumbands); outyear++)
          {
            outcoords[0] = outyear;
            calendar_year_float = calendar_year;
            rtnval = nc_put_var1_float(outfileid, outyearvarid, outcoords, &calendar_year_float); assert(chk_nc(rtnval));
            calendar_year++;
          }
        }
        else 
        { // monthly timestep
          int num_months;
          
          num_months = inband_flag ? innumbands : (last_year - first_year + 1)*12;
          for (outmonth = outband; outmonth<(outband + num_months); outmonth++)
          {
            outcoords[0] = outmonth;
            calendar_year_float = first_year + outmonth/12.0 + 1/24.0;
            rtnval = nc_put_var1_float(outfileid, outyearvarid, outcoords, &calendar_year_float); assert(chk_nc(rtnval));
          } // end of for outmonth
        } // end of if !month_flag else
      }	// end of if !timeless_flag... 		
      
      rtnval = nc_inq_dimid(infileid, "lat", &inlatdimid); 
      if (rtnval!=NC_NOERR) rtnval = nc_inq_dimid(infileid, "y", &inlatdimid); 
      if (rtnval!=NC_NOERR) rtnval = nc_inq_dimid(infileid, "row", &inlatdimid);
      assert(chk_nc(rtnval));
      rtnval = nc_inq_dimlen(infileid, inlatdimid, &nrows_in); assert(chk_nc(rtnval));
      rtnval = nc_inq_dimname(infileid, inlatdimid, latlondimname); assert(chk_nc(rtnval));
      rtnval = nc_inq_varid(infileid, latlondimname, &inlatvarid); if (rtnval!=NC_NOERR) inlatvarid = -1;
      
      rtnval = nc_inq_dimid(infileid, "lon", &inlondimid); 
      if (rtnval!=NC_NOERR) rtnval = nc_inq_dimid(infileid, "x", &inlondimid); 
      if (rtnval!=NC_NOERR) rtnval = nc_inq_dimid(infileid, "col", &inlondimid); 
      assert(chk_nc(rtnval));
      rtnval = nc_inq_dimlen(infileid, inlondimid, &ncols_in); assert(chk_nc(rtnval));
      rtnval = nc_inq_dimname(infileid, inlondimid, latlondimname); assert(chk_nc(rtnval));
      rtnval = nc_inq_varid(infileid, latlondimname, &inlonvarid); if (rtnval!=NC_NOERR) inlonvarid = -1;

      rtnval = nc_inq_varid(infileid, invarname, &invarid); 
      fire_index_flag = multiply_flag = additive_normalize_flag = difference_flag = ratio_flag = 
          neg2zero_flag = neg2FillValue_flag = p100minuslast2_flag = pctchange_flag = sum_flag = FALSE;
      if (rtnval==NC_NOERR)
      {
        rtnval = nc_inq_vartype(infileid, invarid, &invartype); assert(chk_nc(rtnval));
        rtnval = nc_inq_varndims(infileid, invarid, &invarndims); assert(chk_nc(rtnval));
        assert(2<=invarndims && invarndims<=3);
        if (inband_flag && invarndims!=3) 
            err_exit("Input variable has no time or band dimension.  Use 0 for starting band number and 0 for number of bands.");
        if (!inband_flag && invarndims!=2) 
            err_exit("Input variable has a time or band dimension.  Specify a non-zero starting band number and at least 1 for the number of bands.");
        rtnval = nc_inq_vardimid(infileid, invarid, indimids); assert(chk_nc(rtnval));
        in_lat_ndx_pos = in_lon_ndx_pos = in_band_ndx_pos = -1; 
        for (ndx_pos = 0; ndx_pos<invarndims; ndx_pos++) 
        {
          if (indimids[ndx_pos]==inlatdimid) in_lat_ndx_pos = ndx_pos;
          if (indimids[ndx_pos]==inlondimid) in_lon_ndx_pos = ndx_pos;
          if (inband_flag && indimids[ndx_pos]==inbanddimid) in_band_ndx_pos = ndx_pos;
        }
        printf("inlatdimid, inlondimid, inbanddimid = %d, %d, %d\n", inlatdimid, inlondimid, inbanddimid);
        printf("in_lat_ndx_pos, in_lon_ndx_pos, in_band_ndx_pos = %d, %d, %d\n", in_lat_ndx_pos, in_lon_ndx_pos, in_band_ndx_pos);
        assert(in_lat_ndx_pos>=0 && in_lon_ndx_pos>=0);
        if (inband_flag) assert(in_band_ndx_pos>=0);
       
        rtnval = nc_get_att_float(infileid, invarid, "scale_factor", &scale_factor_in);
        if (rtnval==NC_NOERR) 
        {
          printf("scale_factor for %s from file %s = %f\n", invarname, infilename, scale_factor_in);
        }
        else 
        {
          rtnval = nc_inq_att(infileid, invarid, "scaled", &att_type, &attlen);
          if (rtnval==NC_NOERR) 
          {
            if (att_type==NC_FLOAT) 
            {
              rtnval = nc_get_att_float(infileid, invarid, "scaled", &inverse_scale_factor_in); assert(chk_nc(rtnval));
            }
            else if (att_type==NC_CHAR)
            {
              char * scale_factor_str;
              int count;
              
              scale_factor_str = malloc(attlen*nctypelen(att_type) + 2);
              assert(scale_factor_str!=NULL);
              rtnval = nc_get_att_text(infileid, invarid, "scaled", scale_factor_str); assert(chk_nc(rtnval));
              strcat(scale_factor_str,"\0");
              count = sscanf(scale_factor_str, "%f", &inverse_scale_factor_in);
              if (count!=1) err_exit("Error when interpreting scale_factor.");
              assert(inverse_scale_factor_in!=0.);
              free(scale_factor_str);
            }
            else err_exit("'scaled' attribute is not of type NC_FLOAT or NC_CHAR");
            scale_factor_in = 1./inverse_scale_factor_in;
          }        
          else scale_factor_in = 1.0;
        }

        // If this is the first input file, copy attributes from first input file to the output file.
        if (secondline_flag)
        {
          rtnval = nc_redef(outfileid); assert(chk_nc(rtnval));      
          
          // Global attributes: grid_specs, time_series, note
          attlen = 0;
          rtnval = nc_inq_attlen(infileid, NC_GLOBAL, "grid_specs", &attlen); assert(attlen<sizeof(info));
          rtnval = nc_get_att_text(infileid, NC_GLOBAL, "grid_specs", info); 
          if (rtnval==NC_NOERR) 
          {
            info[attlen] = 0; // Make sure the string is null-terminated.
            rtnval = nc_put_att_text(outfileid, NC_GLOBAL, "grid_specs_from_input_file", strlen(info), info); assert(chk_nc(rtnval));
          }
          rtnval = nc_inq_attlen(infileid, NC_GLOBAL, "time_series", &attlen); assert(attlen<sizeof(info));
          rtnval = nc_get_att_text(infileid, NC_GLOBAL, "time_series", info); 
          if (rtnval==NC_NOERR) 
          {
            info[attlen] = 0; // Make sure the string is null-terminated.
            rtnval = nc_put_att_text(outfileid, NC_GLOBAL, "time_series_from_input_file", strlen(info), info); assert(chk_nc(rtnval));
          }
          rtnval = nc_inq_attlen(infileid, NC_GLOBAL, "note", &attlen); assert(attlen<sizeof(info));
          rtnval = nc_get_att_text(infileid, NC_GLOBAL, "note", info); 
          if (rtnval==NC_NOERR) 
          {
            info[attlen] = 0; // Make sure the string is null-terminated.
            rtnval = nc_put_att_text(outfileid, NC_GLOBAL, "note_from_input_file", strlen(info), info); assert(chk_nc(rtnval));
          }
          for (i = 0; i<=9; i++)
          {
            note_str[4] = (char)(0x30 + i);
            note_str_out[0] = 0;
            strcpy(note_str_out, note_str);
            strcat(note_str_out, "_from_input_file");
            rtnval = nc_inq_attlen(infileid, NC_GLOBAL, note_str, &attlen); 
            if (rtnval==NC_NOERR && attlen>0)
            {
              assert(attlen<sizeof(info));
              rtnval = nc_get_att_text(infileid, NC_GLOBAL, note_str, info); assert(chk_nc(rtnval));
              info[attlen] = 0; // Make sure the string is null-terminated.
              rtnval = nc_put_att_text(outfileid, NC_GLOBAL, note_str_out, strlen(info), info); assert(chk_nc(rtnval));
            }
          }
          nc_copy_att(infileid, NC_GLOBAL, "esri_pe_string", outfileid, outvarid);
        
          if (outvartype!=NC_FLOAT) 
          {
//            if (invartype!=outvartype) err_exit("outvartype!=NC_FLOAT && invartype!=outvartype");            
            if (invartype!=outvartype) printf("*** mergeNCdata: outvartype!=NC_FLOAT && invartype!=outvartype\n");            
            scale_factor = scale_factor_in;
            rtnval = nc_put_att_float(outfileid, outvarid, "scale_factor", NC_FLOAT, 1, &scale_factor); assert(chk_nc(rtnval));
          }
          if (invartype==outvartype) nc_copy_att(infileid, invarid, "_FillValue", outfileid, outvarid);
        
          // Variable attributes
          nc_copy_att(infileid, invarid, "missing_value", outfileid, outvarid);
          nc_copy_att(infileid, invarid, "units", outfileid, outvarid);
          nc_copy_att(infileid, invarid, "valid_min", outfileid, outvarid);
          nc_copy_att(infileid, invarid, "valid_max", outfileid, outvarid);
          nc_copy_att(infileid, invarid, "long_name", outfileid, outvarid);
          // nc_copy_att(infileid, invarid, "scaled", outfileid, outvarid);
          nc_copy_att(infileid, invarid, "date_created", outfileid, outvarid);
          nc_copy_att(infileid, invarid, "source", outfileid, outvarid);
          nc_copy_att(infileid, invarid, "esri_pe_string", outfileid, outvarid);
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
        
          rtnval = nc_enddef(outfileid); assert(chk_nc(rtnval)); 

          missing_value_attribute_flag = fill_value_attribute_flag = FALSE;
          secondline_flag = FALSE;
        }

        rtnval = nc_inq_atttype(infileid, invarid, "missing_value", &missing_value_type);
        if (rtnval==NC_NOERR && missing_value_type==NC_FLOAT)
        {
          missing_value_attribute_flag = TRUE;
          rtnval = nc_get_att_float(infileid, invarid, "missing_value", &missing_value_float); assert(chk_nc(rtnval));
        }
        else if (rtnval==NC_NOERR && missing_value_type==NC_SHORT)
        {
          missing_value_attribute_flag = TRUE;
          rtnval = nc_get_att_short(infileid, invarid, "missing_value", &missing_value_short); assert(chk_nc(rtnval));
        }
        else if (rtnval==NC_NOERR && missing_value_type==NC_INT)
        {
          missing_value_attribute_flag = TRUE;
          rtnval = nc_get_att_int(infileid, invarid, "missing_value", &missing_value_int); assert(chk_nc(rtnval));
        }

        rtnval = nc_inq_atttype(infileid, invarid, "_FillValue", &fill_value_type);
        if (rtnval==NC_NOERR && fill_value_type==NC_FLOAT)
        {
          fill_value_attribute_flag = TRUE;
          rtnval = nc_get_att_float(infileid, invarid, "_FillValue", &fill_value_float); assert(chk_nc(rtnval));
        }
        else if (rtnval==NC_NOERR && fill_value_type==NC_SHORT)
        {
          fill_value_attribute_flag = TRUE;
          rtnval = nc_get_att_short(infileid, invarid, "_FillValue", &fill_value_short); assert(chk_nc(rtnval));
          if (outvartype==NC_FLOAT) fill_value_float = fill_value_short;
        }
        else if (rtnval==NC_NOERR && fill_value_type==NC_INT)
        {
          fill_value_attribute_flag = TRUE;
          rtnval = nc_get_att_int(infileid, invarid, "_FillValue", &fill_value_int); assert(chk_nc(rtnval));
          if (outvartype==NC_FLOAT) fill_value_float = fill_value_int;
        }
        
        if (outvartype==NC_FLOAT && fill_value_attribute_flag)
        {
          rtnval = nc_redef(outfileid); assert(chk_nc(rtnval)); 
          nc_put_att_float(outfileid, outvarid, "_FillValue", NC_FLOAT, 1, &fill_value_float);     
          rtnval = nc_enddef(outfileid); assert(chk_nc(rtnval)); 
        }
      }
      else if (rtnval==NC_ENOTVAR && strcmp(invarname, "fire_index")==0)
      { // Special case for calculating fire_index
        fire_index_flag = TRUE;
      }
      else if (rtnval==NC_ENOTVAR && strcmp(invarname, "multiply")==0)
      { // Special case for scaling the 2nd previous variable time series by the previous time-invariant variable 
        // For example scaling the time series of NPP by the cell areas
        multiply_flag = TRUE;
      }
      else if (rtnval==NC_ENOTVAR && strcmp(invarname, "additive_normalize")==0)
      { // Special case for normalizing the 2nd previous variable time series by the previous time-invariant variable 
        // For example normalizing the time series of annual temperatures by subtracting off a long-term mean
        additive_normalize_flag = TRUE;
      }
      else if (rtnval==NC_ENOTVAR && strcmp(invarname, "difference")==0)
      { // Special case for calculating the difference of the previous 2 variables, which may be time series
        difference_flag = TRUE;
      }
      else if (rtnval==NC_ENOTVAR && strcmp(invarname, "ratio")==0)
      { // Special case for calculating the ratio of the previous 2 variables, which may be time series
        ratio_flag = TRUE;
      }
      else if (rtnval==NC_ENOTVAR && strcmp(invarname, "neg->0")==0)
      { // Special case for replacing negative values with 0
        neg2zero_flag = TRUE;
      }
      else if (rtnval==NC_ENOTVAR && strcmp(invarname, "100-last2")==0)
      { // Special case for replacing negative values with 0
        p100minuslast2_flag = TRUE;
      }
      else if (rtnval==NC_ENOTVAR && strcmp(invarname, "%change")==0)
      { // Special case for replacing negative values with 0
        pctchange_flag = TRUE;
      }
      else if (rtnval==NC_ENOTVAR && strcmp(invarname, "neg->FillValue")==0)
      { // Special case for replacing negative values with the FillValue
        neg2FillValue_flag = TRUE;
      }
      else if (rtnval==NC_ENOTVAR && strcmp(invarname, "sum")==0)
      { // Special case for calculating the sum of the previous 2 variables, which may be time series
        sum_flag = TRUE;
      }
      else assert(chk_nc(rtnval));
      if (fire_index_flag || multiply_flag || additive_normalize_flag || p100minuslast2_flag || pctchange_flag) 
      { // termB is not expected to be a time series
        termB_band = outband - 1;
        termA_firstband = outband - 1 - (inband_flag ? innumbands : 1);
        assert(infirstband==1);
      }
      else if (difference_flag || ratio_flag || sum_flag)
      { // termB may be a time series
        termB_band = outband - 1;
        termA_firstband = outband - (inband_flag ? 2*innumbands : 2);
        termB_firstband = outband - (inband_flag ? innumbands : 1);
        assert(infirstband==1);
      }
      else if (neg2zero_flag || neg2FillValue_flag) 
      { // There isn't any termA; termB may be a time series.
        termA_firstband = -1;
        termB_band = outband - 1;
        termB_firstband = outband - (inband_flag ? innumbands : 1);
      }
      else termB_band = termA_firstband = -1;

      // Read the input data and write it to the output file.    
      for (inband = infirstband - 1; inband<=(inlastband - 1); inband++)
      {
        if (inlastband>infirstband) printf(COMMAND_NAME " is starting work on band %d now.\n", inband+1);
        if (inband_flag) incoords[in_band_ndx_pos] = inband;
        if (!twoD_flag) outcoords[0] = outband;
        for (outrow = 0; outrow<nrows_out; outrow++)
        {
          gridrow = first_row + outrow;
          inrow = gridrow - row_offset;
          if (inrow<0 || nrows_in<=inrow) continue;
          incoords[in_lat_ndx_pos] = inrow;
          outcoords[twoD_flag ? 0 : 1] = outrow;
          for (outcol = 0; outcol<ncols_out; outcol++)
          {
            gridcol = first_col + outcol;
            incol = gridcol - col_offset;
            if (incol<0 || ncols_in<=incol) continue;

            if (inlatvarid>=0 && inlonvarid>=0 && inband==(infirstband - 1) && 
                !(fire_index_flag || multiply_flag || additive_normalize_flag || difference_flag || ratio_flag || sum_flag))
            {
              // Read the lat,lon of this cell from the input file and compare to the lat,lon in the output file.
              latloncoords[0] = inrow;
              rtnval = nc_get_var1_float(infileid, inlatvarid, latloncoords, &inlat); assert(chk_nc(rtnval));
              latloncoords[0] = incol;
              rtnval = nc_get_var1_float(infileid, inlonvarid, latloncoords, &inlon); assert(chk_nc(rtnval));
              if (inlon>180.) inlon -= 360.;
              latloncoords[0] = outrow;
              rtnval = nc_get_var1_float(outfileid, outlatvarid, latloncoords, &outlat); assert(chk_nc(rtnval));
              latloncoords[0] = outcol;
              rtnval = nc_get_var1_float(outfileid, outlonvarid, latloncoords, &outlon); assert(chk_nc(rtnval));
              if (
                  (fabs(inlat)<1.E+10 && fabs(inlon)<1.E+10 && inlat!=0. && inlon!=0.) 
                  && (!(fabs(inlat - outlat)<0.002) || !(fabs(inlon - outlon)<0.002)) // a tolerance of about 7 arc-seconds
                  && !(inlat==inrow && inlon==incol))
              { 
                printf("*** " COMMAND_NAME ": georegistration error!\n");
                printf("                   inlat, outlat, inlon, outlon = %f, %f, %f, %f\n", inlat, outlat, inlon, outlon);
                printf("                   inrow, incol = %d, %d; outrow, outcol = %d, %d\n", inrow, incol, outrow, outcol);
                printf("                   row_offset, col_offset = %d, %d\n", row_offset, col_offset);
                assert(0);
              } 
            }         

            incoords[in_lon_ndx_pos] = incol;
            outcoords[twoD_flag ? 1 : 2] = outcol;
            if (fire_index_flag || multiply_flag || additive_normalize_flag || difference_flag || ratio_flag
                || neg2zero_flag || neg2FillValue_flag || p100minuslast2_flag || pctchange_flag || sum_flag)
            { // Instead of reading the data from the input file, calculate it from the data already read.
              float consumed, npp, fire_index, product, termB, termA;
              short termAs, termBs;
              int termAint, termBint;
              Boolean missing_value_flag = FALSE;
              
              assert((termA_firstband>=0 && termB_band>0) || ((neg2zero_flag || neg2FillValue_flag) && termB_band>=0));
              assert(outvartype==NC_FLOAT || outvartype==NC_SHORT || outvartype==NC_INT);
              
              if (!(neg2zero_flag || neg2FillValue_flag))
              {
                // Get termA.
                outcoords[0] = termA_firstband + inband;
                switch (outvartype)
                {
                  case NC_FLOAT:
                    rtnval = nc_get_var1_float(outfileid, outvarid, outcoords, &termA); assert(chk_nc(rtnval));
                    missing_value_flag = (missing_value_attribute_flag && termA==missing_value_float) ||
                        (fill_value_attribute_flag && termA==fill_value_float);
                    break;
                  case NC_SHORT:
                    rtnval = nc_get_var1_short(outfileid, outvarid, outcoords, &termAs); assert(chk_nc(rtnval));
                    missing_value_flag = (missing_value_attribute_flag && termAs==missing_value_short) ||
                        (fill_value_attribute_flag && termAs==fill_value_short);
                    termA = termAs;
                  case NC_INT:
                    rtnval = nc_get_var1_int(outfileid, outvarid, outcoords, &termAint); assert(chk_nc(rtnval));
                    missing_value_flag = (missing_value_attribute_flag && termAint==missing_value_int) ||
                        (fill_value_attribute_flag && termAint==fill_value_int);
                    termA = termAint;
                    break;
                  default: assert(0);
                }
              }
              
              // Get termB
              if (fire_index_flag || multiply_flag || additive_normalize_flag || p100minuslast2_flag || pctchange_flag) 
                  outcoords[0] = termB_band;
              else if (difference_flag || ratio_flag || neg2zero_flag || neg2FillValue_flag || sum_flag) outcoords[0] = termB_firstband + inband;
              else assert(0);
              switch (outvartype)
              {
                case NC_FLOAT:
                  rtnval = nc_get_var1_float(outfileid, outvarid, outcoords, &termB); assert(chk_nc(rtnval));
                  missing_value_flag |= missing_value_attribute_flag && termB==missing_value_float;
                  missing_value_flag |= fill_value_attribute_flag && termB==fill_value_float;
                  break;
                case NC_SHORT:
                  rtnval = nc_get_var1_short(outfileid, outvarid, outcoords, &termBs); assert(chk_nc(rtnval));
                  missing_value_flag |= missing_value_attribute_flag && termBs==missing_value_short;
                  missing_value_flag |= fill_value_attribute_flag && termBs==fill_value_short;
                  termB = termBs;
                  break;
                case NC_INT:
                  rtnval = nc_get_var1_int(outfileid, outvarid, outcoords, &termBint); assert(chk_nc(rtnval));
                  missing_value_flag = (missing_value_attribute_flag && termBint==missing_value_int) ||
                      (fill_value_attribute_flag && termBint==fill_value_int);
                  termB = termBint;
                  break;
                default: assert(0);
              }

              // Calculate the output value.
              if (missing_value_flag && !(neg2zero_flag || neg2FillValue_flag)) 
              {
                assert(missing_value_attribute_flag || fill_value_attribute_flag);
                if (missing_value_attribute_flag) { floatval = missing_value_float; shortval = missing_value_short; }
                if (fill_value_attribute_flag) { floatval = fill_value_float; shortval = fill_value_short; }
              }
              else
              {
                if (fire_index_flag)
                {
                  npp = termB;
                  consumed = termA;
                  if (consumed<0. || consumed>1.e6) continue;
                  if (consumed<=0.0) fire_index = 0.;
                  else
                  {
                    if (fabs(npp)>1.e6) continue;
                    fire_index = npp>0. ? consumed/npp : 1.;
                  }
                  floatval = fire_index;
                }
                else if (additive_normalize_flag || difference_flag) floatval = termA - termB;
                else if (sum_flag) floatval = termA + termB;
                else if (multiply_flag)
                {
                  /*
                  // Here termB is a constant (e.g. varea) for each cell, while termA may be a time series.
                  // if (termA>1.01 || termA<0.) continue; // This is specific to termA being part_burn.
                  outcoords[0] = outband - (inband_flag ? inband - infirstband + 1 : 1);
                  printf("normalize: inband(termB) = %d, outband = %d\n", inband, outband);
                  rtnval = nc_get_var1_float(outfileid, outvarid, outcoords, &termB); 
                  if (!chk_nc(rtnval))
                  {
                    printf("outband = %d; outcoords[0:2] = %ld, %ld, %ld\n", outband, outcoords[0], outcoords[1], outcoords[2]);
                    assert(0);
                  }
                  */
                  product = termA*termB;
                  floatval = product;
                }
                else if (ratio_flag) floatval = termB!=0.0 ? termA/termB : (termA==termB ? 1. : 0.0);
                else if (neg2zero_flag) floatval = termB<0.0 ? 0.0 : termB;
                else if (neg2FillValue_flag) floatval = termB<0.0 ? fill_value_float : termB;
                else if (p100minuslast2_flag) 
                { // for making sand% from silt% and clay%
                  floatval = 100. - termA - termB;
                  if (floatval<1.) floatval = 1.;
                }
                else if (pctchange_flag) floatval = termA!=0. ? 100.*(termB - termA)/termA : 
                    (termB==0 ? 0. : (termB>0. ? 100. : -100.));
                
                shortval = (short)floatval;
                intval = (int)floatval;
                byteval = (char)floatval;
              }
              
              // Restore outcoords[0], which is the band index in the output file where the output value will be placed.
              outcoords[0] = outband;
              switch (outvartype)
              {
                case NC_FLOAT: rtnval = nc_put_var1_float(outfileid, outvarid, outcoords, &floatval); assert(chk_nc(rtnval)); break;
                case NC_SHORT: rtnval = nc_put_var1_short(outfileid, outvarid, outcoords, &shortval); assert(chk_nc(rtnval)); break;
                case NC_INT:   rtnval = nc_put_var1_int(outfileid, outvarid, outcoords, &intval); assert(chk_nc(rtnval)); break;
                case NC_BYTE:  rtnval = nc_put_var1_uchar(outfileid, outvarid, outcoords, &byteval); assert(chk_nc(rtnval)); break;
                default: assert(0); break;
              }
            }
            else 
            {
              Boolean missing_value_flag = FALSE;
              switch (invartype)
              {
                case NC_SHORT:
                  rtnval = nc_get_var1_short(infileid, invarid, incoords, &shortval); assert(chk_nc(rtnval));
                  missing_value_flag = (missing_value_attribute_flag && shortval==missing_value_short) ||
                      (fill_value_attribute_flag && shortval==fill_value_short);
                  floatval = missing_value_flag ? shortval : shortval*scale_factor_in;
                  break;
                case NC_INT:
                  rtnval = nc_get_var1_int(infileid, invarid, incoords, &intval); 
                  if (!chk_nc(rtnval))
                      assert(chk_nc(rtnval));
                  missing_value_flag = (missing_value_attribute_flag && intval==missing_value_int) ||
                      (fill_value_attribute_flag && intval==fill_value_int);
                  floatval = missing_value_flag ? intval : intval*scale_factor_in;
                  break;
                case NC_BYTE: // unsigned char
                  rtnval = nc_get_var1_uchar(infileid, invarid, incoords, &byteval); 
                  /* { size_t testincoords[3];
                    testincoords[0] = incoords[1];
                    testincoords[1] = incoords[0];
                    testincoords[2] = incoords[2];
                    rtnval = nc_get_var1_uchar(infileid, invarid, testincoords, &byteval); 
                  } */
                  if (!chk_nc(rtnval))
                      assert(0);
//                  floatval = byteval + 128.f; assert(0.<=floatval && floatval<=255.f);
                  if (byteval!=0) cell_count++;
                  floatval = byteval*scale_factor_in;
                  // if (byteval==255) continue;
                  break;
                case NC_FLOAT:
                  rtnval = nc_get_var1_float(infileid, invarid, incoords, &floatval); assert(chk_nc(rtnval));
                  missing_value_flag = (missing_value_attribute_flag && floatval==missing_value_float) ||
                      (fill_value_attribute_flag && floatval==fill_value_float);
                  floatval = missing_value_flag ? floatval : floatval*scale_factor_in;
                  break;
                default: assert(0); break;
              }
              switch (outvartype)
              {
                case NC_SHORT:
                  shortval = missing_value_flag ? floatval : floatval/scale_factor;
                  rtnval = nc_put_var1_short(outfileid, outvarid, outcoords, &shortval); assert(chk_nc(rtnval));
                  break;
                case NC_INT:
                  intval = missing_value_flag ? floatval : floatval/scale_factor;
                  rtnval = nc_put_var1_int(outfileid, outvarid, outcoords, &intval); assert(chk_nc(rtnval));
                  break;
                case NC_BYTE:
                  byteval = floatval/scale_factor;
                  rtnval = nc_put_var1_uchar(outfileid, outvarid, outcoords, &byteval); assert(chk_nc(rtnval));
                  break;
                case NC_FLOAT:
                  floatval = missing_value_flag ? floatval : floatval/scale_factor;
                  rtnval = nc_put_var1_float(outfileid, outvarid, outcoords, &floatval); assert(chk_nc(rtnval));
                  break;
                default: assert(0); break;
              }
            }
          } // end of outcol loop
        } // end of outrow loop
        outband++;
        rtnval = nc_sync(outfileid); assert(chk_nc(rtnval));
      } // end of inband loop
    
      rtnval = nc_close(infileid); assert(chk_nc(rtnval));
    
    } // of loop to read lines from the command file associated with the creation of a single output file

    if (cmdfile_line_num<2) err_exit("premature exit (cmdfile_line_num<2)");
    
    rtnval = nc_close(outfileid); assert(chk_nc(rtnval));
    if (cell_count>0) printf("for output file %s, cell_count = %d (count of non-zero byte values read from input file or files)\n", 
        outfilename, cell_count);
  
  } // end of while (continueFlag)
  
  printf("\n" COMMAND_NAME " has exited normally.\n");
  
} // end of main()
          
