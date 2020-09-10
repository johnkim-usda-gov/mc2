/*
 *  makeNCfile.c
 *  MC1
 *
 */
#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <assert.h>

#include "/opt/local/include/netcdf.h"


#define COMMAND_NAME "makeNCfile"
#include "gridspecs.h" 
#include "makeNCfile.h"


int main(int argc, char * argv[])
{
  char time_name[STRING_LEN + 1];
  Boolean timeless_flag;
  char grid_name[STRING_LEN + 1];
  char units[STRING_LEN + 1];
  double latitude0, longitude0, cell_spacing;
  int nrows_grid, ncols_grid;
  char prefix[STRING_LEN + 1];
  char suffix[STRING_LEN + 1];
  char infilename[STRING_LEN + 1];
  char varname[STRING_LEN + 1];
  char outfilename[STRING_LEN + 1];
  int fileid, timeid, rowid, colid, varid;
  short nc_short;
  double nc_double;
  int timedim, rowdim, coldim;
  int dimids[4];
  // char time_str[10];
  FILE * infile;
  int irow, jcol;
  // nclong val;
  int row_current, col_current, grid_row, grid_col;
  int io_rtnval; 
  float fval;
  // short short_val;
  long coords[3];	/*netCDF coords*/
  // float nc_time;
  int year, month;
  int nchars;
  float divisor, scale_factor;
  float valid_min, valid_max;
  int row_offset, col_offset, nrows_input, ncols_input;
  int nrows_output, ncols_output, cols_read_so_far;
  float latitude, longitude;
  char info[STRING_LEN + 1];
  int iarg;
  short * vals; // Dynamically allocated for short vals[nrows_output][ncols_output]
  int itime;
  size_t start[3];
  size_t count[3];
  // ptrdiff_t stride[3];
  fpos_t * file_pos; // current position in each input file; dynamically allocated for fpos_t file_pos[(last_year - first_year + 1)*12]
  int first_year, last_year, nmonths;
  int hdr_len;
  float filesize_as_float;
  Boolean huge_file_flag;
  Boolean month_flag, twoD_flag;
    
  if (argc!=NARGS) instructions();  
  
  // Get the time frame and initialize the time variables.
  strncpy(time_name, argv[TIME_ARG], STRING_LEN - 1); time_name[STRING_LEN] = 0;
  getTime(time_name, &first_year, &last_year, &timeless_flag, &month_flag, &twoD_flag);  
  
  nmonths = timeless_flag ? 1 : 12;
  
  // Dynamic allocation for fpos_t file_pos[(last_year - first_year + 1)*12]
  file_pos = (fpos_t *)malloc((last_year - first_year + 1)*12*sizeof(fpos_t));
  if (file_pos==NULL) err_exit("Failed to allocate memory for file_pos[(last_year - first_year + 1)*12].");

  // Get the grid name and initialize the grid specs variables.
  strncpy(grid_name, argv[GRID_ARG], STRING_LEN - 1); grid_name[STRING_LEN] = 0;
  getGrid(grid_name, &latitude0, &longitude0, &cell_spacing, &nrows_grid, &ncols_grid);

  // Get the variable name and construct the name of the netCDF file.
  strncpy(varname, argv[VAR_ARG], STRING_LEN - 1); varname[STRING_LEN] = 0;
  strncpy(outfilename, varname, STRING_LEN - 3); outfilename[STRING_LEN - 3] = 0;
  strcat(outfilename, ".nc"); outfilename[STRING_LEN] = 0;
  
  // Get the prefix and suffix for the text file name.
  strncpy(prefix, argv[PREFIX_ARG], STRING_LEN - 1); prefix[STRING_LEN] = 0;
  if (strcmp(prefix, "NO_PREFIX")==0) prefix[0] = 0;
  strncpy(suffix, argv[SUFFIX_ARG], STRING_LEN - 1); suffix[STRING_LEN] = 0;
  if (strcmp(suffix, "NO_SUFFIX")==0) suffix[0] = 0;

  // Get hdr_len, row_offset, col_offset, nrows_input, and ncols_input.
  if (sscanf(argv[HDR_LEN_ARG], "%d", &hdr_len)!=1) err_exit("error when reading <hdr_len>.");
  if (hdr_len<0) err_exit("hdr_len<0.");
  if (sscanf(argv[ROW_OFFSET_ARG], "%d", &row_offset)!=1) err_exit("error when reading <row_offset>.");
  if (sscanf(argv[COL_OFFSET_ARG], "%d", &col_offset)!=1) err_exit("error when reading <col_offset>.");
  if (sscanf(argv[NROWS_INPUT_ARG], "%d", &nrows_input)!=1) err_exit("error when reading <nrows_input>.");
  if (sscanf(argv[NCOLS_INPUT_ARG], "%d", &ncols_input)!=1) err_exit("error when reading <ncols_input>.");  
  if (row_offset>=nrows_grid || col_offset>=ncols_grid) err_exit("Error in row_offset or col_offset.");
  if (nrows_input<1) err_exit("Error in row_offset or nrows_input.");
  if (ncols_input<1) err_exit("Error in col_offset or ncols_input.");
  nrows_output = nrows_input + MIN(row_offset, 0);
  nrows_output = MIN(nrows_output, nrows_grid);
  ncols_output = ncols_input + MIN(col_offset, 0);
  ncols_output = MIN(ncols_output, ncols_grid);
  if (nrows_output<1 || ncols_output<1) err_exit("Input grid and output grid are disjoint.");
  
  // Dynamic allocation for short vals[nrows_output][ncols_output];
  vals = (short *)malloc(nrows_output*ncols_output*sizeof(short));
  if (vals==NULL) err_exit("Failed to allocate memory for vals[nrows_output][ncols_output].");

  // Get divisor.
  if (sscanf(argv[DIVISOR_ARG], "%f", &divisor)!=1) err_exit("error when reading <divisor>.");
  if (divisor==0) err_exit("divisor is 0.");
  
  // Get valid_min and valid_max.
  if (sscanf(argv[VALID_MIN_ARG], "%f", &valid_min)!=1) err_exit("error when reading <valid_min>.");
  if (sscanf(argv[VALID_MAX_ARG], "%f", &valid_max)!=1) err_exit("error when reading <valid_max>.");
   
  // Get units string.
  strncpy(units, argv[UNITS_ARG], STRING_LEN - 1); units[STRING_LEN] = 0;
  
  // Get scale_factor.
  if (sscanf(argv[SCALE_FACTOR_ARG], "%f", &scale_factor)!=1) err_exit("error when reading <scale_factor>.");
  if (scale_factor==0) err_exit("scale_factor is 0.");
    
  printf("Creating full-grid output netCDF file: %s\n", outfilename);	
  filesize_as_float = ((float)(nrows_output*ncols_output))*(last_year - first_year + 1)*12*sizeof(short);
  huge_file_flag = filesize_as_float > 1.e9;
  io_rtnval = nc_create(outfilename, huge_file_flag ? NC_CLOBBER | NC_64BIT_OFFSET : NC_CLOBBER, &fileid); 
  assert(chk_nc(io_rtnval));

  /*  Global attributes for row and column start offsets */
  nc_short = MAX(row_offset, 0);
  ncattput(fileid, NC_GLOBAL, "row_offset", NC_SHORT, 1, &(nc_short) );
  nc_short = MAX(col_offset, 0);
  ncattput(fileid, NC_GLOBAL, "col_offset", NC_SHORT, 1, &(nc_short) );

  if (!timeless_flag)
  { /*  Time dimension */
    io_rtnval = nc_def_dim(fileid, "month", huge_file_flag ? NC_UNLIMITED : (last_year - first_year + 1)*12, 
        &timedim); assert(chk_nc(io_rtnval));
    printf("\nTime dim for output file defined.\n"); 
  }

  /*  Row, col dimensions */
  rowdim = ncdimdef(fileid, "lat", nrows_output);
  coldim = ncdimdef(fileid, "lon", ncols_output);
  printf("\nRow/col dims for output file defined as %d, %d.\n", nrows_output, ncols_output); 

  /* row coordinate variable */
  dimids[0] = rowdim;
  rowid = ncvardef(fileid, "lat", NC_FLOAT, 1, dimids); 
  ncattput(fileid, rowid, "long_name", NC_CHAR, strlen("latitude"), "latitude");
  ncattput(fileid, rowid, "standard_name", NC_CHAR, strlen("latitude"), "latitude");
  ncattput(fileid, rowid, "units", NC_CHAR, strlen("degrees_north"), "degrees_north");

  /* column coordinate variable */
  dimids[0] = coldim;
  colid = ncvardef(fileid, "lon", NC_FLOAT, 1, dimids); 
  ncattput(fileid, colid, "long_name", NC_CHAR, strlen("longitude"), "longitude");
  ncattput(fileid, colid, "standard_name", NC_CHAR, strlen("longitude"), "longitude");
  ncattput(fileid, colid, "units", NC_CHAR, strlen("degrees_east"), "degrees_east");
  
  if (!timeless_flag)
  { 
    if (!huge_file_flag)
    { /* time coordinate variable */
      dimids[0] = timedim;
      timeid = ncvardef(fileid, "month", NC_DOUBLE, 1, dimids);
      ncattput(fileid, timeid, "long_name", NC_CHAR, strlen("year and month"), "year and month");
      ncattput(fileid, timeid, "standard_name", NC_CHAR, strlen("year"), "year");
      ncattput(fileid, timeid, "units", NC_CHAR, strlen("year"), "year");
    }
    
    /* target variable */
    dimids[0] = timedim;
    dimids[1] = rowdim;
    dimids[2] = coldim;
    varid = ncvardef(fileid, argv[VAR_ARG], NC_SHORT, 3, dimids);
  }
  else
  {
    /* target variable */
    dimids[0] = rowdim;
    dimids[1] = coldim;
    varid = ncvardef(fileid, argv[VAR_ARG], NC_SHORT, 2, dimids);
  }
  
  nc_short = MISSING_VALUE;
  ncattput(fileid, varid, "missing_value", NC_SHORT, 1, &nc_short);
  ncattput(fileid, varid, "_FillValue", NC_SHORT, 1, &nc_short);
  ncattput(fileid, varid, "units", NC_CHAR, strlen(argv[UNITS_ARG]), argv[UNITS_ARG]);
  // Get scale_factor.
  if (sscanf(argv[SCALE_FACTOR_ARG], "%f", &scale_factor)!=1) err_exit("error when reading <scale_factor>.");
  if (scale_factor==0) 
  {
    printf("*** makeNCfile: scale_factor is 0 in command line.  Setting scale_factor to 1 and continuing.\n");
	scale_factor = 1.;
  }
  ncattput(fileid, varid, "scale_factor", NC_FLOAT, 1, &scale_factor); 
  // Rescale valid_min and valid_max before writing them out as attributes so that they can be compared directly 
  // with the 16 bit integers in the netCDF file.
  nc_short = valid_min/scale_factor;
  ncattput(fileid, varid, "valid_min", NC_SHORT, 1, &(nc_short));
  nc_short = valid_max/scale_factor;
  ncattput(fileid, varid, "valid_max", NC_SHORT, 1, &(nc_short));
 
  info[0] = 0;
  for (iarg=0; iarg<(argc - 1); iarg++)
  {
    if (iarg==NOTE_ARG) continue;
    strncat(info, argv[iarg], STRING_LEN - strlen(info));
	strncat(info, " ", STRING_LEN - strlen(info));
  }
  
  ncattput(fileid, NC_GLOBAL, "cmdline", NC_CHAR, strlen(info), info);
  sprintf(info, "Grid name is %s, origin is %f latitude, %f longitude.  Cell spacing is %f.", 
      grid_name, latitude0, longitude0, cell_spacing);
  ncattput(fileid, NC_GLOBAL, "grid_specs", NC_CHAR, strlen(info), info);
  if (!timeless_flag)
  {
    sprintf(info, "Time series is monthly for %d thru %d." , first_year, last_year);
    ncattput(fileid, NC_GLOBAL, "time_series", NC_CHAR, strlen(info), info);
  }
  ncattput(fileid, NC_GLOBAL, "note", NC_CHAR, strlen(argv[NOTE_ARG]), argv[NOTE_ARG]);
  
  ncendef(fileid);
  ncsync(fileid);

  /* Write out the coordinate variable arrays. */
  if (!timeless_flag && !huge_file_flag)
  {
    itime = 0;
    for (year = first_year; year<=last_year; year++)
    {
      for (month = 1; month<=12; month++)
	  {
	    coords[0] = itime;
	    nc_double = year + ((float)month - 0.5)/12.;
	    ncvarput1(fileid, timeid, coords, (void *)(&nc_double));
	    itime++;
	  } // end of month loop
    } // end of year loop
  } // if (!timeless_flag && !huge_file_flag) ...
  	
  grid_row = MAX(0, row_offset);
  for (irow = 0; irow<nrows_output; irow++) 
  {
	coords[0] = irow;
	latitude = latitude0 - ((float)grid_row + 0.5)*cell_spacing;
	ncvarput1( fileid, rowid, coords, (void *)(&latitude));
	grid_row++;
  } // end of row loop
  
  grid_col = MAX(0, col_offset);
  for (jcol = 0; jcol<ncols_output; jcol++)
  {
	coords[0] = jcol;
	longitude = longitude0 + ((float)grid_col + 0.5)*cell_spacing;
	ncvarput1( fileid, colid, coords, (void *)(&longitude));
	grid_col++;
  } // end of column loop

  // Loop thru all the months
  itime = 0;
  year = first_year;
  int bad_header_count = 0;
  while (year<=last_year)
  {
    month = 1;
	while (month<=nmonths)
	{
	  // Make the text file name for this year and month
	  if (timeless_flag) nchars = sprintf(infilename, "%s%s%s", prefix, varname, suffix);
	  else nchars = sprintf(infilename, "%s_%d%02d.txt", argv[VAR_ARG], year, month);
	  if (nchars>(sizeof(infilename) - 1)) err_exit("input file name is too long.");
	  if (nchars<=0) err_exit("error while constructing input file name.");
  
    // Try to open the text file
	  if ((infile = fopen(infilename,"r"))==NULL) 
    {
      // Text file didn't open, so try to convert the BIL format data for this month
      // to a text file.
      char cmd_text[256];
      cmd_text[0] = 0;
      nchars = sprintf(cmd_text, "gdal_translate -of AAIGrid -ot Float32 -co DECIMAL_PRECISION=3 %d/%s%04d%02d%s %s_%04d%02d.txt\n", 
          year, prefix, year, month, suffix, argv[VAR_ARG], year, month);
      if (nchars>(sizeof(cmd_text) - 1)) err_exit("cmd_text is too long.");
      system(cmd_text);
    
      // Try again to open the text file
      if ((infile = fopen(infilename,"r"))==NULL) 
      {
        printf("Error while opening input file %s\n", infilename);
        err_exit("...aborting now.");
      }
    }
    printf("Opened file %s.\n", infilename);
	  
    Boolean bad_header_flag = FALSE;
	  if (itime==0) read_header(hdr_len, infile);
    else bad_header_flag = !check_and_skip_header(&hdr_len, infile); // Read past the header.
    if (bad_header_flag) bad_header_count++;
	  if (row_offset<0)
	  {
	    skip_rows(-row_offset, infile);
		row_current = 0;
	  }
	  else row_current = row_offset;
	  for (irow = 0; irow<nrows_output; irow++)
	  {
		if (col_offset<0)
		{
		  skip_cols(-col_offset, infile);
		  cols_read_so_far = -col_offset;
		  col_current = 0;
		}
		else 
		{
		  cols_read_so_far = 0;
		  col_current = col_offset;
		}
		for (jcol = 0; jcol<ncols_output; jcol++)
		{
		  io_rtnval = fscanf(infile, "%f,", &fval);
		  if (io_rtnval<=0) 
		  {
		    printf("*** makeNCfile: jcol = %d, ncols_output = %d; irow = %d, nrows_output = %d\n",
			    jcol, ncols_output, irow, nrows_output);
		    err_exit("Error while reading columns.");
		  }
	  
      // In the code below, read "*(vals + irow*ncols_output + jcol)" as "vals[irow][jcol]".
      fval = fval!=MISSING_VALUE ? (float)fval/divisor : MISSING_VALUE;
			assert(-32768.<=fval && fval<=32767.); // Confirm that the data is within the range of a short (16 bit) integer.
      if (fval>=0.) *(vals + irow*ncols_output + jcol) = (short)(fval + 0.5);
      else *(vals + irow*ncols_output + jcol) = (short)(fval - 0.5);
		  
		  col_current++;
		} // end of column loop	
		cols_read_so_far += ncols_output;
		if (cols_read_so_far<ncols_input) skip_rows(1, infile);
		
		row_current++;
      } // end of row loop

	  // Write out one complete spatial grid.
	  if (timeless_flag)
	  {
	    start[0] = start[1] = start[2] = 0;
	    count[0] = nrows_output; count[1] = ncols_output;
		printf("makeNCfile: count[] = %lu, %lu, %lu\n", count[0], count[1], count[2]);
	  }
	  else
	  {
	    start[0] = itime; start[1] = 0; start[2] = 0;
	    count[0] = 1; count[1] = nrows_output; count[2] = ncols_output;
	  }
	  io_rtnval = nc_put_vara_short(fileid, varid, start, count, vals); 
	  if (io_rtnval!=NC_NOERR) { printf("io_rtnval = %d\n", io_rtnval); assert(0); }
	  
	  ncsync(fileid);
	  fclose(infile);  
	  month++;
	  itime++;
	} // end of month loop
	printf("Year %d completed.\n", year);
	year++; 
  } // end of year loop
  if (bad_header_count>0) printf("*** makeNCfile: bad_header_count = %d\n", bad_header_count);
     
  ncclose(fileid);
  free(vals);
  free(file_pos);
  
  return(0);
} // end of main()


void instructions()
{
  printf("\n*** " COMMAND_NAME
      ": Invoke as '" COMMAND_NAME 
	  " <time_name> <grid_name> <variable name> <input_file_prefix> <input_file_suffix> "
	  "<hdr_len> <row_offset> <col_offset> <nrows_input> <ncols_input> <divisor> <valid_min> <valid_max> <units> <scale_factor> <note>'.\n"
	  "  <time_name> is either 'TIMELESS', '1895_2006', '1895_2008', or '1895_1896'.\n"
	  "  <grid_name> is 'USeast4km', 'US800m', 'CA800m', 'PNW800m', 'Yellowstone800m', or 'CA12km'.\n"
      "  <input_file_suffix> should be specified as 'NO_SUFFIX' if there is no suffix.\n"
	  "  <hrd_len> is the number of lines of header in the input file which must be skipped to get to the data itself.\n"
	  "  row_offset is the offset of the first row of the input file relative to the first row of the grid.\n"
	  "  Similarly, col_offset.\n"
	  "  <divisor> is used to scale down the input data before it is written out to the netCDF file.\n"
	  "  <valid min> and <valid max> are written into the netCDF file as attributes of the variable, for use by ncview, etc.\n"
	  "  <valid min> and <valid max> should be specified in the units given in the next argument.\n"
	  "  <units> are the units string for the units of the variable after multiplying by <scale_factor>\n"
	  "  For conventions on units strings see http://www.unidata.ucar.edu/software/udunits/udunits-1/udunits.txt\n"
	  "  <scale_factor> is a number which is written into the netCDF file as the 'scale_factor' attribute of the variable.\n"
	  "  The scale_factor should be chosen so that when you multiply the 16 bit integer form of the data used in the netCDF\n"
	  "  file by the scale_factor, the resulting value is in the units specified in the previous argument.\n"
	  "  <note> is a string describing the data.\n"
	  "  Example: makeNCfile TIMELESS CA800M mask sclass .txt 6 457 524 102 121 1. 0 11 sclass 1. 'sclass data from Bill Kuhn on 6/18/08,"
	  " converted to .nc file on 6/24/08 by Dave Conklin'\n");
  printf("NetCDF version: %s\n", nc_inq_libvers()); 
  exit(-1);
} // end of instructions()



