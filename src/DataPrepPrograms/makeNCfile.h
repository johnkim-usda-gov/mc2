/*
 *  makeNCfile.h
 *  MC1
 *
 *
 */

/* #include "gridspecs.h" */

#ifndef MAKENCFILE_H
#define	MAKENCFILE_H

#define Boolean int
#define TRUE 1
#define FALSE 0
#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)

#define STRING_LEN 200
#define MAX_HDR_LINES 12
#define MAX_HDR_LINE_LEN STRING_LEN

// A command line utility to create a netCDF file with one variable.
// Sets the file up to contain a monthly climate variable for the Yosemite 800m grid for 1895-2006.
 
// From the command line, invoke as "makeNCfile" to get instructions. 

#define MISSING_VALUE -9999

// argv[0] = command
// argv[TIME_ARG] = <TIMELESS | 1895_2006>
// argv[GRID_ARG] = <CA800m | PNW800m | Yellowstone800m>
// argv[VAR_ARG] = <variable name> 
// argv[PREFIX_ARG] = <input_file_prefix> 
// argv[SUFFIX_ARG] = <input_file_suffix>
// argv[HDR_LEN_ARG] = <number of header records>
// argv[ROW_OFFSET_ARG] = <row_offset> // of the text input files
// argv[COL_OFFSET_ARG] = <col_offset> // of the text input files
// argv[NROWS_INPUT_ARG] = <nrows_input> // of the text input files
// argv[NCOLS_INPUT_ARG] = <ncols_input> // of the text input files
// argv[DIVISOR_ARG] = <divisor>  Used to downscale the input data before writing it out to the .nc file
// argv[VALID_MIN_ARG] = <valid_min>  Applies to data after dividing by divisor, and before multiplying by scale_factor.
// argv[VALID_MAX_ARG] = <valid_max>  Applies to data after dividing by divisor, and before multiplying by scale_factor.
// argv[UNITS_ARG] = <units>  A string.  These are the units that apply to valid_min and valid_max, and to the data in the .nc file
//     after it has been multiplied by scale_factor.
// argv[SCALE_FACTOR_ARG] = <scale_factor>
// argv[NOTE_ARG] = <note>  A string describing the data, especially its source.
#define TIME_ARG 1
#define GRID_ARG TIME_ARG+1
#define VAR_ARG GRID_ARG+1
#define PREFIX_ARG VAR_ARG+1
#define SUFFIX_ARG PREFIX_ARG+1
#define HDR_LEN_ARG SUFFIX_ARG+1
#define ROW_OFFSET_ARG HDR_LEN_ARG+1
#define COL_OFFSET_ARG ROW_OFFSET_ARG+1
#define NROWS_INPUT_ARG COL_OFFSET_ARG+1
#define NCOLS_INPUT_ARG NROWS_INPUT_ARG+1
#define DIVISOR_ARG NCOLS_INPUT_ARG+1
#define VALID_MIN_ARG DIVISOR_ARG+1
#define VALID_MAX_ARG VALID_MIN_ARG+1
#define UNITS_ARG VALID_MAX_ARG+1
#define SCALE_FACTOR_ARG UNITS_ARG+1
#define NOTE_ARG SCALE_FACTOR_ARG+1
#define NARGS NOTE_ARG+1
  
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


void instructions();


void skip_cols(int cols_to_skip, FILE * infile)
{
  int icol;
  int val;
  
  icol = 0;
  while (icol<cols_to_skip && fscanf(infile, "%d,", &val)>0) icol++;
  if (icol<cols_to_skip) err_exit("Error while skipping columns.");
} // end of skip_cols()


void skip_rows(int rows_to_skip, FILE * infile)
{
  int irow;
  char line[6961]; // 1392 columns x 5 chars per column + 1
  
  irow = 0;
  while (irow<rows_to_skip && fgets(line, sizeof(line) - 1, infile)!=NULL) irow++;  
  if (irow<rows_to_skip)
  {
	printf("\n*** " COMMAND_NAME ": Not enough rows in file.\n\n");
	fclose(infile);
	exit(-1); 
  }
} // end of skip_rows()

char header_text[MAX_HDR_LINES*(MAX_HDR_LINE_LEN+1)];
char six_row_header_text[MAX_HDR_LINES*(MAX_HDR_LINE_LEN+1)];
char ten_row_header_text[MAX_HDR_LINES*(MAX_HDR_LINE_LEN+1)];
// Boolean first_time_flag = TRUE;

void read_header(int hdr_len, FILE * infile)
{
  int irow;
  char line[6961]; // 1392 columns x 5 chars per column + 1
  Boolean read_OK_flag;
  
  assert(hdr_len<=MAX_HDR_LINES);
  irow = 0;
  read_OK_flag = TRUE;
  while (irow<hdr_len && read_OK_flag) 
  {
    read_OK_flag = fgets(line, sizeof(line) - 1, infile)!=NULL;
    if (read_OK_flag)
    {
      assert(strlen(line)<=MAX_HDR_LINE_LEN);
      printf(" header line %d: %s", irow + 1, line);
      strcpy(&header_text[irow*(MAX_HDR_LINE_LEN + 1)], line);
      irow++;  
    }
  }
  
  if (irow<hdr_len)
  {
	printf("\n*** " COMMAND_NAME ": Fewer rows in file than specified for header length.\n\n");
	fclose(infile);
	exit(-1); 
  }
} // end of read_header()

  
Boolean check_and_skip_header(int * hdr_lenP, FILE * infile)
{
  int irow; // , alt_irow;
  char line[6961]; // 1392 columns x 5 chars per column + 1
  Boolean read_OK_flag, header_OK_flag;
  int orig_hdr_len;
  fpos_t orig_file_pos;
  int rtnval;
  
  rtnval = fgetpos(infile, &orig_file_pos);
  assert(rtnval==0);
  orig_hdr_len = *hdr_lenP;
  assert(orig_hdr_len<=MAX_HDR_LINES);
  irow = 0;
  header_OK_flag = read_OK_flag = TRUE;
  while (irow<orig_hdr_len && read_OK_flag) 
  {
    read_OK_flag = fgets(line, sizeof(line) - 1, infile)!=NULL;
    if (read_OK_flag)
    {
      if (strcmp(line, &header_text[irow*(MAX_HDR_LINE_LEN + 1)])!=0)
       {
        header_OK_flag = FALSE;
        printf("\n*** " COMMAND_NAME ": Header mismatch.\n"
            "    header line #%d\n"
            "    orig: %s\n"
            "    new:  %s\n"
            "    How many lines are there in the new header?\n",
            irow + 1, &header_text[irow*(MAX_HDR_LINE_LEN + 1)], line);
/*
        if (first_time_flag)
        {
          first_time_flag = FALSE;
          assert(irow==0);
          assert(orig_hdr_len==6);
          for (alt_irow = 0; alt_irow<orig_hdr_len; alt_irow++)
              strcpy(&six_row_header_text[alt_irow*(MAX_HDR_LINE_LEN + 1)], 
              &header_text[alt_irow*(MAX_HDR_LINE_LEN + 1)]);
          fsetpos(infile, &orig_file_pos);
          *hdr_lenP = 10;
          printf("    New header length is %d\n", *hdr_lenP);
          read_header(*hdr_lenP, infile);
          for (alt_irow = 0; alt_irow<*hdr_lenP; alt_irow++)
              strcpy(&ten_row_header_text[alt_irow*(MAX_HDR_LINE_LEN + 1)], 
              &header_text[alt_irow*(MAX_HDR_LINE_LEN + 1)]);
          return(header_OK_flag);
        }
        else
        {
          assert(irow==0);
          *hdr_lenP = orig_hdr_len==6 ? 10 : 6;
          printf("    New header length is %d\n", *hdr_lenP);
          for (alt_irow = 0; alt_irow<*hdr_lenP; alt_irow++)
              strcpy(&header_text[alt_irow*(MAX_HDR_LINE_LEN + 1)], 
              *hdr_lenP==6 ? 
              &six_row_header_text[alt_irow*(MAX_HDR_LINE_LEN + 1)] :
              &ten_row_header_text[alt_irow*(MAX_HDR_LINE_LEN + 1)]);
          orig_hdr_len = *hdr_lenP;
        }
*/
      } // end of if (strcmp(line, &header_text[irow*(MAX_HDR_LINE_LEN + 1)])!=0)
      irow++;
    } // end of if (read_OK_flag)
    else header_OK_flag = FALSE;
  } // end of while (irow<orig_hdr_len && read_OK_flag)
  
  if (irow<orig_hdr_len)
  {
    printf("\n*** " COMMAND_NAME ": Fewer rows in file than specified for header length.\n\n");
    fclose(infile);
    header_OK_flag = FALSE; 
  }
  
  return(header_OK_flag);
} // end of check_and_skip_header()

  
#endif	/* MAKENCFILE_H */
 

