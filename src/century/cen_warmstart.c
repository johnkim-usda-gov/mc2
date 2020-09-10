// cen_warmstart.c

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "netcdf.h"
#include "category_bgc.h"


#define COMMAND_NAME "mc2"

void err_exit(const char * msg)
   {
   printf("\n*** " COMMAND_NAME ": %s\n", msg);
   exit(-1);
   } // end of err_exit


void assert_(int flag)
{
  if (!flag)
  {
    assert(0);
  }
} /* end of assert_() */


void write_zeroes_(float * adrStart, float * adrEnd)
{
  memset(adrStart, 0, (adrEnd - adrStart)*sizeof(float));
} /* end of write_zeroes() */


void cpp_read_(int * fcnP, char * fileName, float * valP, char * valName, int * maxLenP, int * returnValP)
// fcn = 1 => open file
// fcn = 2 => read value and name
// fcn = 3 => read value only
// fcn = 4 => close file
// maxLen is the size of the valName character string in FORTRAN
// returnVal = 1 => normal return
// returnVal = 0 => couldn't open file or couldn't read file or file is at EOF
   {
   static FILE * fileP = NULL;
   char headerStr[81];
   char localName[81];
   
   if (*fcnP==1)
      { // First call includes the file name.
      assert((fileName!=NULL && strlen(fileName)>0));
      fileP = fopen(fileName, "r");
      fgets(headerStr, sizeof(headerStr) - 1, fileP);
      *returnValP = 1;
      }
   else if (*fcnP==2)
      { // Read a value and a name.
      *valP = NC_FILL_FLOAT;
      localName[0] = 0;
      fscanf(fileP, "%f%s", valP, localName);
      
      // Trim off the whitespace before and after the name, and put it in the calling program.
      // First, trim off any whitespace on the right end.
      int remaining_length = strlen(localName);
      while (remaining_length>0 && isspace(localName[remaining_length - 1])) remaining_length--;
      localName[remaining_length] = 0;
      // Now, count the characters of whitespace on the left end.
      int whitespace_chars = 0;
      while (whitespace_chars<remaining_length && isspace(localName[whitespace_chars])) whitespace_chars++;
      // Finally, copy the name to the string in the calling program
      int net_length = remaining_length - whitespace_chars;
      assert(net_length>=0);
      int i = 0;
      while (i<net_length && i<(*maxLenP)) 
         {
         *(valName + i) = localName[whitespace_chars + i];
         i++;
         }
      while (i<(*maxLenP)) 
         { // pad out the FORTRAN character string with blanks
         *(valName + i) = ' ';
         i++;
         }
      }
   else if (*fcnP==3)
      { // Read just a value.
      *valP = NC_FILL_FLOAT;
      fscanf(fileP, "%f", valP);
      }
   else if (*fcnP==4) fclose(fileP);
   else assert(0);
      
   return;
   } // end of cpp_read_()


int setCenturyParameter(CenParamStruct paramArray[], const char * paramNameIn, int valNdx, float paramVal)
   {
   int arrayNdx = 0;
   while (strcmp(paramArray[arrayNdx].paramName, "END")!=0)
      {
      if (strcmp(paramArray[arrayNdx].paramName, paramNameIn)==0)
         {
         paramArray[arrayNdx].paramVals[valNdx] = paramVal;
         return(1);
         }
      arrayNdx++;
      }
   return(0);
   } // end of setCenturyParameter()


int getCenturyParameter(CenParamStruct paramArray[], const char * paramNameIn, int valNdx, float * paramValP)
   {
   int arrayNdx = 0;
   while (strcmp(paramArray[arrayNdx].paramName, "END")!=0)
      {
      if (strcmp(paramArray[arrayNdx].paramName, paramNameIn)==0)
         {
         *paramValP = paramArray[arrayNdx].paramVals[valNdx];
         return(1);
         }
      arrayNdx++;
      }
   return(0);
   } // end of getCenturyParameter()


int treeLineNum;
int cropLineNum;      

void reset_tree_params_() { treeLineNum = 0; }
void reset_crop_params_() { cropLineNum = 0; }

void cen_init_parameters()
   {
   reset_tree_params_();
   reset_crop_params_();
   } // end of cen_init_parameters()
   
   
void cpy2FORTRANstr(char * FORTRAN_string, char * C_string, int FORTRAN_string_size)
   {

   int C_string_len = strlen(C_string);
   int len = C_string_len<FORTRAN_string_size ? C_string_len : FORTRAN_string_size;
   memcpy(FORTRAN_string, C_string, len);
   while (len<FORTRAN_string_size)
      { // pad out the FORTRAN string with blanks
      FORTRAN_string[len] = ' ';
      len++;
      }

   } // end of copy2FORTRANstr()


CenParamStruct pprdwcParams[] = { // coefficients for pprdwc() function
{"pprdwc_grass", {0.0, 1.0, 0.6, NC_FILL_FLOAT, NC_FILL_FLOAT}},
{"pprdwc_tree", {0.5, 1.0, 0.9, NC_FILL_FLOAT, NC_FILL_FLOAT}}};


void get_pprdwc3_(float * val1P, float * val2P, float * val3P, int * recordNdxP) 
   {
   int recordNdx = *recordNdxP;
   assert(recordNdx>=0 && recordNdx<(sizeof(pprdwcParams)/sizeof(CenParamStruct)));
   *val1P = pprdwcParams[recordNdx].paramVals[0];
   *val2P = pprdwcParams[recordNdx].paramVals[1];
   *val3P = pprdwcParams[recordNdx].paramVals[2];
   } // end of get_pprdwc3_()


CenParamStruct cropParams[] = { // data from crop.100
{"SUPRG    -C3-     -C4-",},
{"PRDX(1)", {175,     150.,     200.}}, //         PRDX(1)       
{"PPDF(1)", {30.,      18.,      30.}}, //          PPDF(1)       
{"PPDF(2)", {45.,      32.,      45.}}, //          PPDF(2)       
{"PPDF(3)", {1.,       1.2,      1.0}}, //          PPDF(3)       
{"PPDF(4)", {3.0,      3.0,      3.0}}, //          PPDF(4)       
{"BIOFLG", {1.}}, //                             BIOFLG        
{"BIOK5", {60.}}, //                            BIOK5         
{"PLTMRF", {1.}}, //                             PLTMRF        
{"FULCAN", {100.}}, //                           FULCAN        
{"FRTC(1)", {0.}}, //                             FRTC(1)       
{"FRTC(2)", {0.}}, //                             FRTC(2)   
{"FRTC(3)", {0.}}, //                             FRTC(3)       
{"BIOMAX", {400.}}, //                           BIOMAX        
{"PRAMN(1,1)", {20.,      20.,      20.}}, //          PRAMN(1,1)    
{"PRAMN(2,1)", {0.,       0.,       0.}}, //           PRAMN(3,1)    
{"PRAMN(3,1)", {0.,       0.,       0.}}, //           PRAMN(2,1)    
{"PRAMN(1,2)", {30.,      30.,      30.}}, //           PRAMN(1,2)    
{"PRAMN(2,2)", {0.,       0.,       0.}}, //           PRAMN(2,2)    
{"PRAMN(3,2)", {0.,       0.,       0.}}, //           PRAMN(3,2)    
{"PRAMX(1,1)", {30.,      30.,      30.}}, //          PRAMX(1,1)    
{"PRAMX(2,1)", {0.,       0.,       0.}}, //           PRAMX(2,1)    
{"PRAMX(3,1)", {0.,       0.,       0.}}, //           PRAMX(3,1)    
{"PRAMX(1,2)", {40.,      40.,      80.}}, //          PRAMX(1,2)    
{"PRAMX(2,2)", {0.,       0.,       0.}}, //           PRAMX(2,2)    
{"PRAMX(3,2)", {0.,       0.,       0.}}, //           PRAMX(3,2)    
{"PRBMN(1,1)", {40.,      40.,      60.}}, //          PRBMN(1,1)    
{"PRBMN(2,1)", {0.,       0.,       0.}}, //           PRBMN(2,1)    
{"PRBMN(3,1)", {0.,       0.,       0.}}, //           PRBMN(3,1)    
{"PRBMN(1,2)", {0.,       0.,       0.}}, //           PRBMN(1,2)    
{"PRBMN(2,2)", {0.,       0.,       0.}}, //           PRBMN(2,2)    
{"PRBMN(3,2)", {0.,       0.,       0.}}, //           PRBMN(3,2)    
{"PRBMX(1,1)", {50.,      50.,      80.}}, //          PRBMX(1,1)    
{"PRBMX(2,1)", {0.,       0.,       0.}}, //           PRBMX(2,1)    
{"PRBMX(3,1)", {0.,       0.,       0.}}, //           PRBMX(3,1)    
{"PRBMX(1,2)", {0.,       0.,       0.}}, //           PRBMX(1,2)    
{"PRBMX(2,2)", {0.,       0.,       0.}}, //           PRBMX(2,2)    
{"PRBMX(3,2)", {0.,       0.,       0.}}, //           PRBMX(3,2)    
{"FLIGNI(1,1)", {0.02}}, //                           FLIGNI(1,1)   
{"FLIGNI(2,1)", {0.0012}}, //                         FLIGNI(2,1)   
{"FLIGNI(1,2)", {0.26}}, //                           FLIGNI(1,2)   
{"FLIGNI(2,2)", {-0.0015}}, //                        FLIGNI(2,2)   
{"HIMAX", {0.}}, //                             HIMAX         
{"HIWSF", {0.}}, //                             HIWSF         
{"HIMON(1)", {0.}}, //                             HIMON(1)      
{"HIMON(2)", {0.}}, //                             HIMON(2)      
{"EFRGRN(1)", {0.5}}, //                            EFRGRN(1)     
{"EFRGRN(2)", {0.5}}, //                            EFRGRN(2)     
{"EFRGRN(3)", {0.5}}, //                            EFRGRN(3)     
{"VLOSSP", {0.04}}, //                           VLOSSP        
{"FSDETH(1)", {0.2}}, //                            FSDETH(1)     
{"FSDETH(2)", {0.95}}, //                           FSDETH(2)     
{"FSDETH(3)", {0.2}}, //                            FSDETH(3)     
{"FSDETH(4)", {150.}}, //                           FSDETH(4)     
{"FALLRT", {0.2}}, //                            FALLRT        
{"RDR", {0.05}}, //                           RDR           
{"RTDTMP", {2.}}, //                             RTDTMP        
{"CRPRTF(1)", {0.5}}, //                            CRPRTF(1)     
{"CRPRTF(2)", {0.}}, //                             CRPRTF(2)     
{"CRPRTF(3)", {0.}}, //                             CRPRTF(3)     
{"SNFXMX(1)", {0.000}}, //                          SNFXMX(1)     
{"DEL13C", {0.}}, //                             DEL13C        
{"CO2IPR", {1.25}}, //                           CO2IPR        
{"CO2ITR", {0.75}}, //                           CO2ITR        
{"CO2ICE(1,1,1)", {1.25}}, //                           CO2ICE(1,1,1)  
{"CO2ICE(1,1,2)", {1.}}, //                             CO2ICE(1,1,2)  
{"CO2ICE(1,1,3)", {1.}}, //                             CO2ICE(1,1,3)  
{"CO2ICE(1,2,1)", {1.25}}, //                           CO2ICE(1,2,1)  
{"CO2ICE(1,2,2)", {1.}}, //                             CO2ICE(1,2,2)  
{"CO2ICE(1,2,3)", {1.}}, //                             CO2ICE(1,2,3)  
{"CO2IRS", {1.}}, //                             CO2IRS        
{"END", },
};
   

CenParamStruct treeParams[] = { // data from tree.100
{"SUPRT	-DN-	-EN-	-DB-	-EB-",},	
{"DECID", {1.0}},					
{"PRDX(3)", {10000.,	10000.,	10000.,	10000.,	10000.}}, // PRDX3 gDM/month max gross production	
{"PRDX(4)", {250.,	250.,	250.,	250.,	250.}}, //	'PRDX(4)' gC/month max net production
{"PPDF(1)", {22.0,	15.,	15.,	22.,	30.}}, //	'PPDF(1)'
{"PPDF(2)", {38.0,	30.,	30.,	35.,	45.}}, //	'PPDF(2)'
{"PPDF(3)", {1.1,	.50,	.50,	.20,	.50}}, //	'PPDF(3)'
{"PPDF(4)", {3.0,	5.0,	5.0,	4.5,     4}}, //	'PPDF(4)'
{"CERFOR(1,1,1)", {60.0,	100.0,	100.0,	20.0,	20.0}}, //	'CERFOR(1,1,1)'
{"CERFOR(1,1,2)", {0.0,	0.0,	0.0,	0.0,	0.0}}, //	'CERFOR(1,1,2)'
{"CERFOR(1,1,3)", {0.0,	0.0,	0.0,	0.0,	0.0}}, //	'CERFOR(1,1,3)'
{"CERFOR(1,2,1)", {45.0,	50.0,	50.0,	35.0,	35.0}}, //	'CERFOR(1,2,1)'
{"CERFOR(1,2,2)", {0.0,	0.0,	0.0,	0.0,	0.0}}, //	'CERFOR(1,2,2)'
{"CERFOR(1,2,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(1,2,3)'
{"CERFOR(1,3,1)", {204.,	310.,	310.,	80.0,	120.0}}, //	'CERFOR(1,3,1)'
{"CERFOR(1,3,2)", {0.0,	0.0,	0.0,	0.0,	0.0}}, //	'CERFOR(1,3,2)'
{"CERFOR(1,3,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(1,3,3)'
{"CERFOR(1,4,1)", {260.0,	900.0,	900.0,	140.0,	150.0}}, //	'CERFOR(1,4,1)'
{"CERFOR(1,4,2)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(1,4,2)'
{"CERFOR(1,4,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(1,4,3)'
{"CERFOR(1,5,1)", {127.,	600.,	600.,	83.,	150.}}, //	'CERFOR(1,5,1)'
{"CERFOR(1,5,2)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(1,5,2)'
{"CERFOR(1,5,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(1,5,3)'
{"CERFOR(2,1,1)", {70.0,	100.0,	100.0,	40.0,	40.0}}, //	'CERFOR(2,1,1)'
{"CERFOR(2,1,2)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(2,1,2)'
{"CERFOR(2,1,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(2,1,3)'
{"CERFOR(2,2,1)", {50.0,	81.0,	81.0,	50.0,	60.0}}, //	'CERFOR(2,2,1)'
{"CERFOR(2,2,2)", {0.0,	0.0,	0.0,	0.0,	0.0}},  //	'CERFOR(2,2,2)'
{"CERFOR(2,2,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(2,2,3)'
{"CERFOR(2,3,1)", {204.,	310.,	310.,	99.,	180.}}, //	'CERFOR(2,3,1)'
{"CERFOR(2,3,2)", {0.0,	0.0,	0.0,	0.0,	0.0}}, //	'CERFOR(2,3,2)'
{"CERFOR(2,3,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(2,3,3)'
{"CERFOR(2,4,1)", {260.0,	800.0,	800.0,	140.0,	300.0}}, //	'CERFOR(2,4,1)'
{"CERFOR(2,4,2)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(2,4,2)'
{"CERFOR(2,4,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(2,4,3)'
{"CERFOR(2,5,1)", {126.,	80.,	80.,	500.,	300.}}, //	'CERFOR(2,5,1)'
{"CERFOR(2,5,2)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(2,5,2)'
{"CERFOR(2,5,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(2,5,3)'
{"CERFOR(3,1,1)", {80.0,	90.0,	90.0,	40.0,	40.0}}, //	'CERFOR(3,1,1)'
{"CERFOR(3,1,2)", {0.0,	0.0,	0.0,	0.0,	0.0}}, //	'CERFOR(3,1,2)'
{"CERFOR(3,1,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(3,1,3)'
{"CERFOR(3,2,1)", {50.0,	80.0,	80.0,	50.0,	76.0}}, //	'CERFOR(3,2,1)'
{"CERFOR(3,2,2)", {0.0,	0.0,	0.0,	0.0,	0.0}}, //	'CERFOR(3,2,2)'
{"CERFOR(3,2,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(3,2,3)'
{"CERFOR(3,3,1)", {204.,	300.,	300.,	80.,	84.}}, //	'CERFOR(3,3,1)'
{"CERFOR(3,3,2)", {0.0,	0.0,	0.0,	0.0,	0.0}}, //	'CERFOR(3,3,2)'
{"CERFOR(3,3,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(3,3,3)'
{"CERFOR(3,4,1)", {260.0,	900.0,	900.0,	140.0,	155.0}}, //	'CERFOR(3,4,1)'
{"CERFOR(3,4,2)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(3,4,2)'
{"CERFOR(3,4,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(3,4,3)'
{"CERFOR(3,5,1)", {126.,	550.,	550.,	80.0,	155.0}}, //	'CERFOR(3,5,1)'
{"CERFOR(3,5,2)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(3,5,2)'
{"CERFOR(3,5,3)", {0.,	0.,	0.,	0.,	0.}}, //	'CERFOR(3,5,3)'
{"DECW1", {0.9}}, //					'DECW1'
{"DECW2", {0.4}}, //					'DECW2'
{"DECW3", {0.4}}, //					'DECW3'
{"FCFRAC(1,1)", {0.36,	0.37,	0.37,	0.34,	0.25}}, //	'FCFRAC(1,1)'
{"FCFRAC(2,1)", {0.37,	0.34,	0.34,	0.40,	0.25}}, //	'FCFRAC(2,1)'
{"FCFRAC(3,1)", {0.09,	0.10,	0.10,	0.09,	0.10}}, //	'FCFRAC(3,1)'
{"FCFRAC(4,1)", {0.16,	0.18,	0.18,	0.15,	0.30}}, //	'FCFRAC(4,1)'
{"FCFRAC(5,1)", {0.02,	0.01,	0.01,	0.02,	0.10}}, //	'FCFRAC(5,1)'
{"FCFRAC(1,2)", {0.33,	0.37,	0.37,	0.34,	0.34}}, //	'FCFRAC(1,2)'
{"FCFRAC(2,2)", {0.3,	0.34,	0.34,	0.40,	0.25}}, //	'FCFRAC(2,2)'
{"FCFRAC(3,2)", {0.09,	0.10,	0.10,	0.09,	0.11}}, //	'FCFRAC(3,2)'
{"FCFRAC(4,2)", {0.25,	0.18,	0.18,	0.15,	0.22}}, //	'FCFRAC(4,2)'
{"FCFRAC(5,2)", {0.03,	0.01,	0.01,	0.02,	0.08}}, //	'FCFRAC(5,2)'
{"LEAFDR(1)", {0.02,	0.00,	0.03,	0.00,	0.07}}, //	'LEAFDR(1)'
{"LEAFDR(2)", {0.02,	0.00,	0.03,	0.00,	0.07}}, //	'LEAFDR(2)'
{"LEAFDR(3)", {0.02,	0.00,	0.03,	0.00,	0.07}}, //	'LEAFDR(3)'
{"LEAFDR(4)", {0.02,	0.00,	0.03,	0.00,	0.07}}, //	'LEAFDR(4)'
{"LEAFDR(5)", {0.02,	0.00,	0.03,	0.00,	0.07}}, //	'LEAFDR(5)'
{"LEAFDR(6)", {0.02,	0.00,	0.03,	0.00,	0.07}}, //	'LEAFDR(6)'
{"LEAFDR(7)", {0.02,	0.00,	0.03,	0.00,	0.07}}, //	'LEAFDR(7)'
{"LEAFDR(8)", {0.02,	0.00,	0.03,	0.00,	0.07}}, //	'LEAFDR(8)'
{"LEAFDR(9)", {0.02,	0.00,	0.03,	0.00,	0.07}}, //	'LEAFDR(9)'
{"LEAFDR(10)", {0.02,	0.00,	0.03,	0.00,	0.07}}, //	'LEAFDR(10)'
{"LEAFDR(11)", {0.02,	0.00,	0.03,	0.00,	0.07}}, //	'LEAFDR(11)'
{"LEAFDR(12)", {0.02,	0.00,	0.10,	0.00,	0.07}}, //	'LEAFDR(12)'
{"BTOLAI", {0.008,	0.012,	0.004,	0.012,	0.007}}, //	'BTOLAI'
{"KLAI", {2000.0,	2000.0,	2000.0,	1000.0,	1000.0}}, //	'KLAI'
{"LAITOP", {-0.47000}}, //				'LAITOP'
{"MAXLAI", {8.,	8.,	8.,	8.,	8.}}, //	'MAXLAI'
{"MAXLDR", {1.0}}, //					'MAXLDR'
{"FORRTF(1)", {0.45000}}, //					'FORRTF(1)'
{"FORRTF(2)", {0.0}}, //					'FORRTF(2)'
{"FORRTF(3)", {0.0}}, //					'FORRTF(3)'
{"SAPK", {1500.00}}, //					'SAPK'
{"SWOLD", {0.0}}, //					'SWOLD'
{"WDLIG(1)", {0.21000}}, //					'WDLIG(1)'
{"WDLIG(2)", {0.22000}}, //					'WDLIG(2)'
{"WDLIG(3)", {0.25000}}, //					'WDLIG(3)'
{"WDLIG(4)", {0.30000}}, //					'WDLIG(4)'
{"WDLIG(5)", {0.30000}}, //					'WDLIG(5)'
{"WOODDR(1)", {1.0,	1.0,	0.0,	1.0,	0.0}}, //	'WOODDR(1)'
{"WOODDR(2)", {0.032,	0.05,	0.05,	0.04,	0.03}}, //	'WOODDR(2)'
{"WOODDR(3)", {0.040,	0.01,	0.01,	0.01,	0.01}}, //	'WOODDR(3)'
{"WOODDR(4)", {0.0012,	0.0008,	0.0008,	0.002,	0.002}}, //	'WOODDR(4)'
{"WOODDR(5)", {0.0021,	0.001,	0.001,	0.004,	0.004}}, //	'WOODDR(5)'
{"SNFXMX(2)", {0.000}}, //   				'SNFXMX(2)'
{"DEL13C", {0.0}}, //					'DEL13C'
{"CO2IPR", {1.25}}, //					'CO2IPR'
{"CO2ITR", {0.75}}, //					'CO2ITR'
{"CO2ICE(1,1,1)", {1.25}}, //					'CO2ICE(1,1,1)'
{"CO2ICE(1,1,2)", {1.0}}, //					'CO2ICE(1,1,2)'
{"CO2ICE(1,1,3)", {1.0}}, //					'CO2ICE(1,1,3)'
{"CO2ICE(1,2,1)", {1.25}}, //					'CO2ICE(1,2,1)'
{"CO2ICE(1,2,2)", {1.0}}, //					'CO2ICE(1,2,2)'
{"CO2ICE(1,2,3)", {1.0}}, //					'CO2ICE(1,2,3)'
{"CO2IRS", {1.0}}, //					'CO2IRS'
{"BASFC2", {1.0}}, //					'BASFC2'
{"BASFCT", {400.0}}, //					'BASFCT'
{"SITPOT", {2400.0,	4800.0,	4800.0,	2400.0,	2400.0}}, //	'SITPOT'
{"END", },
};
   
   
void get_tree1_(float * valP, char * name) 
   {
   assert(treeLineNum>=0 && treeLineNum<(sizeof(treeParams)/sizeof(CenParamStruct)));
   *valP = treeParams[treeLineNum].paramVals[0];
   cpy2FORTRANstr(name, treeParams[treeLineNum].paramName, 20);
   
   treeLineNum++;
   } // end of getTree1_()


void get_tree5_(float * val1P, float * val2P, float * val3P, float * val4P, float * val5P, char * name) 
   {
   assert(treeLineNum>=0 && treeLineNum<(sizeof(treeParams)/sizeof(CenParamStruct)));
   *val1P = treeParams[treeLineNum].paramVals[0];
   *val2P = treeParams[treeLineNum].paramVals[1];
   *val3P = treeParams[treeLineNum].paramVals[2];
   *val4P = treeParams[treeLineNum].paramVals[3];
   *val5P = treeParams[treeLineNum].paramVals[4];
   cpy2FORTRANstr(name, treeParams[treeLineNum].paramName, 20);
   
   treeLineNum++;
   } // end of getTree5_()
   
   
void get_crop1_(float * valP, char * name) 
   {
   assert(cropLineNum>=0 && cropLineNum<(sizeof(cropParams)/sizeof(CenParamStruct)));
   *valP = cropParams[cropLineNum].paramVals[0];
   cpy2FORTRANstr(name, cropParams[cropLineNum].paramName, 6);
   
   cropLineNum++;
   } // end of get_crop1_()


void get_crop3_(float * val1P, float * val2P, float * val3P, char * name) 
   {
   assert(cropLineNum>=0 && cropLineNum<(sizeof(cropParams)/sizeof(CenParamStruct)));
   *val1P = cropParams[cropLineNum].paramVals[0];
   *val2P = cropParams[cropLineNum].paramVals[1];
   *val3P = cropParams[cropLineNum].paramVals[2];
   cpy2FORTRANstr(name, cropParams[cropLineNum].paramName, 6);
   
   cropLineNum++;
   } // end of get_crop3_()
   
   
CO2Struct atmosCO2conc;

void get_co2_(int * yrP, float * CO2concP)
   {
   int yrNdx = *yrP - 1; // FORTRAN starts at 1, C starts at 0
   assert(0<=yrNdx && yrNdx<(sizeof(atmosCO2conc.yearlyCO2conc)/sizeof(float)));
   *CO2concP = atmosCO2conc.yearlyCO2conc[yrNdx];
   } // end of get_co2_()
   

void set_constant_CO2(float CO2conc)
   {
   int ramp_len = sizeof(atmosCO2conc.yearlyCO2conc)/sizeof(float);
   int yrNdx;
   
   for (yrNdx = 0; yrNdx<ramp_len; yrNdx++) atmosCO2conc.yearlyCO2conc[yrNdx] = CO2conc;
         
   } // end of set_constant_CO2()
   
   
int read_co2ramp_C(char * co2ramp_name, int first_calendar_year, int min_years, int years_offset_into_input_data)
// Will read CO2 ramp files with 2 columns: calendar year, CO2 concentration
// or, for backward compatibility, CO2 ramp files with just 1 column: CO2 concentration.
// Eventually, lines which do not begin with a digit are treated as comments.
   { 
   FILE * fileP;
   char skippedLine[81];
   // char headerStr[81];
   // printf("read_co2ramp_C: co2ramp_name, first_calendar_year, min_years, years_offset_into_input_data = %s, %d, %d, %d\n",
   //      co2ramp_name, first_calendar_year, min_years, years_offset_into_input_data);
   
   float inVal;
   int ramp_len = sizeof(atmosCO2conc.yearlyCO2conc)/sizeof(float);

   atmosCO2conc.first_calendar_year = first_calendar_year;
   
   if (co2ramp_name==NULL) err_exit("co2ramp_name is NULL");
   else if (strlen(co2ramp_name)<=0) err_exit("strlen(co2ramp_name) = 0");
   else printf("co2ramp_name = %s\n", co2ramp_name);
   assert((co2ramp_name!=NULL && strlen(co2ramp_name)>0));
   
   fileP = fopen(co2ramp_name, "r");
   if (fileP==NULL) 
      {
      printf("read_co2ramp_C: couldn't open %s\n", co2ramp_name);
      return(0);
      }
      
   // fgets(headerStr, sizeof(headerStr), fileP);
   
   int lineCount = 0;
   while (lineCount<years_offset_into_input_data) 
      {
      fgets(skippedLine, sizeof(skippedLine), fileP); // Works as long as the lines in the file aren't longer than sizeof(skippedLine) - 1
      // printf("lineCount, skippedLine = %d, %s\n", lineCount, skippedLine);
      lineCount++;
      }
   if (feof(fileP))
      {
      printf("read_co2_ramp_C: encountered eof while skipping lines, at linecount = %d.  This can happen if the CO2ramp file has Windows-style endlines instead of Unix-style endlines.\n", lineCount);
      fclose(fileP);
      return(0);
      }

   int yrNdx = 0;
   while (!feof(fileP) && yrNdx<ramp_len)
      { int count;
      inVal = NC_FILL_FLOAT;
      count = fscanf(fileP, "%f", &inVal);
      // printf("read_co2ramp_C: inVal = %f\n", inVal);
      // If we couldn't read a number, return an error
      if (!feof(fileP) && (count<1 || inVal==NC_FILL_FLOAT)) 
         {
         fclose(fileP);
         return(0);
         }
      // If the value looks like the right calendar year, read the next value
      if (!feof(fileP) && inVal==(first_calendar_year + yrNdx))
         {
         count = fscanf(fileP, "%f", &inVal);
         // printf("read_co2ramp_C 2nd read: inVal = %f\n", inVal);
         }
      // If the value doesn't look like a CO2 concentration in parts per million, return an error.
      if (!feof(fileP) && (count!=1 || inVal<0. || inVal>1500.))
         {
         fclose(fileP);
         return(0);
         }
      if (count==1)
         {
         atmosCO2conc.yearlyCO2conc[yrNdx] = inVal;
         yrNdx++;
         }
      }
      
   fclose(fileP);
   
   if (yrNdx<min_years) return(0);
   
   assert(yrNdx>=1);
   float lastVal = atmosCO2conc.yearlyCO2conc[yrNdx - 1];
   while (yrNdx<ramp_len)
      { // Pad out the array with the last value read.
      atmosCO2conc.yearlyCO2conc[yrNdx] = lastVal;
      yrNdx++;
      }
      
   return(1);
   } // end of read_co2ramp_C()

   
   
   
   
   
