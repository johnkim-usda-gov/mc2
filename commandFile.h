/*
 *  commandFile.h
 *  MC2
 */

#define STRLEN				1024
#define PATHNAME_LENGTH	1024

typedef enum {rpUnknown, rpString, rpInt, rpFloat, rpCodeFlag} rpFcnEnum; 
typedef struct
{
	const char * keyword;
	rpFcnEnum fcnType;
	int num_valid_occur;
} ParamStruct;


typedef struct
{
	const char * name;
	int max_num;
	double default_value;
	float lower_limit;
	float upper_limit;
	bool inclusive_lower_limit;
	bool inclusive_upper_limit;
	const char * msg;  
} FloatParamInfoStruct;


typedef struct
{
	const char * name;
	int max_num;
	double default_value;
	float lower_limit;
	float upper_limit;
	const char * msg;  
} IntParamInfoStruct;


typedef enum { sptArbitrary, sptOnOff, sptPath, sptFileName, sptVarName } StringParamType;
typedef struct
{
	const char * name;
	int max_length;
	const char * default_value;
	StringParamType param_type;
	const char * msg;  
} StringParamInfoStruct;


/* rpFcn - function for processing a run parameter
   char * rpFcn(keyword, num_valid_occur, input_buf, RunDataP)
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
//typedef char * (*rpFcn)(char *, int, char *);


