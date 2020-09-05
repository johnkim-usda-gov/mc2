/*
 *  MC2.h
 */

#define MONTH_INTERVAL 1
#define YEAR_INTERVAL 2 
#define MULTIYR_INTERVAL 3 
#define SINGLE_INTERVAL 4 /* use for constant input data, e.g. elevation */

bool chk_nc(int io_rtnval);
void err_exit(const char * msg);

typedef struct 
{
	int	fileid;
	int	row_dimid, row_varid, rowNdxPos;
	int	col_dimid, col_varid, colNdxPos;
	int	time_dimid, time_varid, timeNdxPos;
	int  ws_dimid, ws_varid, wsNdxPos;
	size_t row_dimlen, col_dimlen, time_dimlen, ws_dimlen; // ws_dimlen is usually NUM_WS_OUTVARS
	int vardimids[4];
	short row_offset, col_offset;
	int	input_data_varid;
	nc_type row_vartype, col_vartype, time_vartype, input_data_vartype;
	char * file_name_and_path;
	char * var_name;
	bool scale_factor_flag;
	float scale_factor;
	bool missing_value_flag;
	float missing_value;
	bool fill_value_flag;
	float fill_value;
	char * units;
} ncid_type;


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
typedef char * (*rpFcn)(char * keyword, int times_called, char * buf);

class OutVarDescription 
{
	public:
		OutVarDescription();
		void save(float inVal);   
		void save(float inVal, int interval);

	public:
		const char * descrip;
		nc_type type;  
		float scale;
		const char * name;
		int interval; /* use enum instead ??*/
		bool show;
		const char * units;
		int var_id;
		bool categoricalFlag;
		int max_categories;
		int moOutNdx, yrOutNdx, multiyrOutNdx, singleOutNdx, dictNdx;

		Simulation * m_pS;

}; // end of class OutVarDescription


class Simulation
{
	// Methods
	public:
		Simulation(int argc, char * argv[]);
		// Calls interpretCalibrationName(), getModelParameters(), getRunParameters().
		~Simulation(); 

		// entry point from C++ main()
		void run(); // Calls openInputFiles(), initializeOutputFiles(), maskedOut(), and saveOutputData().
		// Iteratively calls the runModelAndCollectData() method for each process model.

		// routines called by the constructor
		BaseCalibrationEnum interpretCalibrationName(char * calibration_name); // 1X
		bool getModelParameters(BaseCalibrationEnum modelPS); // Calls ModelParamsClass::interpretMAPSSparameterSet(), 
		// interpretSwitch(). 1X
		bool getRunParameters(char * cmd_file_name); // Calls rpStringFcn(), rpIntFcn(), rpFloatFcn(),
		// rpCodeFlagFcn(), process_command_file(), RunParamsClass::interpretRunMode(), 
		// RunParamsClass::interpretGridName(). 1x

		// routines called by run()
		bool openInputFiles(); // Calls openInputDataNCfile(), openWarmstartNCfile(). 1X
		void initializeOutputFiles(); // Calls initializeMAPSSoutputFile(), initializeEQoutputFile(),
		// initializeWarmstartOutputFile(). 1X
		// bool maskedOut(int row_index, int col_index);  // 1X
		void saveOutputData(); // Closes output netCDF files. 1X

		// routine called by the runModelAndCollectData() methods in the process models
		bool getInputData(int month, int grid_row, int grid_col, bool ckLatLonFlag, InputDataClass * pInputData); 
		// Calls north_south_centroid(), east_west_centroid(), getEarthData(), getClimateVar().

		// routines called by other methods
		void access_file(FILE ** file, char * filename, const char * status); // 2X
		bool addVarToList(char * name);
		void addVarToMoList(int varNdx);
		void addVarToMultiyrList(int varNdx);
		void addVarToSingleList(int varNdx);
		void addVarToYrList(int varNdx);
		void changeFRI(int vtype, int minFRI, int maxFRI);
		void defineEQvars(ncid_type * ncidP, OutVarDescription * OutputList, int filter_length); // Calls describeVar(). 2X
		void describeFSvars(OutVarDescription * FSoutputList); // 1X
		int defineNCvar(int fileid, int dimids[], const char * name, const char * description, const char * units, nc_type type); // 3X
		void defineNCvarsInList(ncid_type * ncidP, OutVarDescription * list, int filter_length);
		OutVarDescription describeVar(const char * name, const char * description, const char * units, nc_type type); // Called many times.
		double east_west_centroid(int col); // 1X
		void getAttInfo(const char * att_name, ncid_type * pNcid, bool * pFlag, float * pValue); // 2X
		bool getAttVal(ncid_type * ncidP, const char * att_name, short * short_attValP);
		bool getClimateVar(ncid_type * ncidP, float * varP, int row, int col, int num_months, bool ckLatLonFlag, double lat, double lon); // 6X
		bool getClimateVar(ncid_type * ncidP, float * varP, int grid_row, int grid_col, int num_months, bool ckLatLonFlag, double cell_lat, double cell_lon, int months_offset); // Calls verifyLatLon(). 2X
		bool getColInfo(ncid_type * ncidP); // 5X
		void getEarthData(InputDataClass * pInputData, SoilDataStruct * soil_dataP, int row, int col, bool ckLatLonFlag, double lat, double lon); // Calls verifyLatLon(). 1X
		bool getOffsets(ncid_type * ncidP);
		bool getRowInfo(ncid_type * ncidP); // 5X
		bool getTimeInfo(ncid_type * ncidP); // 5X
		char * getUnits(ncid_type * pNcid); // 1X
		void getVarInfo(ncid_type * ncidP); // Calls getAttInfo(), getUnits(). 1X
		int initializeActiveCellArray(); // 1X
		void initializeDiscretionaryOutputFile(ncid_type * ncidP, const char * suffix, OutVarDescription outList[], int num_vars,
				char * cmd_line, char * cmd_file_name); // 3X
		void initializeDiscretionaryOutputFiles(char * cmd_line, char * cmd_file_name); // 1X
		void initializeEQoutputFile(char * cmd_line, char * cmd_file_name); // Calls initializeEQoutputList(), makePath(), 
		// writeStdAttributes(), writeCommands(), writeRunParameters(), writeModelParameters(), getRowInfo(), 
		// getColInfo(), getTimeInfo(), defineEQvars(). 1X
		void initializeEQoutputList(OutVarDescription f[], int filter_length); // 3X
		void initializeMAPSSoutputFile(char * cmd_line, char * cmd_file_name); // Calls makePath(), writeStdAttributes(),
		// writeCommands(), writeRunParameters(), writeModelParameters(), getRowInfo(), getColInfo(), getTimeInfo(),
		// describeVar(). 1X
		void initializeMAPSSoutputList(); // Calls makeOutputListEntry(). 1X
		void initializeMonthlyOutputFile(ncid_type * ncidP, const char * suffix, OutVarDescription outList[], int num_vars,
				char * cmd_line, char * cmd_file_name, int num_months); // 1X
		void initializeVarDict(); // 1X
		void initializeWSoutputFile(char * cmd_line, char * cmd_file_name); // Calls initializeEQoutputList(), makePath(), 
		// writeStdAttributes(), writeCommands(), writeRunParameters(), writeModelParameters(), getRowInfo(), 
		// getColInfo(), getTimeInfo(), defineEQvars(), defineFSvars(). 2X
		void initializeWSoutputList(OutVarDescription f[], int filter_length); // 2X
		bool interpretSwitch(char * switchStr);  // 4X
		void makeCategoricalVarDictEntry(int varNdx, const char * var_name, const char * var_descrip, const char * var_units, int max_categories); // Called once per categorical output variable; 
		// calls makeVarDictEntry().
		void makeOutputListEntry(OutVarDescription * f, const char * name, const char * description, const char * units, nc_type type); // 4X
		char * makePath(char * prefix, const char * suffix); // 3X
		void makeVarDictEntry (int varNdx, const char * name, const char * description, const char * units, nc_type type, int interval);
		void MC2instructions();
		double north_south_centroid(int row); // 1X
		bool openInputDataNCfile(ncid_type * ncidP, char * directory, const char * file_name); 
		//  Calls same with 5 args. 10X
		bool openInputDataNCfile(ncid_type * ncidP, char * directory, const char * file_name, const char * var_name);  
		// Calls getRowInfo(), getColInfo(), getTimeInfo(), getVarInfo(). 3X
		bool openWarmstartNCfile(ncid_type * ncidP, char * file_name_and_path); // Calls getRowInfo(), getColInfo(), 
		// getTimeInfo(), initializeMAPSSoutputList(), initializeEQoutputList(). 3X
		int process_command_file(char * file_name_and_path, int num_valid_occur[]); // Calls read_param(). 1X
		void read_param(int num_valid_occur[], char * name, char * buf); // Calls rpStringFcn(), rpIntFcn(),
		// rpFloatFcn(), rpCodeFlagFcn(). 1X
		const char * rpCodeFlagFcn(const char * keyword, int n, char * buf); // 3X
		const char * rpFloatFcn(const char * keyword, int times_called, char * buf); // 3X
		const char * rpIntFcn(const char * keyword, int times_called, char * buf);  // 3X
		const char * rpStringFcn(const char * keyword, int times_called, char * buf);  // 3X
		void saveOutputFile(ncid_type * ncidP);
		void saveMultiyrOutputData();
		void saveOneMonthOneCellOutputData(int moNdx);
		void saveOneYearOneCellOutputData(int yrNdx);
		void setNdxPos(ncid_type * ncidP);
		void verifyLatLon(ncid_type * ncidP, int row, int col, double cell_lat, double cell_lon); 
		// Calls ScienceFcns::close_enough(). 4X
		void writeCommands(ncid_type * ncidP, char * cmd_line, char * cmd_file_name); // Calls access_file(). 3X
		void writeModelParameters(ncid_type * ncidP); // 3X
		void writeMonthlyData(size_t coords[]); // 2X?
		void writeMultiyrData(size_t coords[]); // 2X?
		void writeOutputVarList(const char * file_name, OutVarDescription List[], int list_len);
		void writeRunParameters(ncid_type * ncidP); // 3X
		void writeSingleData(size_t coords[]);
		void writeStdAttributes(ncid_type * ncidP); // 4X
		void writeStdAttributes(ncid_type * ncidP, size_t time_length); // 2X
		void writeYearlyData(size_t coords[]); // 2X?


		// Member Data
	public:
		char * m_cmdFileName;
		char * m_cmdLine;
		int m_numProcessModels;
		BaseCalibrationEnum m_baseCalibration;
		RunParamsClass runParams;
		ModelParamsClass modelParams; 
		ncid_type MAPSSoutFile_ncid;
		int mclassVarid;
		ncid_type EQoutFile_ncid;
		OutVarDescription MAPSSoutputList[NUM_MAPSS_OUTVARS];
		OutVarDescription EQandWSoutputList[NUM_WS_OUTVARS];
		OutVarDescription EQandWSinputList[NUM_WS_OUTVARS];
		OutVarDescription VarDict[NUM_OPTIONAL_OUTVARS];

		OutVarDescription m_moOutList[NUM_OPTIONAL_OUTVARS];
		int m_num_monthly_vars;
		ncid_type m_moOutFile_ncid;
		float m_moOutVars[NUM_OPTIONAL_OUTVARS];
		float * m_moBufP;

		OutVarDescription m_yrOutList[NUM_OPTIONAL_OUTVARS];
		int m_num_yearly_vars;
		ncid_type m_yrOutFile_ncid;
		float m_yrOutVars[NUM_OPTIONAL_OUTVARS];
		float * m_yrBufP;

		OutVarDescription m_multiyrOutList[NUM_OPTIONAL_OUTVARS];
		int m_num_multiyr_vars;
		ncid_type m_multiyrOutFile_ncid;
		float m_multiyrOutVars[NUM_OPTIONAL_OUTVARS];
		float * m_multiyrBufP;
		short * m_multiyrIntervalBufP;

		OutVarDescription m_singleOutList[NUM_OPTIONAL_OUTVARS];
		int m_num_single_vars;
		ncid_type m_singleOutFile_ncid;
		float m_singleOutVars[NUM_OPTIONAL_OUTVARS];

		ncid_type mask_ncid;
		ncid_type WSoutFile_ncid;
		ncid_type soilData_ncid;
		ncid_type bd_ncid;
		ncid_type elev_ncid;
		ncid_type ppt_ncid;
		ncid_type tmin_ncid;
		ncid_type tmax_ncid;
		ncid_type tmp_ncid;
		ncid_type vpr_ncid;
		ncid_type tdmean_ncid;
		ncid_type wnd_ncid;

		// Next 6 are for W_WA base calibration
		ncid_type fog_ncid; // fog layer 
		ncid_type topomoist_ncid; // topographic moisture layer 
		ncid_type deltaTsl_ncid; // difference between mean air temperature and mean air temperature at sea level layer
		ncid_type aspect_ncid; // aspect layer
		ncid_type sw_ncid; // shortwave layer
		ncid_type cad_ncid; // cold air drainage layer

		ncid_type MAPSSdata_ncid;
		ncid_type EQdata_ncid;
		ncid_type WSdata_ncid;
		int * m_pActiveCellArray;
		int m_maskVal;

		// enum MAPSSparameterSetName {defaultParams, ORNLparams, US10kmFAOparams, US10kmSCSparams, WiesYoseParams };
		char * MpSN[5]; // = {"defaultParams", "ORNLparams", "US10kmFAOparams", "US10kmSCSparams", "WiesYoseParams"}; 

		int m_num_vpr_tmp_issues;
		int m_cmd_file_len;
		int m_num_target_cells;
		int m_num_active_cells;
		int m_num_failed_cells;

	private:
		std::vector<ProcessModel *> * m_pProcessModelList;

}; // end of class Simulation




