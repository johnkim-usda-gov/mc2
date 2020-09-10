//
//  PNVmodel.h
//  MC2project
//
//  Created by David R Conklin on 7/18/14.
//


class PNVmodel: public ProcessModel
{
	PNVmodel() {} 
	public:
	PNVmodel(Simulation * sP, BiogeographyInputData * inValsP); 
	~PNVmodel() {}

	bool runModelAndCollectData(const int year_index, const int row_index, const int col_index); 

	PNVcode pnv(BiogeographyInputData * inValsP, ModelParamsClass * mP);
	PNVcode SubregionSpecificCode(BiogeographyInputData * inValsP);
	PNVcode Subregion10010103();
	PNVcode Subregion10010104(); 
	PNVcode Subregion10010203();
	PNVcode Subregion10010204();
	PNVcode Subregion10010205();
	PNVcode Subregion10010206();
	PNVcode Subregion10010207();
	PNVcode Subregion10010208(); 
	PNVcode Subregion10010210();
	PNVcode Subregion10010214();
	PNVcode Subregion10010215();
	PNVcode Subregion10020101();
	PNVcode Subregion10020110();
	PNVcode Subregion10020202();

	int con(bool testFlag, int true_answer, int false_answer) { return(testFlag ? true_answer : false_answer); }

	// inputs to the subregion specific routines
	float prcpgrid; // precipitation at sea level, in H2O
	float foggrid; // fog effect, each unit = 20" H2O
	float cadgrid; // cold air drainage effect, ft of elevation
	float msltgrid; // mean air temperature at sea level, deg F
	float tmgrid; // topographic moisture effect, ft of elevation
	float singrid; // aspect 
	float solargrid; // shortwave, mean daily clear-sky shortwave radiation, in kilojoules per square meter per day 
	float elevtmp; // elevation, ft

	// outputs from the subregion specific routines
	float vgtmp2; // alpine zone lower boundary, ft
	float vgtmp3; // parkland zone lower boundary, ft
	float vgtmp4; // mountain hemlock zone lower boundary, ft
	float vgtmp5; // Pacific silver fir zone lower boundary, ft
	float vgtmp6; // Sitka spruce zone upper boundary, ft
	float vgtmp7; // subalpine fir zone lower boundary, ft
	float vgtmp8; // subalpine fir zone upper boundary, ft
	float vgtmp9; // Douglas-fir zone upper boundary, ft
	float vgtmp10; // Douglas-fir zone lower boundary, ft
	float vgtmp11; // mountain hemlock upper boundary, ft
	float vgtmp12; // western hemlock lower boundary, ft
	float vgtmp13; // Sitka spruce zone upper boundary, ft in 10010206
	float vgtmp14; // Grand fir zone lower boundary, ft

}; // end of class PNVmodel

