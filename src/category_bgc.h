
/*****************************************************************
  $RCSfile: category_bgc.h,v $
  $Revision: 1.1 $
  $Date: 2006/10/13 18:38:31 $
  $Locker: wellsj $

  $Log: category_bgc.h,v $
  Revision 1.1  2006/10/13 18:38:31  wellsj
  Initial revision

  Revision 1.2  2004/09/24 17:54:12  bachelet
  header

 ******************************************************************/

/**********
  These are all the VEMAP veg types.  These are the vveg.v1 type
  names taken straight from the email from Nan Rosenbloom on Mon, 19
  Dec 94 16:18:08 MST.
 ***********/
typedef enum {
	VUnknown				= '\0',

	VTundra					= 1,
	VBorealConiferousForest			= 2,
	VMaritimeTemperateConiferousForest	= 3,
	VContinentalTemperateConiferousForest	= 4,
	VCoolTemperateMixedForest		= 5,
	VWarmTemperateSubtropicalMixedForest	= 6,
	VTemperateDeciduousForest		= 7,
	VTropicalDeciduousForest		= 8,
	VTropicalEvergreenForest		= 9,
	VTemperateMixedXeromorphicWoodland	= 10,
	VTemperateConiferXeromorphicWoodland	= 11,
	VTropicalThornWoodland			= 12,
	VTemperateSubtropicalDeciduousSavanna	= 13,
	VWarmTemperateSubtropicalMixedSavanna	= 14,
	VTemperateConiferSavanna		= 15,
	VTropicalDeciduousSavanna		= 16,
	VC3Grasslands				= 17,
	VC4Grasslands				= 18,
	VMediterraneanShrubland			= 19,
	VTemperateAridShrubland			= 20,
	VSubtropicalAridShrubland		= 21,
	VTaiga					= 22,
	VBorealLarchForest			= 23,
	VIce					= 90,
	VInlandWaterBodies			= 91,
	VWetlands				= 92
} VCategory;


/** D. Bachelet 9-25-97 Need Ice when doing the world **/

/** Defines for fix.100 file names **/

/**Although fix.100 file names are dimensioned 15, you will
  see only 14 spaces here.  15th is reserved for null terminator.
  Important when doing strcpy in vetofix. **/

#  define FixTundra             "arcfix.100    "
#  define FixBoreal             "borfix.100    "
#  define FixDryForest          "dryffix.100   "
#  define FixDryGrass           "drygfix.100   "
#  define FixDryTropical        "drytrpfix.100 "
#  define FixForest             "ffix.100      "
#  define FixGrassland          "gfix.100      "
#  define FixTropical           "trpfix.100    "
#  define FixIce		"arcfix.100    "

typedef struct { char * paramName; float paramVals[5]; } CenParamStruct;

typedef struct { int first_calendar_year; float yearlyCO2conc[5000]; } CO2Struct;


