// MAPSSvegClasses.h

#ifndef MAPSSVEGCLASSES_H
# define MAPSSVEGCLASSES_H

// MAPSS veg types
typedef enum {

	MissingSoilsData	= -7,					/* something missing */
	MisClass			= -5,					/* could not find class */
	Unknown				= '\0',					/* Ocean */

	NaturalBarren        = 15,

	ForestDeciduousBroadleaf				= 100,	/* Forest Rules */
	ForestMixedWarm_DEB						= 101,
	ForestMixedCool							= 102,
	ForestDeciduousNeedle					= 103,	
	ForestEvergreenBroadleafTropical		= 105,
	ForestEvergreenNeedleTaiga				= 107,
	ForestMixedWarm_EN						= 108,
	ForestHardwoodCool						= 111,		

	ForestEvergreenNeedleMaritime			= 112,
	ForestEvergreenNeedleContinental		= 113,
	ForestSubalpine      = 114,

	TreeSavannaDeciduousBroadleaf 			= 200,	/* TreeSavanna Rules */
	TreeSavannaMixedWarm_DEB				= 201,
	/* TreeSavannaDeciduousNeedle				= 202,		 not used */
	ForestSeasonalTropical_ED				= 109,
	ForestSavannaDryTropical_ED				= 110,
	TreeSavannaMixedCool_EN					= 205,
	TreeSavannaMixedWarm_EN					= 206,

	TreeSavannaEvergreenNeedleMaritime		= 207,
	TreeSavannaEvergreenNeedleContinental	= 208,

	TreeSavannaPJContinental				= 209,	
	TreeSavannaPJMaritime					= 210,
	TreeSavannaPJXericContinental			= 211,
	TreeSavannaSubalpine = 214,

	OpenShrublandNoGrass					= 301,

	ShrubSavannaDeciduousBroadleaf			= 302,	/* Shrub Rules */
	ShrubSavannaMixedWarm_DEB				= 303,
	/* ShrubSavannaDeciduousMicro				= 304,		 not used */
	ShrubSavannaTropical_EB					= 305,
	ShrubSavannaMixedCool_EN				= 307,
	ShrubSavannaEvergreenMicro				= 308,

	ShrubSavannaSubTropicalMixed			= 310,	
	ShrublandSubTropicalXeromorphic			= 311,	
	ShrublandSubTropicalMediterranean		= 312,	

	/* Grassland Rules */
	GrassTallC3								= 414,
	GrassMidC3								= 415,
	GrassShortC3							= 416,

	GrassTallC3C4							= 417,
	GrassMidC3C4							= 418,
	GrassShortC3C4							= 419,

	GrassTallC4								= 420,
	GrassMidC4								= 421,
	GrassShortC4							= 422,

	GrassSemiDesertC3						= 423,
	GrassSemiDesertC3C4						= 424,
	GrassSemiDesertC4						= 425,


	DesertBoreal							= 500,	/* Desert Rules */
	DesertTemperate							= 501,
	DesertSubtropical						= 502,
	DesertTropical							= 503,
	DesertExtreme							= 504,

	TaigaTundra								= 600,	/* Heat Limited */
	Tundra 									= 601,
	Ice										= 602,

} MAPSSvegClass; // end of typedef enum MAPSSvegClass


// aggregated veg types from Bachelet et al. 2001 Ecology 4: 164-185
typedef enum  
{
	aggvUnknown = 0,
	aggvTundra = 1,
	aggvTaigaTundra = 2,
	aggvConiferForest = 3,
	aggvNEmixedForest = 4,
	aggvTemperateDeciduousForest = 5,
	aggvSEmixedForest = 6,
	aggvTropicalBroadleafForest = 7,
	aggvSavannaWoodland = 8,
	aggvShrubWoodland = 9,
	aggvGrassland = 10,
	aggvAridLand = 11
}  aggVegClass; // end of typedef enum aggVegClass



#endif /* MAPSSVEGCLASSES_H */
