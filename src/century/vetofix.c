
/*****************************************************************
        $RCSfile: vetofix.c,v $
        $Revision: 1.1 $
        $Date: 2006/10/13 18:22:26 $
        $Locker: wellsj $

$Log: vetofix.c,v $
Revision 1.1  2006/10/13 18:22:26  wellsj
Initial revision

Revision 1.3  2005/03/26 18:11:34  lenihan
reassigned vclass types to different fixed types

Revision 1.2  2004/09/24 19:54:16  bachelet
hedaer

******************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>


#include "category_bgc.h"

/* VCategory	*vclass; */

void
vetofix_(vclass, fixname)
  int     *vclass;
	char		*fixname;
{

VCategory vegclass;

vegclass = (VCategory)(*vclass);


	switch (vegclass) {
/** D. Bachelet 9-25-97 Need Ice to do the world  **/
                case VIce :
            (void) strcpy(fixname, FixIce);
                        break;
		case VTundra :
            (void) strcpy(fixname, FixTundra);
			break;
                case VTaiga :
            (void) strcpy(fixname, FixTundra);
                        break;
                case VBorealLarchForest :
            (void) strcpy(fixname, FixBoreal);
                        break;
		case VBorealConiferousForest :
            (void) strcpy(fixname, FixBoreal);
			break;
		case VMaritimeTemperateConiferousForest :
            (void) strcpy(fixname, FixForest);
			break;
		case VContinentalTemperateConiferousForest :
            (void) strcpy(fixname, FixForest);
			break;
		case VCoolTemperateMixedForest :
            (void) strcpy(fixname, FixForest);
			break;
		case VWarmTemperateSubtropicalMixedForest :
            (void) strcpy(fixname, FixForest);
			break;
		case VTemperateDeciduousForest :
            (void) strcpy(fixname, FixForest);
			break;
		case VTropicalDeciduousForest :
            (void) strcpy(fixname, FixDryTropical);
			break;
		case VTropicalEvergreenForest :
            (void) strcpy(fixname, FixTropical);
			break;
		case VTemperateMixedXeromorphicWoodland :
            (void) strcpy(fixname, FixDryGrass);
			break;
		case VTemperateConiferXeromorphicWoodland :
            (void) strcpy(fixname, FixDryGrass);
			break;
		case VTropicalThornWoodland :
            (void) strcpy(fixname, FixDryTropical);
			break;
		case VTemperateSubtropicalDeciduousSavanna :
            (void) strcpy(fixname, FixDryForest);
			break;
		case VWarmTemperateSubtropicalMixedSavanna :
            (void) strcpy(fixname, FixDryForest);
			break;
		case VTemperateConiferSavanna :
            (void) strcpy(fixname, FixDryGrass);
            		break;
		case VTropicalDeciduousSavanna :
            (void) strcpy(fixname, FixDryTropical);
			break;
		case VC3Grasslands :
            (void) strcpy(fixname, FixGrassland);
			break;
		case VC4Grasslands :
            (void) strcpy(fixname, FixGrassland);
			break;
		case VMediterraneanShrubland :
            (void) strcpy(fixname, FixDryGrass);
			break;
		case VTemperateAridShrubland :
            (void) strcpy(fixname, FixDryGrass);
			break;
		case VSubtropicalAridShrubland :
            (void) strcpy(fixname, FixDryTropical);
			break;
/*
		case VTemperateBroadleafEvergreenForest :
            (void) strcpy(fixname, FixForest);
			break;   */
        default:
            /**********
                Garbage in, garbage out.
            ***********/
            (void) strcpy(fixname, "BOGUS");

            (void) fprintf(stderr,
                "Hi from vetofix.c bogometer reading found,\n"
                "\tunknown VEMAP veg type entered VEMAP -> fixname translator.\n"
                "\tLook for math error problems.\n"
				"vegclass = %d\n", vegclass);
			assert(0);
            break;
	}
}

