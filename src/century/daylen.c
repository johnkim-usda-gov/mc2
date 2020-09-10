
/*****************************************************************
  $RCSfile: daylen.c,v $
  $Revision: 1.1 $
  $Date: 2006/10/13 18:22:26 $
  $Locker: wellsj $

  $Log: daylen.c,v $
  Revision 1.1  2006/10/13 18:22:26  wellsj
  Initial revision

  Revision 1.2  2004/09/24 19:49:05  bachelet
  header

 ******************************************************************/


/*              Copyright 1993 Colorado State University     */
/*                      All Rights Reserved                  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/* #define M_PI 3.1415926536 */

	float
daylen_(int *month, float *sitlat)
{

	/*  Author: Melannie Hartman	*/
	/*  Changed 1/96 for inclusion in Gridded Century  -*/
	/*  Compute daylength		*/
	/*  rlat:	latitude in radians	*/

	double adelt;
	double ahou;
	float  dylngth;	/* daylength in hours	*/
	double rlat;	/* latitude in radians	*/
	double temp1, temp2;
	double jday[] = {	1.	/* Jan 1 julian day	*/,
		32.	/* Feb 1 julian day	*/,
		61.	/* Mar 1 julian day	*/,
		92.	/* Apr 1 julian day	*/,
		122.	/* May 1 julian day	*/,
		153.	/* Jun 1 julian day	*/,
		183.	/* Jul 1 julian day 	*/,
		214.	/* Aug 1 julian day	*/,
		245.	/* Sep 1 julian day	*/,
		275.	/* Oct 1 julian day	*/,
		306.	/* Nov 1 julian day	*/,
		337.	/* Dec 1 julian day	*/};

	/* NOTE: JDAY isn't exact, but should be close enough. */

	/* Convert sitlat which is in degrees to radians
	   rlat = ((double)*sitlat * M_PI / 180.0);     */
	rlat = ((double)*sitlat * 3.1415926536 / 180.0);

	/*    temp1 = 2.0 * M_PI * (jday[*month] - 77.0) / 365.0;  */
	temp1 = 2.0 * 3.1415926536 * (jday[*month] - 77.0) / 365.0;

	/*    fprintf(stderr, "temp1: %lf  month: %d \n",temp1,*month); */
	adelt = 0.4014 * sin((double)temp1);
	temp1 = 1.0 - pow((-tan(rlat) * (adelt)), 2.0);
	if (temp1 < 0.0)
	{
		/**	printf("error in daylen - temp1 = %lf\n", temp1);
		  printf("rlat = %lf\n", rlat);
		  printf("month = %d\n", month);
		  printf("julian day = %lf\n", jday[*month]); **/
		temp1 = 0.0;
	}
	temp1 = sqrt(temp1);
	temp2 = -tan(rlat) * tan(adelt);
	ahou = atan2(temp1,temp2);


	/*     dylngth = (float)((ahou / M_PI) * 24.0);*/
	dylngth = (float)((ahou / 3.1415926536) * 24.0);
	return dylngth;
}



