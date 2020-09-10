/*
 *  flowUtilities.c
 *  MC1
 *
 *  Created by Dave Conklin on 4/12/10.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


/*	from old flowd.c file
 *	Global data file for flow routines in century.
 */

#define LENFST 500	/* Length of flowstack */

struct stack{
	float *from;	/* Source */
	float *to;	/* Destination */
	float when;	/* Time to flow */
	float amt;	/* Amount */
} flowstack[LENFST+1];

int nflows;		/* Number of flows */

void flow_err(int error_num, float when);

/*	from old flow.c file
 *      Schedules a flow by storing the arguments in a stack which
 *      will be accessed by the function 'flowup' which performs
 *      the actual flow for all stack elements where 'when' is
 *      less than or equal to time.
 *
 *	Called from the Fortran version, the underscore in flow_
 *	allows for correct referencing during linkage.  No 
 *	wrapper - this routine was written to be called directly
 *	from Fortran.
 */

	void
flow_(float *from, float *to, float *when, float *howmuch)
{
	if (isnan(*from) || isnan(*to) || isnan(*howmuch))
	{
		printf("flow_(%f, %f, %f, %f)\n", *from, *to, *when, *howmuch);
		assert(0);
	}

	/* Increment the number of flows stored in the stack */

	nflows += 1;

	if (nflows > LENFST)
	{
		/* Stack Overflow */
		flow_err(1, *when);
		exit(1);
	}


	else
		/* Store the arguments in the stack */
	{
		flowstack[nflows].from = from;
		flowstack[nflows].to = to;
		flowstack[nflows].when = *when;
		flowstack[nflows].amt = *howmuch;
	}

	return;
}


/*	from old flowup.c file
 *	Perform the flows for all stack elements whose 'when'
 *	value is less than or equal to time.
 *
 *	Called from the Fortran version, the underscore in flowup_
 *	allows for correct referencing during linkage.  No 
 *	wrapper - this routine was written to be called directly
 *	from Fortran.
 */

	void
flowup_(float *time)
{
	int FlowsDone, FlowsNotDone, i;

	FlowsNotDone = 0;
	FlowsDone = 0;

	if (nflows <= 0.0)
		return;

	/* If there are any flows in the stack, determine which
	 * need to go now and do it.
	 */
	else
		for (i=1; i<=nflows; i++)
		{
			if (*time < flowstack[i].when)
			{
				FlowsNotDone +=1;
				/* This one doesn't need to be done yet; 
				 * move it down the stack if other flows
				 * have been performed already.
				 */
				if (FlowsDone > 0)
					flowstack[FlowsNotDone] = flowstack[i];
			}
			else
			{
				if (flowstack[i].amt != 0.0)
				{
					*(flowstack[i].from) -= flowstack[i].amt;
					*(flowstack[i].to) += flowstack[i].amt;
				}
				FlowsDone += 1;
			}

		}

	nflows = FlowsNotDone;
	return;
}

// from old file ferr.c
/*
 *	error messages
 */

static char *err_msg[] = {
	"No Error.",					/* 0 */
	"Stack Overflow",				/* 1 */
	"Value of source variable is negative",		/* 2 */
	"Amount to flow has been limited",		/* 3 */
};

	void
flow_err(int error_num, float when)
{
	fprintf(stderr,
			"\n Century flow error occurred at time = %10.4f\n%s\n\n",
			when, err_msg[error_num]);
	return;
}
