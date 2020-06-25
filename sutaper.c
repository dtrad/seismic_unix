/* Copyright (c) Colorado School of Mines, 2000.*/
/* All rights reserved.                       */

/* SUTAPER: $Revision: 1.10 $ ; $Date: 1996/09/13 21:49:05 $	*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" SUTAPER - Taper the edge traces of a data panel to zero.	",
" 								",
" sutaper <diskfile >stdout [ntaper=5]				",
" 								",
" The taper is an \"ntaper\" point sine-squared taper 		",
" symmetrically applied at each end of the data set.		",
" 								",
NULL};

/* Credits:
 *
 *	CWP: Chris, Jack
 *
 * Trace header fields accessed: ns
 */
/**************** end self doc ***********************************/


#define NTAPER	5

segy tr;

int
main(int argc, char **argv)
{
	int nt;		/* number of sample points on traces	*/
	int itr, ntr;	/* trace counter and total 		*/
	float *taper;	/* vector of taper weights		*/
	int ntaper;	/* number of taper weights		*/


	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);


	/* Get info from first trace */
	if (!(ntr = gettra(&tr, 0)))  err("can't get first trace");
	nt = (int) tr.ns;	/* disaster to pass &ushort */


	/* Get parameter */
	if (!getparint("ntaper", &ntaper))	ntaper = NTAPER;
	if (ntaper > ntr/2)
	    err("taper vector length, %d, exceeds ntr/2 = %d", ntaper, ntr/2);


	/* Set up taper weights */
	taper = ealloc1float(ntaper);
	{ register int k;
	  for (k = 0; k < ntaper; ++k) {
		float s = sin(k*PI/(2*ntaper));
		taper[k] = s*s;
	  }
	}
						


	/* Loop over the traces, tapering those at the ends */

	/* Take care of the trace already read (taper[0] = 0.0) */
	{ register int i;
	  for (i = 0; i < nt; ++i)  tr.data[i] = 0.0;
	}
	puttr(&tr);

	/* Taper at the left end of the data set */
	for (itr = 1; itr < ntaper; ++itr) {
		register int i;

		gettr(&tr);
		for (i = 0; i < nt; ++i)  tr.data[i] *= taper[itr];
		puttr(&tr);
	}

	/* Pass on traces in the middle of the data set */
	for ( ; itr < ntr - ntaper; ++itr) {
		gettr(&tr);
		puttr(&tr);
	}
	
	/* Taper at the right end of the data set */
	for ( ; itr < ntr; ++itr) {
		register int i;

		gettr(&tr);
		for (i = 0; i < nt; ++i)  tr.data[i] *= taper[ntr - itr - 1];
		puttr(&tr);
	}


	return EXIT_SUCCESS;
}
