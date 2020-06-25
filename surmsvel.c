/* Copyright (c) Colorado School of Mines, 1999.*/
/* All rights reserved.                       */

/* SUINTVEL: $Revision: 1.8 $ ; $Date: 1996/09/13 21:49:05 $		*/

#include "su.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"									",
" SURMSVEL - convert interval velocity model to stacking velocity model	",
"									",
" surmsvel v=  h=  outpar=/dev/tty					",
"									",
" Required parameters:					        	",
"				h=layer thicknesses vector		",
"				v=interval velocities vector		",
"									",
" Optional parameters:							",
"	outpar=/dev/tty		output parameter file in the form:	",
"	v =	stacking velocities 					",
"	h =	normal incidence times		 			",
"									",
" Examples:								",
"    surmsvel vs=5000,5523,6339,7264 t0=.4,.8,1.125,1.425 outpar=intpar	",
"									",
"    surmsvel par=intpar outpar=stkpar					",
"									",
" Note: suintvel does not have standard su syntax since it does not	",
"      operate on seismic data.  Hence stdin and stdout are not used.	",
"									",
" Note: may go away in favor of par program, velconv, by Dave		",
"									",
NULL};

/* Credits:
 *	CWP: Jack 
 *
 * Technical Reference:
 *	The Common Depth Point Stack
 *	William A. Schneider
 *	Proc. IEEE, v. 72, n. 10, p. 1238-1254
 *	1984
 *
 * Formulas:
 *    	Note: All sums on i are from 1 to k
 *
 *	From Schneider:
 *	Let h[i] be the ith layer thickness measured at the cmp and
 *	v[i] the ith interval velocity.
 *	Set:
 *		t[i] = h[i]/v[i]
 *	Define:
 *		t0by2[k] = 0.5 * t0[k] = Sum h[i]/v[i]
 *		vh[k] = vs[k]*vs[k]*t0by2[k] = Sum v[i]*h[i]
 *	Then:
 *		dt[i] = h[i]/v[i] = t0by2[i] - t0by2[i-1]
 *		dvh[i] = h[i]*v[i] = vh[i] - vh[i-1]
 *		h[i] = sqrt(dvh[i] * dt[i])
 *		v[i] = sqrt(dvh[i] / dt[i])
 *
 *
 */
/**************** end self doc *******************************************/

int
main(int argc, char **argv)
{
	float *v=NULL;		/* interval velocities		*/
	float *h=NULL;		/* layer thicknesses at the cmp	*/
	float *vs=NULL;	/* stacking velocities		*/
	float *t=NULL;	/* zero incidence times		*/
	float *dt=NULL;	/* zero incidence times		*/
	register int i;		/* counter				*/
	int n;			/* number of layers			*/
	float ttot;
	float sum;
	cwp_String outpar;	/* name of file holding output parfile	*/
	FILE *outparfp;		/* ... its file pointer			*/
	

	/* Initialize */
	initargs(argc, argv);
	requestdoc(0);


	outpar = "/dev/tty" ;	getparstring("outpar", &outpar);
	outparfp = efopen(outpar, "w");


	/* Allocate space for the model */
	if ((n = countparval("h"))) {
		v  = ealloc1float(n);
		h  = ealloc1float(n);
		dt  = ealloc1float(n);
		vs = ealloc1float(n);
		t = ealloc1float(n);
	} else err("no h's specified");

	/* Get the normal incidence times and stacking velocities */
	if (n != getparfloat("h", h))
		err("expected %d intervals", n);
	if (n != getparfloat("v", v))
		err("expected %d velocities", n);

	/* Check that vs's and t0's are positive */
	for (i = 0; i < n; i++) {
		if (v[i] <= 0.0)
			err("v's must be positive: v[%d] = %f", i, v[i]);
		if (h[i] <= 0.0)
			err("h's must be positive: h[%d] = %f", i, h[i]);
		dt[i]=h[i]/v[i];
	}

	/* Compute h(i), v(i) */
	vs[0] = v[0];
	ttot = dt[0];
	t[0]=ttot;
	sum=v[0]*v[0]*dt[0];
	for (i = 1; i < n; i++) {
	  ttot=ttot+dt[i];
	  sum=sum+v[i]*v[i]*dt[i];
	  vs[i]=sqrt(sum/ttot);
	  t[i]=ttot;
	}

	/* Make par file */
	fprintf(outparfp, "t=");
	for (i = 0; i < n - 1; i++) {
		fprintf(outparfp, "%g,", t[i]);
	}
	fprintf(outparfp, "%g\n", t[n-1]);

	fprintf(outparfp, "vs=");
	for (i = 0; i < n - 1; i++) {
		fprintf(outparfp, "%g,", vs[i]);
	}
	fprintf(outparfp, "%g\n", vs[n-1]);


	return EXIT_SUCCESS;
}



