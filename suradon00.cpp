/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SURADON0:  $Date: September 2000  */

#include "su.h"
#include "segy.h"
#include "clibrarytd.h"
#include "header.h"
#include <signal.h>
#include <math.h>
#include <time.h>
#include "eomig.h"
#include "clibrary.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADFREQ - Frequency domain Parabolic Radon Transform              ",
  "             (for interpolation)                                     ",
  " 	   								",
  " suradfreq < stdin > stdout [optional parameters]          		",
  "                                                                     ",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset
/**************** end self doc ***********************************/
static void closefiles(void);
/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
FILE *offsetfile; 
///////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  segy tr,trf,traux;
  register int it;
  int j, i, k, ih, ix; // General counters 
  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;

  int nh;   // number of offset traces
  int nt;   // number of time samples 
  int nx;    // number of midpoints
  int nh2;   // number of output offset traces
  float dt; // Time sampling
  float dh; // offset interval
  float dx;  // midpoint interval

  ///////////////////////////////////////////////////////////////

  float *d;      /* single trace */
  float **d;   /* Common Scattering gather */
  float  **m;    /* Final migrated trace   */
  float  *t;     // time axis for input and output
  float  *h;      // halfoffset
  float *h2;     /* axis for output offset if resample==1 */
  float t0;      // First useful time
  const double  pi=acos(-1.); 
  int verbose;		/* flag for echoing info		*/
  float fmax;





  if (istmpdir) eremove(headerfile);
  if (istmpdir) eremove(tracefile);  

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}


/* for graceful interrupt termination */
static void closefiles(void)
{
	efclose(headerfp);
	efclose(tracefp);
	eremove(headerfile);
	eremove(tracefile);
	exit(EXIT_FAILURE);
}






























