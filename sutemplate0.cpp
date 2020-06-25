/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUTEMPLATE0:  $Date: September 2000  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include <signal.h>
#include <math.h>
#include <time.h>
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
 * Trace header fields accessed: ns, cdp, dt, offset             */
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
  segy tr,trf;
  register int it;
  int j, i, k, ih, ix; // General counters 
  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;

  int nh;   // number of offset traces
  int nt;   // number of time samples 
  int nq;    // number of model traces
  int nh2;   // number of output offset traces
  float dt; // Time sampling
  float dh; // offset interval
  float dx;  // midpoint interval
  unsigned int ntr;     // Total number of input traces
  int cdp; // cdp number
  ///////////////////////////////////////////////////////////////

  float **d;   /* Gather */
  float  **m;    /* Model   */
  float  *t;     // time axis for input and output
  float  *h;     // offset axis
  float *h2;     /* secondary offser axis for output offset */
  float t0;      // First useful time
  const double  pi=acos(-1.); 
  int verbose;		/* flag for echoing info		*/
  float fmax;
  char *tmpdir;		/* directory path for tmp files		*/
  cwp_Bool istmpdir=cwp_false;/* true for user-given path		*/ 

  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
  // Register starting time
  start=time(0);


  ///////////////////////////////////////////////////////////////////
  fprintf(stderr,"*****************\n");
  if (!getparint("verbose", &verbose)) verbose = 0; 
  if (!getparint("nq", &nq))  nq = 20;
  if (!getparfloat("t0",&t0)) t0=0;
  
  /* Look for user-supplied tmpdir */
  
  if (!getparstring("tmpdir",&tmpdir) &&
      !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
  if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
    err("you can't write in %s (or it doesn't exist)", tmpdir);  
 
    // Store traces in the disk for later use
  if (STREQ(tmpdir,"")){
    tracefp = etmpfile();
    headerfp = etmpfile();
    if (verbose) warn("using tmpfile() call");
  }else{ /* user-supplied tmpdir */
    char directory[BUFSIZ];
    strcpy(directory, tmpdir);
    strcpy(tracefile, temporary_filename(directory));
    strcpy(headerfile, temporary_filename(directory));
    /* Trap signals so can remove temp files */
    signal(SIGINT,  (void (*) (int)) closefiles);
    signal(SIGQUIT, (void (*) (int)) closefiles);
    signal(SIGHUP,  (void (*) (int)) closefiles);
    signal(SIGTERM, (void (*) (int)) closefiles);
    tracefp = efopen(tracefile, "w+");
    headerfp = efopen(headerfile, "w+");
    istmpdir=cwp_true;		
    if (verbose) warn("putting temporary files in %s", directory);
  }


  // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  if (!tr.ntr) err("ntr header field must be set");

  dt   = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  nh= (int) tr.ntr;

  fprintf(stderr,"nt=%d,dt=%f\n",nt,dt);  
  fprintf(stderr,"nh=%d\n",nh);
  
  
  ih=0;
  do {
    ih++;
    cdp=(int) tr.cdp;
    if (verbose) fprintf(stderr,"tr.cdp=%d\n",tr.cdp);
    efwrite(&tr,HDRBYTES,1,headerfp);
    efwrite(tr.data,FSIZE, nt, tracefp);
  } while (gettr(&tr));
  
  nh=ih;
  fprintf(stderr,"nh=%d\n",nh);

  erewind(tracefp);
  erewind(headerfp);


  d=ealloc2float(nt,nh);
  m=ealloc2float(nt,nq);


  for (ih=0;ih<nh;ih++){ 
    efread(&tr,HDRBYTES,1,headerfp);
    efread(tr.data,FSIZE,nt,tracefp);
    memcpy((void *) d[ih],(const void *) tr.data,nt*sizeof(float));    
    if (verbose) fprintf(stderr,"Original data ===> tr.offset=%d\n",tr.offset);  
    //tr.ntr=nhcmp[ix];
    puttr(&tr);
  }

  erewind(tracefp);
  erewind(headerfp);


  for (ih=0;ih<nh;ih++){ 
    efread(&tr,HDRBYTES,1,headerfp);
    memcpy((void *) tr.data,(const void *) d[ih],nt*sizeof(float));
    //if (verbose) fprintf(stderr,"Original data ===> tr.offset=%d\n",tr.offset);  
    //tr.ntr=nhcmp[ix];
    puttr(&tr);
  }

  if (istmpdir) eremove(headerfile);
  if (istmpdir) eremove(tracefile);  

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  free2float(d);
  free2float(m);


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






























