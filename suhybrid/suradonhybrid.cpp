/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SURADON00FORW  $Date: Septem 2000  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "radonhybrid.h"
#include <signal.h>
#include <math.h>
#include <time.h>
#include "Complex.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADONHYBRID - Frequency domain Radon Transform  using two         ",
  " 	            different operators.	       			",
  " 	   								",
  " suradfreq < stdin > stdout [optional parameters]          		",
  "                                                                     ",
  " Input must be sorted by offset.                                     ",
  "                                                                     ",
  " par=            file with stacking velocities and time as obtained  ",
  "                 from Velan (PARFILE)                                ",
  " nq=100          Size of Radon EOM                                   ",
  " nq1=nq/2+20     Size of model1                                      ",
  " nq2=nq-nq1      Size of model2                                      ",
  " smute=2         stretch greater than smutex100 % is muted		",
  " nmofactor=1.9   nmofactor * equiv offset is used for nmo		",
  " filter=0        =1 mutes model1, =2 mutes model2,                   ",
  "      	    =3 apply polygonal filter on 1                      ",
  "		    =4 apply polygonal filter on 2                      ",
  " rtmethod1= 1         						",
  " rtmethod2= 3         						",
  " fmax1 = 20        							",
  " fmax1 = 70        							",
  " depth1 = 1       For the pseudo offset of Foster and Mosher		",
  " depth2 = 3       Only used if rtmethod1 = 3 or rtmethod =3 		",
  " symmetric = 0    If = 1 defines q1 axis to be symmetric without	",
  "                  actually passing for zero                          ",
  " 	This program reads a data gather sorted by offset and produces  ",
  "     in the file = modelfile a Radon transform panel with two differents",
  "     model spaces, corresponding to two different Radon operators    ",
  "     The nature of each operator is defined by the corresponding     ",
  "     parameters rtmethod 1 and 2, depth 1 and 2, nq 1 and 2, etc.  	",
  "     The model obeys the following problem d=L1 m1 + L2 m2           ",
  "     The inverse RT gather (reconstructed) is output to stdoutput    ", 
  "     The filter parameter is used to mute the radon space before	",
  "     reconstruction.							",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset            */
/**************** end self doc ***********************************/
static void closefiles(void);
/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
FILE *modelfilep;       /* fp for model file output  		*/
FILE *offsetfile; 
///////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  segy tr;
  cwp_String modelfile=""; /* output sufile for the model */ 
  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;

  int nh;   // number of offset traces
  int nt;   // number of time samples 
  float dt; // Time sampling
  int testadj; // =1 test adjoint
  int smooth; // =1 smoothing filter
  ///////////////////////////////////////////////////////////////
  int  k, ih, iq; // General counters 
  register int it;
  float *d=0;      /* single trace */
  float **data=0;   /* Common Scattering gather */
  float  *m=0;    /* Final modeled trace   */
  float  *t=0;     // time axis for input and output
  float  *h=0;      // halfoffset
  float t0=0;      // First useful time
  float *velint=0;/* array[nt] of vel for a particular trace */
  float **ovv=0;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  float *cdp=0;	        /* array[ncdp] of cdps */  
  int ncdp;	/* number of cdps specified */
  char *tmpdir=0;		/* directory path for tmp files		*/

  cwp_Bool istmpdir=cwp_false;/* true for user-given path		*/ 

  int verbose;		/* flag for echoing info		*/

  int nq;  
  float *q;

  inv_par inv; // struct defined in inversion_par.h 
  float **model=0;            // Temporal array for the model
  float *Wm=0;   // Factorization of inv model covariance matrix 
  float *Wd=0;   // Factorization of inv residual covariance matrix 
  float smute;
  float nmofactor;
  int ntr;
  /// For radon_beam
  float qmin1;
  float qmin2;
  int nq1;
  int nq2;
  int rtmethod1;
  int rtmethod2;
  float depth1;
  float depth2;
  float factor1;
  float factor2;
  float fmax1;
  float fmax2;
  int filter;
  float *ffilter=0;
  float *amps=0;
  int symmetric;


  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
  // Register starting time
  start=time(0);
   // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  
  //if (!tr.offset) err("offset header field must be set");
  
  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  if (!(ntr=tr.ntr)) err("***ntr must be set\n");

  nh=ntr;

  fprintf(stderr,"nt=%d,dt=%f\n",nt,dt);  

  // Get parameters 
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparfloat("eps", &inv.eps))  inv.eps = 1e-7;
  if (!getparfloat("eps1", &inv.eps1))  inv.eps1 = 5e-2;
  if (!getparfloat("eps2", &inv.eps2))  inv.eps2 = 5e-2;
  if (!getparfloat("step", &inv.step))  inv.step = 0.9;
  if (!getparint("itercg", &inv.itercg))  inv.itercg = 5;
  if (!getparint("iter_end", &inv.iter_end))  inv.iter_end = 1;
  if (!getparint("norm",&inv.norm)) inv.norm=0;
  if (!getparint("testadj", &testadj))  testadj = 0;
  if (!getparint("smooth", &smooth))  smooth = 0;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=0;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";

  // The following are parameters use for radon_beam with two operators
  if (!getparfloat("depth1",&depth1)) depth1=1000;
  if (!getparfloat("depth2",&depth2)) depth2=1000;
  if (!getparfloat("qmin1",&qmin1)) qmin1=-1e-4;
  if (!getparfloat("qmin2",&qmin2)) qmin2=-5e-4;
  if (!getparint("nq1", &nq1))  nq1 = 100;
  if (!getparint("nq2", &nq2))  nq2 = 100;nq=nq1+nq2;
  if (!getparfloat("factor1",&factor1)) factor1=4;
  if (!getparfloat("factor2",&factor2)) factor2=4;
  if (!getparint("rtmethod1", &rtmethod1))  rtmethod1 = 1;
  if (!getparint("rtmethod2", &rtmethod2))  rtmethod2 = 2;
  if (!getparint("filter",&filter)) filter=2;
  if (!getparfloat("fmax1",&fmax1)) fmax1=30;
  if (!getparfloat("fmax2",&fmax2)) fmax2=70;
  if (!getparint("symmetric",&symmetric))  symmetric = 1;
  ////////////////////////////////////////////////////////////////////////
  /* Get frequencies that define the filter */

  int npoly;
  int namps;
  if ((npoly = countparval("ffilter"))!=0) {
    ffilter = ealloc1float(npoly);
    getparfloat("ffilter",ffilter);
  } else npoly = 0;

  if ((namps = countparval("amps"))!=0) {
    amps = ealloc1float(namps);
    getparfloat("amps",amps);
  } else  namps = 0;

  if (!(namps==npoly)) 	err("number of f values must = number of amps values");

  for (k=0;k<npoly;k++) fprintf(stderr,"ffilter[%d]=%f\n",k,ffilter[k]);
  for (k=0;k<namps;k++) fprintf(stderr,"amps[%d]=%f\n",k,amps[k]);

  ////////////////////////////////////////////////////////////////////////

  if ((modelfilep=fopen(modelfile,"w"))==NULL)
    err("cannot open file=%s\n",modelfile);

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
  
  // Allocate memory for data and model
  d=ealloc1float(nt);
  data=ealloc2float(nt,nh);
  model=alloc2float(nt,nq);
  m=ealloc1float(nt);
  t=ealloc1float(nt);
  q=ealloc1float(nq);
  h=ealloc1float(nh);
  Wm=alloc1float(nt*nq);
  Wd=alloc1float(nt*nh);

    
  // arrays for nmo model
  // Velint will contain velocity law for one cmp
  if ((velint=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Space for velint could not be allocated\n");

  /* Create axis for output cpds, equivalent offset and time */
  for(it=0;it<nt;it++) t[it]=0+it*dt;

  ///////////////////////////////////////////////////////////////////////
  
  memset( (void *) data[0], (int) '\0', nh * nt *FSIZE);
  ntr = 0;
  do {
    ntr++;
    if (verbose) fprintf(stderr,"tr.cdp=%d\n",tr.cdp);
    efwrite(&tr,HDRBYTES,1,headerfp);
    efwrite(tr.data,FSIZE, nt, tracefp);
  } while (gettr(&tr));
  erewind(tracefp);
  erewind(headerfp);  
  fprintf(stderr,"ntr=%d\n",ntr);
  ih=0;
  nh=ntr;
  
  for (ih=0;ih<nh;ih++){
    efread(tr.data, FSIZE, nt, tracefp);
    memcpy((void *) data[ih],(const void *) tr.data,nt*sizeof(float));
    efread(&tr,HDRBYTES,1,headerfp);
    h[ih]=tr.offset;
  }
  int cdpgather=tr.cdp;
  /* compute new square slowness and anis function */
  
  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  getvelocities(dt,nt,ncdp,cdp,ovv);
  /* compute new square slowness and anis function */
  interpovv(nt,ncdp,cdp,ovv,cdpgather,velint);
  
  radonhybrid(data,h,nh,t,nt,dt,model,q,nq,velint,inv,qmin1,qmin2,nq1,nq2,factor1,factor2,smute,nmofactor,rtmethod1,rtmethod2,depth1,depth2,fmax1,fmax2,filter,npoly,ffilter,amps,symmetric);

  // Output trace
  erewind(tracefp);
  erewind(headerfp);  
    
  for (ih=0;ih<nh;ih++){ 
    efread(&tr,HDRBYTES,1,headerfp);
    memcpy((void *) tr.data,(const void *) data[ih],nt*sizeof(float));
    puttr(&tr);
  }

  if (0){
    save_gather(model,nq,nt,0.004,"model2");
    system("sufft  < model2 | suamp | suxwigb title=\"model2\" perc=99 & ");
  }

  
  for (iq=0;iq<nq;iq++){
    memcpy((void *) tr.data,(const void *) model[iq],nt*sizeof(float));
    tr.dt=(int) (dt*1e6);       
    tr.ntr=nq;
    tr.ns=nt;
    tr.tracl=iq+1;
    tr.tracr=iq+1;
    tr.f2=q[iq];
    fputtr(modelfilep,&tr);    
  }

  free1float(cdp);
  free2float(ovv);
  free1float(Wd);
  free1float(Wm);
  free1float(q);
  free1float(velint);
  free1float(t);
  free1float(h);
  free1float(m);
  free2float(data);
  free2float(model);
  free1float(d);
  efclose(headerfp);
  efclose(tracefp);
  efclose(modelfilep);

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





























