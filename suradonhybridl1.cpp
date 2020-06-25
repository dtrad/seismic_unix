/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SURADON00FORW  $Date: Septem 2000  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "radon.h"
#include <signal.h>
#include <math.h>
#include <time.h>

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
  " nq=80           Size of Radon EOM                                   ",
  " testadj=0       =1 Test adjoint with random numbers                 ",
  " smooth=0        =1 Pass triangular filter to smooth  output         ",
  " eps1=1e-3       Numerator hyperparameter for Wm                     ",
  " verbose=0       =1  Extra Information                               ",
  " smute=2         stretch greater than smutex100 % is muted		",
  " nmofactor=1.9   nmofactor * equiv offset is used for nmo		",
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
  float eps; // small tolerance number
  int testadj; // =1 test adjoint
  int smooth; // =1 smoothing filter
  ///////////////////////////////////////////////////////////////
  int  i,k, ih, iq; // General counters 
  register int it;
  float *d;      /* single trace */
  float **data;   /* Common Scattering gather */
  float  *m;    /* Final modeled trace   */
  float  *t;     // time axis for input and output
  float  *h;      // halfoffset
  float t0;      // First useful time
  float *velint;/* array[nt] of vel for a particular trace */
  float *oa1t;	/* array[nt] of anis1 for a particular trace */
  float *oa2t;	/* array[nt] of anis2 for a particular trace */

  char *tmpdir;		/* directory path for tmp files		*/
  cwp_Bool istmpdir=cwp_false;/* true for user-given path		*/ 

  int verbose;		/* flag for echoing info		*/
  float fmax;
  int nq;  
  float *q;
  float dq;
  float **model;            // Temporal array for the model
  float step;
  float eps1;
  float eps2;
  int itercg;
  int iter_end;
  int norm;
  float *Wm;   // Factorization of inv model covariance matrix 
  float *Wd;   // Factorization of inv residual covariance matrix 
  float factor;
  float smute;
  float nmofactor;
  float quantil;
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
  int mute1;
  int mute2;
  float fmax1;
  float fmax2;
  float *ffilter;
  float *amps;

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
  fmax=0.8/(2*dt);
  nh=ntr;

  fprintf(stderr,"nt=%d,dt=%f\n",nt,dt);  

  // Get parameters 
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparfloat("eps", &eps))  eps = 1e-7;
  if (!getparfloat("eps1", &eps1))  eps1 = 5e-2;
  if (!getparfloat("eps2", &eps2))  eps2 = 5e-2;
  if (!getparfloat("step", &step))  step = 0.9;
  if (!getparint("itercg", &itercg))  itercg = 5;
  if (!getparint("iter_end", &iter_end))  iter_end = 1;
  if (!getparint("testadj", &testadj))  testadj = 0;
  if (!getparint("smooth", &smooth))  smooth = 0;
  if (!getparint("nq", &nq))  nq = 100;
  if (!getparint("norm",&norm)) norm=0;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=2;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
  if (!getparfloat("quantil",&quantil)) quantil=1;
 //////////////////////////////////////////////////////////////////////
  // The following are parameters use for radon_beam with two operators
  if (!getparfloat("depth1",&depth1)) depth1=1000;
  if (!getparfloat("depth2",&depth2)) depth2=1000;
  if (!getparfloat("qmin1",&qmin1)) qmin1=-1e-4;
  if (!getparfloat("qmin2",&qmin2)) qmin2=-5e-4;
  if (!getparint("nq1", &nq1))  nq1 = nq/2+20;
  if (!getparint("nq2", &nq2))  nq2 = nq-nq1;
  if (!getparfloat("factor1",&factor1)) factor1=4;
  if (!getparfloat("factor2",&factor2)) factor2=4;
  if (!getparint("rtmethod1", &rtmethod1))  rtmethod1 = 1;
  if (!getparint("rtmethod2", &rtmethod2))  rtmethod2 = 3;
  if (!getparint("mute1",&mute1)) mute1=0;
  if (!getparint("mute2",&mute2)) mute2=0;
  if (!getparfloat("fmax1",&fmax1)) fmax1=20;
  if (!getparfloat("fmax2",&fmax2)) fmax2=70;
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

  ///////////////////////////////////////////////////
  /* compute new square slowness and anis function */
#include "getvelocities.h" 
  velint=ealloc1float(nt);
  oa1t = ealloc1float(nt);
  oa2t = ealloc1float(nt);

  /* compute new square slowness and anis function */
  interpovv(nt,ncdp,cdp,ovv,oa1,oa2,cdpgather,velint,oa1t,oa2t);
  //radonhybrid(datah,nh,t,nt,dt,model,q,nq,velint,eps1,eps2,eps,mm,q,nq,itercg,iter_end,step,lsmethod,factor,smute,nmofactor,mute,parmute,he,nhcmp[ix],quantil,qmin1,qmin2,nq1,nq2,rtmethod1,rtmethod2, depth1,depth2, factor1, factor2, mute1, mute2, fmax1,fmax2,ffilter,amps);  
  radonhybrid(data,h,nh,t,nt,dt,model,q,nq,velint,itercg,iter_end,step,eps2,eps1,quantil,norm,qmin1,qmin2,nq1,nq2,factor1,factor2,smute,nmofactor,rtmethod1,rtmethod2,depth1,depth2,fmax1,fmax2,ffilter,amps);

  // Output trace
  erewind(tracefp);
  erewind(headerfp);  
    
  for (ih=0;ih<nh;ih++){ 
    efread(&tr,HDRBYTES,1,headerfp);
    memcpy((void *) tr.data,(const void *) data[ih],nt*sizeof(float));
    puttr(&tr);
  }

  if (1){
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
  free2float(oa1);
  free2float(oa2);
  free1float(oa1t);
  free1float(oa2t);
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






























