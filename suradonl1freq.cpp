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
#include "dan.h"
#include "radonl1freq.h"

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
FILE *modelfilep;       /* fp for model file output  		*/
FILE *offsetfile; 
///////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{
  segy tr;
  op_param op2; /* struct defined in radonl1freq.h */
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
  int  k, ih, iq; // General counters 
  register int it;
  float *d;      /* single trace */
  float **data;   /* Common Scattering gather */
  float  *m;    /* Final modeled trace   */
  float  *t;     // time axis for input and output
  float  *h;      // halfoffset
  float t0;      // First useful time

  int verbose;		/* flag for echoing info		*/
  float fmax;
  int nq;  
  float *q;
  float **model;            // Temporal array for the model
  float eps1;
  float eps2;
  int itercg;
  int iter_end;
  int norm;
  float smute;
  float nmofactor;
  int ntr;
  /// For radon_beam
  float *ffilter=0;
  float *amps=0;
   /* Velocity */
  float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  float *cdp;	        /* array[ncdp] of cdps */  
  int ncdp;	/* number of cdps specified */
  float *velint;/* array[nt] of vel for a particular trace */


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

  fprintf(stderr,"nt=%d,dt=%f,nh=%d\n",nt,dt,nh);  

  // Get parameters 
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparfloat("eps", &eps))  eps = 1e-7;
  if (!getparfloat("eps1", &eps1))  eps1 = 5e-2;
  if (!getparfloat("eps2", &eps2))  eps2 = 5e-2;
  if (!getparint("itercg", &itercg))  itercg = 5;
  if (!getparint("iter_end", &iter_end))  iter_end = 1;
  if (!getparint("testadj", &testadj))  testadj = 0;
  if (!getparint("smooth", &smooth))  smooth = 0;
  if (!getparint("norm",&norm)) norm=0;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=2;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
 //////////////////////////////////////////////////////////////////////
  // The following are parameters use for radon_beam with two operators
  if (!getparfloat("depth1",&op2.depth1)) op2.depth1=1000;
  if (!getparfloat("depth2",&op2.depth2)) op2.depth2=1000;
  if (!getparfloat("qmin1",&op2.qmin1)) op2.qmin1=-1e-4;
  if (!getparfloat("qmin2",&op2.qmin2)) op2.qmin2=-5e-4;
  if (!getparint("nq1", &op2.nq1))  op2.nq1 = 50;
  if (!getparint("nq2", &op2.nq2))  op2.nq2 = 50;
  if (!getparfloat("factor1",&op2.factor1)) op2.factor1=4;
  if (!getparfloat("factor2",&op2.factor2)) op2.factor2=4;
  if (!getparint("rtmethod1", &op2.rtmethod1))  op2.rtmethod1 = 1;
  if (!getparint("rtmethod2", &op2.rtmethod2))  op2.rtmethod2 = 3;
  if (!getparint("mute1",&op2.mute1)) op2.mute1=0;
  if (!getparint("mute2",&op2.mute2)) op2.mute2=0;
  if (!getparfloat("fmax1",&op2.fmax1)) op2.fmax1=20;
  if (!getparfloat("fmax2",&op2.fmax2)) op2.fmax2=70;
  ////////////////////////////////////////////////////////////////////////
  /* Get frequencies that define the filter */
  nq=op2.nq1+op2.nq2;

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

  
  // Allocate memory for data and model
  d=ealloc1float(nt);
  data=ealloc2float(nt,nh);
  model=alloc2float(nt,nq);
  m=ealloc1float(nt);
  t=ealloc1float(nt);
  q=ealloc1float(nq);
  h=ealloc1float(nh);

  TRACE;
  memset( (void *) data[0], (int) '\0', nh * nt *FSIZE);
  ih = 0;
  do {
    if (1) fprintf(stderr,"tr.cdp=%d,ih=%d\n",tr.cdp,ih);
    memcpy((void *) data[ih],(const void *) tr.data,nt*sizeof(float));
    h[ih]=tr.offset;
    ih++;
  } while (gettr(&tr));
  TRACE;
  erewind(stdin);
  nh=ih;
  if (ntr!=nh) fprintf(stderr,"ntr=%d is different from nh=%d \n",ntr,nh);
  if (nh>ntr) err("Change ntr in header\n");
  int cdpgather=tr.cdp;
  TRACE;
  ///////////////////////////////////////////////////
  /* compute new square slowness and anis function */
  //#include "getvelocities.h" 
  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  getvelocities(dt,nt,ncdp,cdp,ovv);
  /* compute new square slowness and anis function */
  interpovv(nt,ncdp,cdp,ovv,cdpgather,velint);

  /* Create axis for output cpds, equivalent offset and time */
  for(it=0;it<nt;it++) t[it]=0+it*dt;

  radonl1freq(data,h,nh,t,nt,dt,model,q,velint,itercg,iter_end,op2,smute,
	      nmofactor,ffilter,amps);

  // Output trace
    
  for (ih=0;ih<nh;ih++){
    gettr(&tr);
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
  free1float(velint);
  free1float(q);
  free1float(t);
  free1float(h);
  free1float(m);
  free2float(data);
  free2float(model);
  free1float(d);
  efclose(modelfilep);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}






























