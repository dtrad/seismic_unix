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

int testfunction();

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADFREQ - Frequency domain Parabolic Radon Transform              ",
  " 	   								",
  " suradfreq < stdin > stdout [optional parameters]          		",
  "                                                                     ",
  " Input must be sort by cdp and offset, for example                   ",
  "                   susort cdp offset < input | .....                 ",
  "                                                                     ",
  " par=            file with stacking velocities and time as obtained  ",
  "                 from Velan (PARFILE)                                ",
  " offsetfile      ascii file that contains the desired output offset  ",
  "                 IMPORTANT: nhe has to be larger than the number of  ",
  "                 new offsets.                                        ",
  "                                                                     ", 
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
/* Globals (so can trap signal) defining temporary disk files */
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
  int  i, ih, iq; // General counters 
  register int it;
  float *d;      /* single trace */
  float **data;   /* Common Scattering gather */
  float  *m;    /* Final modeled trace   */
  float  *t;     // time axis for input and output
  float  *h;      // halfoffset
  float t0;      // First useful time

  /* Velocity */
  float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  float *cdp;	        /* array[ncdp] of cdps */  
  int ncdp;	/* number of cdps specified */
  float *velint;/* array[nt] of vel for a particular trace */

  int verbose;		/* flag for echoing info		*/
  float fmax;
  int nq;  
  float qmin;
  float qmax;
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
  float depth;   // depth for pseudohyperbolic
  int ntr;
  int rtmethod;

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
  if (!getparfloat("eps", &eps))  eps = 1e-7;
  if (!getparfloat("eps1", &eps1))  eps1 = 5e-2;
  if (!getparfloat("eps2", &eps2))  eps2 = 5e-2;
  if (!getparfloat("step", &step))  step = 0.9;
  if (!getparint("itercg", &itercg))  itercg = 5;
  if (!getparint("iter_end", &iter_end))  iter_end = 1;
  if (!getparint("testadj", &testadj))  testadj = 0;
  if (!getparint("smooth", &smooth))  smooth = 0;
  if (!getparfloat("qmin", &qmin))  qmin = -5e-9;
  if (!getparfloat("qmax", &qmax))  qmax = 5e-9;
  if (!getparint("nq", &nq))  nq = 100;
  if (!getparint("norm",&norm)) norm=0;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparfloat("factor",&factor)) factor=0.8;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=2;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
  if (!getparfloat("quantil",&quantil)) quantil=1;
  if (!getparint("rtmethod",&rtmethod)) rtmethod=2;
  if (!getparfloat("depth",&depth)) depth=2000;
  if (!getparfloat("fmax",&fmax)) fmax=0.8/(2*dt);

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
  Wm=alloc1float(nt*nq);
  Wd=alloc1float(nt*nh);

  // Velocity axis for Radon space

  dq = (qmax-qmin)/(nq-1);
  for (i=0;i<nq;i++) q[i] = qmin+i*dq;

  // arrays for nmo model
  // Velint will contain velocity law for one cmp
  if ((velint=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Space for velint could not be allocated\n");

  /* Create axis for output cpds, equivalent offset and time */
  for(it=0;it<nt;it++) t[it]=0+it*dt;

  ///////////////////////////////////////////////////////////////////////
  
  memset( (void *) data[0], (int) '\0', nh * nt *FSIZE);
  ih = 0;
  do {
    if (verbose) fprintf(stderr,"tr.cdp=%d\n",tr.cdp);
    memcpy((void *) data[ih],(const void *) tr.data,nt*sizeof(float));
    h[ih]=tr.offset;
    ih++;
  } while (gettr(&tr));
  erewind(stdin);
  nh=ih;
  if (ntr!=nh) fprintf(stderr,"ntr=%d is different from nh=%d \n",ntr,nh);
  if (nh>ntr) err("Change ntr in header\n");


  int cdpgather=tr.cdp;

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
  
  radontoepf(data,h,nh,t,nt,dt,model,q,nq,velint,eps1,factor,smute,nmofactor,rtmethod,depth,fmax);
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

  if (1) plotgather(model,nq,nt,0.004);

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

  free1float(Wd);
  free1float(Wm);
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






























