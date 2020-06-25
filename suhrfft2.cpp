/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SURADON00FORW  $Date: Septem 2000  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "hrfft2.h"
#include <signal.h>
#include <math.h>
#include <time.h>

int testfunction();

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUHRFFT2 - High resolution 2D FFT                                   ",
  " 	   								",
  " suhrfft2  < stdin > stdout [optional parameters]          		",
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
  " nk=80           Size of Radon EOM                                   ",
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

///////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  segy tr;
  cwp_String modelfile=""; /* output sufile for the model */ 
  cwp_String offsetfile=""; /* output sufile for the model */ 
  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;

  int nh;   // number of input offset traces
  int nh2;   // number of output offset traces
  int nt;   // number of time samples 
  float dt; // Time sampling
  float eps; // small tolerance number
  int testadj; // =1 test adjoint
  int smooth; // =1 smoothing filter
  ///////////////////////////////////////////////////////////////
  int  ih, ik; // General counters 
  register int it;
  float *d;      /* single trace */
  float **data;   /*  Input gather */
  float **data2;   /* Output gather */
  float  *m;    /* Final modeled trace   */
  float  *t;     // time axis for input and output
  float  *h=0;      // input offset
  float  *h2=0;      // output offset
  float t0;      // First useful time

  /* Velocity */
  float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  float *cdp;	        /* array[ncdp] of cdps */  
  int ncdp;	/* number of cdps specified */
  float *velint;/* array[nt] of vel for a particular trace */

  int verbose;		/* flag for echoing info		*/
  float fmax;
  int nk;  
  float kmin;
  float kmax;
  float *k;
  float dk;
  complex **model;            // Temporal array for the model
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
  int ntr;

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
  if (!getparfloat("kmin", &kmin))  kmin = -0.5;
  if (!getparfloat("kmax", &kmax))  kmax = 0.5;
  if (!getparint("nk", &nk))  nk = 100;
  if (!getparint("norm",&norm)) norm=0;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparfloat("factor",&factor)) factor=0.8;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=2;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
  if (!getparfloat("fmax",&fmax)) fmax=0.8/(2*dt);
  if (!getparstring("offsetfile",&offsetfile)) offsetfile=0;   

  ////////////////////////////////////////////////////////////////////////

  if ((modelfilep=fopen(modelfile,"w"))==NULL)
    err("cannot open file=%s\n",modelfile);
  
  // Allocate memory for data and model
  d=ealloc1float(nt);
  data=ealloc2float(nt,nh);
  model=alloc2complex(nk,2*nt);
  m=ealloc1float(nt);
  t=ealloc1float(nt);
  k=ealloc1float(nk);
  h=ealloc1float(nh);
  Wm=ealloc1float(nt*nk);
  Wd=ealloc1float(nt*nh);

  
  if (offsetfile){
    h2=ealloc1float(5*nh);
    nh2=read_ascii_file(offsetfile,h2);
    h2=erealloc1float(h2,nh2);
  }
  else{
    nh2=nh;
    h2=h;
  }

  data2=ealloc2float(nt,nh2);
  // Velocity axis for Radon space

  dk = (kmax-kmin)/(nk-1);
  for (ik=0;ik<nk;ik++) k[ik] = kmin+ik*dk;
  fprintf(stderr,"kmin=%f, kmax=%f \n", kmin,kmax);
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
  
  hrfft2_interface(data, model,h,nh,t,nt,k,nk,velint,eps1,smute,nmofactor,fmax,data2,h2,nh2);
  TRACE;
  // Output trace
  // if the number of input traces is equal to the number of output traces preserve the header
  if (nh==nh2){  
    for (ih=0;ih<nh2;ih++){ 
      gettr(&tr);
      memcpy((void *) tr.data,(const void *) data2[ih],nt*sizeof(float));
      puttr(&tr);
    }
  }
  else{
    gettr(&tr);
    for (ih=0;ih<nh2;ih++){ 
      tr.offset=h2[ih];
      tr.ntr=nh2;
      memcpy((void *) tr.data,(const void *) data2[ih],nt*sizeof(float));
      puttr(&tr);
    }
  }

  TRACE;
  /*
  if (0){
    save_gather(model,nk,nt,0.004,"model2");
    system("sufft  < model2 | suamp | suxwigb title=\"model2\" perc=99 & ");
  }
 
  if (1) plotgather(model,nk,nt,0.004);
  */
  TRACE;
  
  for (ik=0;ik<nk;ik++){
    //memcpy((void *) tr.data,(const void *) model[ik],nt*sizeof(float));
    for (it=0;it<nt;it++) tr.data[it]=abs(model[it][ik]);
    tr.dt=(int) (dt*1e6);       
    tr.ntr=nk;
    tr.ns=nt;
    tr.tracl=ik+1;
    tr.tracr=ik+1;
    tr.f2=k[ik];
    fputtr(modelfilep,&tr);    
  }
  
  TRACE;
  if (offsetfile) free1float(h2);
  free1float(cdp);
  free2float(ovv);
  free1float(velint);
  TRACE;
  free1float(Wd);
  free1float(Wm);
  free1float(k);
  TRACE;
  free1float(t);
  free1float(h);
  free1float(m);
  free2float(data);
  free2float(data2);
  TRACE;
  free2complex(model);
  free1float(d);
  efclose(modelfilep);
  TRACE;
  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}






























