/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUEOMIG:  $Date: March 1999  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "radonline_stack.h"
#include <signal.h>
#include <math.h>
#include <time.h>


/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADONLINE_STACK - Fast frequency domain Parabolic Radon Transform ",
  "                     with stack and mute                             ",
  " 	   								",
  " suradonline_stack < stdin > stdout [optional parameters]       		",
  "                                                                     ",
  " IMPORTANT:****************************************                  ",
  " Input must be sort by cdp and offset, for example                   ",
  "                   susort cdp offset < input | .....                 ",
  "                                                                     ",
  " It removes This program is intended to input a whole line of data.  ",
  "                                                                     ",
  " cdpmin=0        Fisrt CDP in meters                                 ",
  " cdpmax=100      Last CDP 						",
  " par=            file with stacking velocities and time as obtained  ",
  "                 from Velan (PARFILE)                                ",
  "                                                                     ", 
  " nhmax=200       Maximum number of traces for cmp                    ",
  " nq=100          Size of Radon EOM                                   ",
  " eps = 1e-1      hyperparameter for LS regularization                ",
  " mute=0          =1 mute multiples                                   ",
  " parnmute=1e-8   if q > parmute the radon trace is tapered           ",
  " verbose=0       =1  verbose output                                  ",
  " t0=0            First useful time                                   ",
  " quantil=1       Quantil to scale the model weight function          ",
  " smute=2         stretch greater than smutex100 % is muted		",
  " nmofactor=1     nmofactor * offset is used for nmo		        ",
  " pseudohyp=0     =1 Foster and Mosher                                ",
  "                                                                     ",
  " tests =0       No tests (fast operation)                            ", 
  "		   =1 Shows Radon space                                 ",
  "		   =2 Simple stack                                      ",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset
/**************** end self doc ***********************************/
segy tr,tr2;

int main(int argc, char **argv)
{

  cwp_String offsetfile=NULL;   /* output sufile for the cmp   */ 

  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;
  
  // Data space
  int nh;   // number of offset traces
  int nt;   // number of time samples 
  int nx;    // number of midpoints
  int nhmax; // max number of offsets.
  float dt; // Time sampling
  float dh; // offset interval
  float dx;  // midpoint interval
  float eps; // small tolerance number
  
  ///////////////////////////////////////////////////////////////
  int j, i, k, ih, ix; // General counters 
  register int it;
  float *d;      /* single trace */
  float **cmp;   /* Common midpoint gather */
  float *x;      /* axis for CDPs */
  float *h;      /* axis for input offset */
  float *t;     // time axis for input and output
  float t0;      // First useful time
  float *trace;  /* trace obtained by stack */
  unsigned int ntrmax;  // Number of traces to process from data file
  complex **oper; // stack+IRT+mute+RT operator
  /*** Arrays required for testing purposes only (option tests=2 ****/
  float *mutevector; /* Vector with 0 and 1 to filter out the Radon domain */
  float *stackvector; /* Vector with 0 and 1 to filter out the Radon domain */
  float **model; /* Radon domain gather required for testing purposes */
  complex ***LSoper;  /* LS Radon operator required for testing purposes */
  complex ***Loper;  /* Radon operator required for testing purposes */

  /* Velocity */
  int ncdp;	/* number of cdps specified */
  float *cdp;	/* array[ncdp] of cdps */
  int nvnmo;	/* number of vnmos specified */
  float *vnmo;	/* array[nvnmo] of vnmos */
  int ntnmo;	/* number of tnmos specified */
  float *tnmo;	/* array[ntnmo] of tnmos */
  float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  float *velint;/* array[nt] of vel for a particular trace */
  char *tmpdir;		/* directory path for tmp files		*/
  cwp_Bool istmpdir=cwp_false;/* true for user-given path		*/ 
  long ntr;     // Total number of input traces
  int itr;
  float  cdpmin;          //limit for the output
  float  cdpmax;          //limit for the output
  float  dxcdp;           // output cdp interval 
  int verbose;		/* flag for echoing info		*/

  // Radon axis  
  int nq;  
  float qmin;
  float *q;
  float dq;
  float dh_max;
  float dh_av;
  float qmaxt;
  float qmax;

  // Frequency axis
  int nfft,maxfreq,nf;
  float df;
  float fmax;  
  
  //  Options 
  int mute;
  float parmute;
  float factor;
  float smute;
  float nmofactor;
  int ngettr; // Number of bytes read by gettr (0 after)
  int pseudohyp; // Apply pseudohyperbolic RT
  float depth;   // depth for pseudohyperbolic
  int rtmethod;
  int tests; /* =0 No tests (fast operation) 
		=1 Shows Radon space 
		=2 Simple stack */
  // Filter
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
  tr2=tr; /* this trace has to be used later in the cdp loop */
  if (verbose) fprintf(stderr,"nt=%d,dt=%f\n",nt,dt);  
  
  // Get parameters 
  if (!getparfloat("cdpmin", &cdpmin))  cdpmin = 0;
  if (!getparfloat("cdpmax", &cdpmax))  cdpmax = 1000;
  if (!getparfloat("dxcdp", &dxcdp))  dxcdp = 1;
  if (!getparint("nhmax", &nhmax))  nhmax = 200;
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparfloat("eps", &eps))  eps = 1e-1;
  if (!getparfloat("qmin", &qmin))  qmin = -5e-9;
  if (!getparint("nq", &nq))  nq = 100;
  if (!getparuint("ntrmax",&ntrmax)) ntrmax = 1000000;
  if (!getparfloat("parmute",&parmute)) parmute=1e-8;
  if (!getparint("mute",&mute)) mute=0;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparfloat("factor",&factor)) factor=0.8;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=1;
  if (!getparstring("offsetfile",&offsetfile)) offsetfile=NULL;
  if (!getparint("pseudohyp",&pseudohyp)) pseudohyp=0;
  if (!getparint("rtmethod",&rtmethod)) rtmethod=2;
  if (!getparfloat("depth",&depth)) depth=2000;
  if (!getparfloat("fmax",&fmax)) fmax=0.8/(2*dt);
  if (!getparint("tests",&tests)) tests=0; 

  if (fmax==0) fmax=1./(2*dt);
  if (verbose) fprintf(stderr,"Maximum frequency ===> %f \n",fmax);
  if (verbose) fprintf(stderr,"Offset file ===> %s \n",offsetfile);
  ///////////////////////////////////////////////////////////////////////
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
  
  nx=(int) floor((cdpmax-cdpmin)/dxcdp + 1);

  if (verbose) fprintf(stderr,"number of ouput cdps = %d\n",nx);

  /* get velocity functions, linearly interpolated in time */
  // Velint will contain velocity law for one cmp
  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  getvelocities(dt,nt,ncdp,cdp,ovv);
  
  cmp=ealloc2float(nt,nhmax);
  model=ealloc2float(nt,nq);
  x=ealloc1float(nx);
  h=ealloc1float(nhmax);
  t=ealloc1float(nt);
  q=ealloc1float(nq);
  trace=ealloc1float(nt);


  /* Create axis for output cpds, offset, radon and time */
  dx=dxcdp;
  for(ix=0;ix<nx;ix++) x[ix]=cdpmin+dx*ix;
  
  for(it=0;it<nt;it++) t[it]=t0+it*dt;
  
  if (offsetfile) nh=read_ascii_file(offsetfile,h); 
  else err("This program needs an ascii file with the offset information");
  h=realloc1float(h,nh);

  if (verbose) fprintf(stderr,"Data dimensions: nh=%d, nt=%d\n",nh,nt);

  q[0]=qmin;
  interval(h,nh,&dh_max,&dh_av);
  fprintf(stderr,"dh_max=%f, dh_av=%f\n", dh_max, dh_av);
  radon_param(fmax,h,nh,dh_av,qmin,&qmaxt,&qmax,&dq,nq,rtmethod,factor);
  if (verbose) fprintf(stderr,"q max=%e,qmax used=%e, dq=%e\n", qmaxt,qmax,dq);
  for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 
  
  ///////////////////////////////////////////////////////////////////////
  fft_parameters(nt,dt,&nfft,&nf,&df);
  maxfreq=(int) (fmax/df);
  if (verbose) 
    fprintf(stderr,"dimensions of operator are maxfreq=%d times nh=%d\n",maxfreq,nh);
  oper=ealloc2complex(nh,maxfreq);
  if (tests==0) stack_rt_operator(oper,maxfreq,df,h,nh,q,nq,t,nt,eps,parmute);
  else if (tests==1){
    Loper=ealloc3complex(nq,nh,nf);
    LSoper=ealloc3complex(nh,nq,nf);
    mutevector=ealloc1float(nq);
    stackvector=ealloc1float(nq);
    stack_rt_LS_operators(LSoper,Loper,maxfreq,df,h,nh,q,nq,t,nt,eps,parmute,mutevector,
			  stackvector);
  }
  else if (tests==2) testop(oper,maxfreq,nh);
  /*************************************************************/
  /* Loop on the CMP gathers  */
  if (verbose) 
    fprintf(stderr,"Starting cdp loop ===> %d cmps to compute \n",nx);
  //rewind(stdin);
  for (ix=0;ix<nx;ix++){
    if (verbose) fprintf(stderr,"cdp=%f\n",x[ix]);   
    for (ih=0;ih<nh;ih++) memset(cmp[ih],(int) '\0',nt*FSIZE);
    /* put back last trace */
    tr=tr2;
    // search for the first cdp  
    while (tr.cdp<x[ix]) if (!(ngettr=gettr(&tr))) break;
    ntr = 0;
    do {
      ntr++;
      ih=ntr-1;
      if ((verbose)&&(1)) fprintf(stderr,"ntr=%d,tr.cdp=%d,tr.offset=%d\n",
			   ntr,tr.cdp,tr.offset);
      h[ih]=tr.offset;
      memcpy(cmp[ih],tr.data,nt*FSIZE);


    } while ( (ngettr=gettr(&tr)) && (tr.cdp == x[ix]));
    fprintf(stderr,"x=%f, ntr=%d\n",x[ix],ntr);
    /* Save last trace for next cdp */
    tr2=tr;
    if (ntr!=nh){
      warn("All cmp need to have the same geometry\n");
      warn("skipping this cdp ....\n");
      continue;
    }

    /* compute new square slowness  */
    interpovv(nt,ncdp,cdp,ovv,x[ix],velint);
    if (tests==0)
      radonline_stack0(cmp,oper,trace,h,nh,t,nt,dt,velint,smute,nmofactor,fmax,nfft,df); 
    else if (tests==1)
      radonline_LS(cmp,LSoper,Loper,model,trace,h,nh,t,nt,q,nq,dt,velint,smute,nmofactor,
		   fmax,nfft,df,mutevector,stackvector); 
    else if (tests==2)
      line_stack0(cmp,oper,trace,h,nh,t,nt,dt,velint,smute,nmofactor,fmax,nfft,df); 

    for (it=0;it<nt;it++)  tr.data[it]=trace[it];
    tr.cdp=(int) x[ix];
    tr.ntr=nx;
    
    fputtr(stdout,&tr);
  }

  /**************************************************/
  free1float(trace);
  free2complex(oper);
  
  if (namps) free1float(amps);
  
  if (npoly) free1float(ffilter);

  if (tests==1){
    free3complex(LSoper);
    free3complex(Loper);
    free1float(mutevector);
    free1float(stackvector);
  }
  free1float(cdp);
  free2float(ovv);
  free1float(q);
  free1float(velint);
  free1float(t);
  free1float(h);
  free1float(x);
  free2float(cmp);
  free2float(model);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);

  return EXIT_SUCCESS;
}































