/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUEOMIG:  $Date: March 1999  */

#include "header.h"
#include "eomig0freq.h"
#include "eomig.h"
#include <signal.h>
#include <time.h>

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUEOMIG - Radon Equivalent Offset migration                         ",
  " Prestack time migration with Radon transform                        ", 
  "	                                                          	",
  " 	   								",
  " sueomig0 < stdin > stdout [optional parameters]          		",
  " 									",
  " cdpmin=0        Fisrt CDP in meters                                 ",
  " cdpmax=100      Last CDP 						",
  " par=            file with stacking velocities and time as obtained  ",
  "                 from Velan                                          ",
  "                                                                     ",
  " output=1        =0 test - Output=input                              ",
  "                 =1 Migrated zero offset section                     ",
  "                 =2 CSP gather                                       ",
  "                 =3 Least square migrated traces			",
  "                 =4 Least square CSP gather                          ",
  "                 =5 Radon EOM domain                                 ",
  "                 =6 Velocity Correction                              ",
  " nh=100          Number of traces for csp                            ",
  " ncsp=10         Number of CSP to compute for every file  reading    ",
  " hmax=1000       Maximum and minimum offset                          ",
  " hmin=1000       Minimum and minimum offset                          ",
  "                 Velocity space to use in Radon                      ",
  " nq=20           Size of Radon EOM                                   ",
  " testhe=0        =1 Computes equiv offset at every time              ",
  "                 =2 Uses one side CSP                                ", 
  " smooth=0        =1 Pass triangular filter to smooth  output         ",
  " beta=90         Maximum angle to migrate                            ",
  " buffer=5        Size of filter along time for velocity corrections  ", 
  " aper=(MAX(fabs(hemax),fabs(hemin))/cdpspace):Number of CMPs for CSP ",
  " mute=0          =1 mute nmute traces at the sides of Velocity space ",
  " parnmute1=-1e-8  if q < parmute1 the radon trace is muted           ",
  " parnmute2=4e-8   if q > parmute2 the radon trace is muted            ",
  " verbose=0       =1  Extra Information                               ",
  " t0=0            First useful time                                   ",
  " keepwm=1        =1 The value for Wm from previous CSP is used       ",
  " sinc=0          =1 uses sinc interpolator                 		",
  " smute=2         stretch greater than smutex100 % is muted		",
  " nmofactor=1     nmofactor * offset is used for nmo	         	",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset */
/**************** end self doc ***********************************/
static void closefiles(void);
/* Globals (so can trap signal) defining temporary disk files */
FILE *modelfilep;       /* fp for model file output  		*/
FILE *cspfilep;         /* fp for csp file output  		*/
///////////////////////////////////////////////////////////////////
segy tr,trf;
int main(int argc, char **argv)
{

  cwp_String modelfile=""; /* output sufile for the model */ 
  cwp_String cspfile="";   /* output sufile for the csp   */ 
  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;

  ///////////////////////////////////////////////////////

  int nt;   // number of time samples 
  int nx;    // number of midpoints
  int nh; // number of equivalent offsets.
  float dt; // Time sampling
  float dhe; // equivalent offset interval
  float dh; // offset interval
  float dx;  // midpoint interval
  float eps; // small tolerance number
  int smooth; // =1 smoothing filter
  int KEEP_WM; // If keepwm==1 the value for Wm from previous CSP is used
  int rtmethod=2; /* =2 Parabolic RT */ 
  ///////////////////////////////////////////////////////////////
  int i, ih, iq, ix; /* General counters (offset, Radon, midpoint) */ 
  register int it;
  float *d;      /* single trace */
  float ***csp;   /* Common Scattering gather */
  unsigned short  ***cspfold;
  float  *trace;    /* Final migrated trace   */
  float  *trace2;    /* input trace to csp  */
  float *x;      /* axis for CDPs */
  float *he;     /* axis for equivalent offset */
  float *h;     /* axis for true offset */
  float  *t;     // time axis for input and output
  float t0;      // First useful time
  unsigned int ntrmax;  // Numebr of traces to process from data file
  unsigned int ncsp;  // Number of CSPs to process per reading loop
  unsigned int icsp;
  unsigned short fold;
  /* Velocity */
  int ncdp;	/* number of cdps specified */
  float *cdp;	/* array[ncdp] of cdps */
  float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  float *velint;/* array[nt] of vel for a particular trace */
  float **velatcsp; /* array[nt] of vel for a particular trace */ 

  unsigned int ntr;     // Total number of input traces
  const double  pi=acos(-1.); 
  float  cdpmin;          //limit for the output
  float  cdpmax;          //limit for the output
  float  dxcdp;           // output cdp interval 
  float  hemin;          //limit for the equivalent offset
  float  hemax;          //limit for the equivalent offset     
  float  hmin;          /* min offset */
  float  hmax;          /* max offset */    

  float  *hec;           /* equivalent offset for every trace */            
  float  *tc;            /* time when equivalent offset bin increases */  
  int verbose;		/* flag for echoing info		*/
  int output;           /* =0 output are original traces, used for test */
			/* =1 output are  migrated traces       */
                        /* =2 output are CSP gathers,           */ 
  float fmax;
  int *cdpt;         /* Vector with cdps number */
  int precise;       /* =1 Use velocity at to instead of t */
  int testhe;        /* =1 Compute he at every t           */
  float beta;        /* max angle (degrees) to include in migration */
  float cdpspace;    /* midpoint sampling interval; */
  float aper;
  
  // deviations from velocity law
  int nq;  
  float qmin;
  float qmax;
  float qmaxt;
  float *q;
  float dq;
  float **model;            // Temporal array for the model

  //  LS mig  
  int itervel;
  int buffer;
  int norm;
  float *Wm;   // Factorization of inv model covariance matrix 
  float *Wd;   // Factorization of inv residual covariance matrix 
  int mute1;
  int mute2;
  float parmute1;
  float parmute2;
  int sinc; //=1 uses sinc interpolation
  float quantil;
  float factor;
  float depth;   // depth for pseudohyperbolic 
  float smute;
  float nmofactor;
  float scalefold;
  float t0mute;
  float scale;
  /* Options */
  int tests;
  int compute; /* =1 compute operator, =0 read saved file */
  int csptest; /* CSP at a particular midpoint location can be saved in the cspfile */ 
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
  // Register starting time
  start=time(0);
   // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  erewind(stdin); /* the first trace is read again later */

  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");

  //if (!tr.offset) err("offset header field must be set");
  
  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  if (!(ntr=(unsigned int) tr.ntr)) err("***ntr must be set\n");
  fprintf(stderr,"nt=%d,dt=%f\n",nt,dt);  
  fprintf(stderr,"ntr=%d\n",(int) ntr);

  // Get parameters 
  
  if (!getparfloat("cdpmin", &cdpmin))  cdpmin = 0;
  if (!getparfloat("cdpmax", &cdpmax))  cdpmax = 1000;
  if (!getparfloat("dxcdp", &dxcdp))  dxcdp = 10;
  if (!getparfloat("hmin", &hmin))  hmin = 0; hmin=fabs(hmin);
  if (!getparfloat("hmax", &hmax))  hmax = 1000; hmax=fabs(hmax);
  if (!getparint("nh", &nh))  nh = 50;
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparint("output", &output)) output = 1;
  if (!getparfloat("eps", &eps))  eps = 1e-7;
  if (!getparint("itervel", &itervel))  itervel = 3;
  if (!getparint("smooth", &smooth))  smooth = 0;
  if (!getparfloat("qmin", &qmin))  qmin =0;
  if (!getparfloat("qmax", &qmax))  qmax =1e-8;
  if (!getparint("nq", &nq))  nq = 20;
  if (!getparint("precise", &precise))  precise = 0;
  if (!getparint("testhe", &testhe))  testhe = 0;
  if (!getparint("norm",&norm)) norm=0;
  if (!getparfloat("beta", &beta))  beta = 45;
  if (!getparint("buffer", &buffer))  buffer = 5;
  if (!getparfloat("cdpspace",&cdpspace)) cdpspace = 1;
  if (!getparuint("ntrmax",&ntrmax)) ntrmax = 1000000;
  if (!getparuint("ncsp",&ncsp)) ncsp = 20 ;
  if (!getparfloat("parmute1",&parmute1)) parmute1=-1e-9;
  if (!getparfloat("parmute2",&parmute2)) parmute2=4e-9;
  if (!getparint("mute1",&mute1)) mute1=0;
  if (!getparint("mute2",&mute2)) mute2=0;
  if (!getparushort("fold",&fold)) fold=1;
  if (!getparfloat("aper",&aper)) aper=(MAX(fabs(hmax/2),fabs(hmin/2))/cdpspace);
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparfloat("t0mute",&t0mute)) t0mute=0;
  if (!getparint("keepwm",&KEEP_WM)) KEEP_WM=0;
  if (!getparint("sinc",&sinc)) sinc=0;
  if (!getparfloat("factor",&factor)) factor=0.8;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=1;
  if (!getparfloat("scalefold",&scalefold)) scalefold=1;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
  if (!getparstring("cspfile",&cspfile)) cspfile="cspfile.su";
  if (!getparfloat("depth",&depth)) depth=2000;
  if (!getparfloat("fmax",&fmax)) fmax=0.8/(2*dt);
  if (!getparfloat("quantil",&quantil)) quantil=1;
  /* Some options */

  if (!getparfloat("scale",&scale)) scale=1;
  if (!getparint("tests",&tests)) tests=0;
  if (!getparint("compute",&compute)) compute=1;
  if (!getparint("csptest",&csptest)) csptest=0;

  if (fmax==0) fmax=1./(2*dt);

  modelfilep=efopen(modelfile,"w");
  cspfilep=efopen(cspfile,"w");

  nx=(int) floor((cdpmax-cdpmin)/dxcdp + 1);
  // Because this implementation computes ncsp csps together
  // we need to modify nx such that mod(nx/ncsp) =0 
  int rest=(int) fmod(nx,ncsp);
  int nxold=nx;
  if (rest) nx+=(ncsp-rest);
  if (1) fprintf(stderr,"rest=%d,ncsp=%d,nx=%d,nxold=%d\n",rest,ncsp,nx,nxold); 
  fprintf(stderr,"Aperture in cdp number = %f\n",aper);

  if (verbose) fprintf(stderr,"***Equivalent Offset Migration***\n");
  if (verbose) fprintf(stderr,"output=%d\n",output);
  if (verbose) fprintf(stderr,"number of ouput cdps = %d\n",nx);
  if (verbose && smooth) fprintf(stderr,"Smoothed Output \n");
  if (verbose && precise) fprintf(stderr,"precise computation of v[t0] \n");
  if (verbose && testhe) fprintf(stderr,"slow method: he at every t[it]\n");
  if (verbose)  fprintf(stderr,"beta=%f\n",beta);

  beta=beta*pi/180.;
  /* get velocity functions, linearly interpolated in time */
  // Velint will contain velocity law for one csp
  // velatcsp will contain velocity law for ncsps 
  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  velatcsp=ealloc2float(nt,ncsp);
  getvelocities(dt,nt,ncdp,cdp,ovv);
  for (ix=0;ix<ncdp;ix++) 
    fprintf(stderr,"cdp[%d]=%f\n",ix,cdp[ix]); 
  cdpt=ealloc1int(ntr+1);

  // Allocate memory for data and model
  d=ealloc1float(nt);
  csp=ealloc3float(nt,nh,ncsp);
  trace=ealloc1float(nt);
  trace2=ealloc1float(nt);
  x=ealloc1float(nx);
  he=ealloc1float(nh); /* equivalent offset */
  h=ealloc1float(nh);  /* true offset */
  t=ealloc1float(nt);
  tc=ealloc1float(nh+2);
  hec=ealloc1float(nh+2);

  q=ealloc1float(nq+1);
  model=alloc2float(nt,nq);
  Wm=alloc1float(nt*nq);
  Wd=alloc1float(nt*nh);
  if ((cspfold=(unsigned short ***) 
       alloc3(nt,nh,ncsp,sizeof(unsigned short)))==NULL) 
    fprintf(stderr,"***Sorry, space for csp could not be allocated\n");
  


  /* Create axis for output cpds, equivalent offset and time */
  //if (testhe==2) hmin=0;
  dx=dxcdp;
  /* he is the equivalent offset  */
  /* h is the offset = 2 * he */
  dh=(hmax-hmin)/(nh-1);dhe=dh/2;      
  for(ix=0;ix<nx;ix++) x[ix]=cdpmin+dx*ix;
  for(it=0;it<nt;it++) t[it]=0+it*dt;
  for(ih=0;ih<nh;ih++) h[ih]=hmin+ih*dh;
  for(ih=0;ih<nh;ih++) he[ih]=h[ih]/2.;
  if (verbose) fprintf(stderr,"hmin=%f,hmax=%f,dh=%f\n",h[0],h[nh-1],dh);
  if (verbose) fprintf(stderr,"hemin=%f,hemax=%f,dhe=%f\n",he[0],he[nh-1],dhe);

  // Velocity axis for Radon space
  q[0]=qmin;
  radon_param(fmax,h,nh,dh,qmin,&qmaxt,&qmax,&dq,nq,rtmethod,factor);
  if (verbose) fprintf(stderr,"q max=%e,qmax used=%e, dq=%e\n", qmaxt,qmax,dq);
  for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 

  ///////////////////////////////////////////////////////////////////////
  // Plot velocity model
  if (0) plotgather(ovv,ncdp,nt,0.004);
  // PLot interpolated velocity model
  float **velgather=NULL;
  if(1){
    velgather=ealloc2float(nt,nx);
    for (ix=0;ix<nx;ix++){
      interpovv(nt,ncdp,cdp,ovv,x[ix],velint);
      memcpy((void *) velgather[ix],(const void *) velint,nt*sizeof(float));
    }
    plotgather(velgather,nx,nt,0.004);
    free2float(velgather);
  }
  ///////////////////////////////////////////////////////////////////////

  /** Create SIRMUR operator  */

  complex **oper=0; // stack+IRT+mute+RT operator
  complex ***LSoper=0;  /* LS Radon operator required for testing purposes */
  complex ***Loper=0;  /* Radon operator required for testing purposes */
  float *mutevector=0; /* Vector with 0 and 1 to filter out the Radon domain */
  float *stackvector=0; /* Vector with 0 and 1 to filter out the Radon domain */
  // Frequency axis
  int nfft,maxfreq,nf;
  float df;


  fft_parameters(nt,dt,&nfft,&nf,&df);
  maxfreq=(int) (fmax/df);
  oper=ealloc2complex(nh,maxfreq);
  if (tests==0){
    if (compute){
      stack_rteom_operator(oper,maxfreq,df,h,nh,q,nq,t,nt,eps,parmute1,parmute2);
      save2dfile(oper,nh,maxfreq,"oper");
    }
    else read2dfile(oper,nh,maxfreq,"oper");
  }
  else if (tests==1){
    Loper=ealloc3complex(nq,nh,nf);
    LSoper=ealloc3complex(nh,nq,nf);
    mutevector=ealloc1float(nq);
    stackvector=ealloc1float(nq);
    if (compute) 
      stack_rteom_LS_operators(LSoper,Loper,maxfreq,df,h,nh,q,nq,t,nt,eps,parmute1,parmute2,
			    mutevector,stackvector);
    /* if repeated test are required with the same operator, save them */
    if (compute){
      save3dfile(LSoper,nh,nq,nf,"LSoper");
      save3dfile(Loper,nq,nh,nf,"Loper");
      //save1dfile(mutevector,nq,"mutevector");
      save1dfile(stackvector,nh,"stackvector");
    }
    /* and read them */
    if (!compute){
      mute_vector(mutevector,q,nq,parmute1,parmute2);
      read3dfile(LSoper,nh,nq,nf,"LSoper");
      read3dfile(Loper,nq,nh,nf,"Loper");
      //read1dfile(mutevector,nq,"mutevector");
      read1dfile(stackvector,nh,"stackvector");
    }

  }


  
  /************************************************************************/
  fprintf(stderr,"Starting cdp loop ****\n");  
  /* Loop to compute ncsp CSP gathers from x[ix] to x[ix+ncsp] */
  for (ix=0;ix<nx;ix+=ncsp){
    fprintf(stderr,"cdp=%f to cdp=%f\n",x[ix],x[ix+ncsp-1]);   
    
    /* Initialize to zero the ouput trace and the csp gather corresp to ix */  
    for (icsp=0;icsp<ncsp;icsp++){
      memset( (void *) cspfold[icsp][0], (int) '\0', nh*nt *sizeof(unsigned short));
      memset( (void *) csp[icsp][0], (int) '\0', nh*nt *FSIZE);
    
      /* compute new square slowness and anis function */
      if (0)  fprintf(stderr,"Compute velocity trend for csp...\n");
      interpovv(nt,ncdp,cdp,ovv,x[ix+icsp],velint);
      
      memcpy((void *) velatcsp[icsp],
	     (const void *) velint,nt*sizeof(float));
      
      for (it=0;it<nt;it++) velatcsp[icsp][it]*=scale;
    }
    gettr(&tr);
    ntr=1;
    
    /* Loop for mapping one trace to the CSP gather */
    fprintf(stderr,"Reading data ... \n");
    do{
      ntr++;

      // Compute the contribution of the read trace to every csp
      for (icsp=0;(icsp<ncsp)&&(x[ix+icsp]<=cdpmax);icsp++){
	/* Map this trace to the current csp 
	   Three options:
	   equiv_offset_test does all computations, (precise and slow)
	   equiv_offset_1side precomputes times, and uses 1 side of offset
	   equiv_offset precomputes times, and uses + and - offsets */
	if (fabs(tr.cdp-x[ix+icsp])<aper){

	  
	  if (testhe==1) 
	    equiv_offset_test(tr,csp[icsp],cspfold[icsp],t,x[ix+icsp],
			      velatcsp[icsp],dt,he,dhe,nt,nh,beta,cdpspace,scalefold);
	  /* sending segy structures gives problems, send the members instead */
	  else if (testhe==2){
	    equiv_offset_1side(tr.data,tr.offset,tr.cdp,csp[icsp],cspfold[icsp],t,
			       x[ix+icsp],velatcsp[icsp],dt,he,dhe,nt,nh,tc,hec,
			       precise,beta,cdpspace,scalefold);
	  }
	  else
	    equiv_offset(tr,csp[icsp],cspfold[icsp],t,x[ix+icsp],
			 velatcsp[icsp],dt,he,dhe,nt,nh,tc,hec,precise,
			 beta,cdpspace,scalefold);
	}
	//else fprintf(stderr,"Out ...\n");
      }
      //fprintf(stderr,"icsp=%d\n",icsp);;
    }while(gettr(&tr) && (ntr < ntrmax));

    
    fprintf(stderr,"Rewind(stdin) ....\n");
    erewind(stdin);
    fprintf(stderr,"calculating RT....\n");


    if(0){
      velgather=ealloc2float(nt,ncsp);
      for (icsp=0;icsp<ncsp;icsp++){
	interpovv(nt,ncdp,cdp,ovv,x[ix+icsp],velint);
	memcpy((void *) velgather[icsp],(const void *) velint,nt*sizeof(float));
      }
      plotgather(velgather,ncsp,nt,0.004);
      free2float(velgather);
    }

    for (icsp=0;(icsp<ncsp)&&(x[ix+icsp]<=cdpmax);icsp++){
      /* compute new square slowness and anis function */
      // interpovv(nt,ncdp,cdp,ovv,oa1,oa2,x[ix+icsp],velint,oa1t,oa2t); 
      if (fold)
	for (ih=0;ih<nh;ih++) 
	  for (it=0;it<nt;it++)
	    if (cspfold[icsp][ih][it]>0.5) 
	      csp[icsp][ih][it]/=(float) cspfold[icsp][ih][it];
      
      if ((output==1)&&(x[ix+icsp]==csptest)){
	for (ih=0;ih<nh;ih++){
	  memcpy((void *) tr.data,
		 (const void *) csp[icsp][ih],nt*sizeof(float));
	  // for (it=0;it<nt;it++) tr.data[it]=csp[icsp][ih][it];
	  tr.cdp=(int) x[ix+icsp]; // front of the CSP
	  tr.dt=(int) (dt*1e6);       
	  tr.ntr=nh;
	  tr.ns=nt;
	  tr.tracl=ih+1;
	  tr.tracr=ih+1;
	  tr.offset=(int) (2*he[ih]);
	  tr.sx=(int) x[ix];
	  tr.gx=(int) x[ix];
	  fputtr(cspfilep,&tr);    
	}
      }
      
      if (output==3 || output==4 || output==5 || (output==2 && tests==1)){
	if (tests==0)
	  radonline_stack0(csp[icsp],oper,trace,h,nh,t,nt,dt,velatcsp[icsp],smute,
			   nmofactor,fmax,nfft,df); 
	else if (tests==1)
	  radonline_LS(csp[icsp],LSoper,Loper,model,trace,h,nh,t,nt,q,nq,dt,velatcsp[icsp],
		       smute, nmofactor,fmax,nfft, df, mutevector, stackvector); 
	else if ((tests==2)&&(0))  /* line_stack0 is only for testing */
	  /* For normal migration used eomig1 at the end of the file */
	  line_stack0(csp[icsp],oper,trace,h,nh,t,nt,dt,velatcsp[icsp],smute,nmofactor,
		      fmax,nfft,df); 
      }
      
      //if (0) radonsolver0_mute(csp[icsp],he,nh,t,nt,dt,model,q,nq,velatcsp[icsp],itercg,iter_end,step,eps2,eps1, quantil,norm,factor,smute,nmofactor,rtmethod,depth,fmax,parmute1,parmute2,mute1,mute2,t0mute); 
      
      // Output trace
      
      if (output==2 || output==4){ /* CSP gather */
	for (ih=0;ih<nh;ih++){
	  memcpy((void *) tr.data,
		 (const void *) csp[icsp][ih],nt*sizeof(float));
	  tr.cdp=(int) x[ix+icsp]; // front of the CSP
	  tr.dt=(int) (dt*1e6);       
	  tr.ntr=nh;
	  tr.ns=nt;
	  tr.tracl=ih+1;
	  tr.tracr=ih+1;
	  tr.offset=(int) (2*he[ih]);
	  tr.sx=(int) x[ix];
	  tr.gx=(int) x[ix];
	  puttr(&tr);    
	}
      }
      
      if (output==5 || output ==7 || output ==4){ /** Radon space **/
	for (iq=0;iq<nq;iq++){
	  for (it=0;it<nt;it++) tr.data[it]=model[iq][it];
	  tr.cdp=(int) x[ix+icsp]; // front of the CSP
	  tr.dt=(int) (dt*1e6);       
	  tr.ntr=nq;
	  tr.ns=nt;
	  tr.tracl=iq+1;
	  tr.tracr=iq+1;
	  tr.offset=(int) (q[iq]);
	  tr.f2=q[iq];
	  tr.sx=(int) x[ix];
	  tr.gx=(int) x[ix];
	  fputtr(modelfilep,&tr);    
	}
      }
      
      if (((output==1)||(output==3)||(output==6))&&(x[ix+icsp]<=cdpmax)){ /*Migrated section*/
	if (tests==2) eomig1(csp[icsp],he,trace,t,velatcsp[icsp],t0,nt,nh);
	memcpy(tr.data, trace, nt*sizeof(float));
	tr.cdp=(int) x[ix+icsp];
	tr.dt=(int) (dt*1e6);
	tr.ntr=nxold;
	tr.ns=nt;
	tr.tracl=ix+1+icsp;
	tr.tracr=ix+1+icsp;
	tr.offset=0;
	tr.sx=(int) x[ix+icsp];
	tr.gx=(int) x[ix+icsp];
	filt(tr.data,nt,dt,fmax,0,50,trf.data);
	puttr(&tr);
	fprintf(stderr,"tr.cdp=%d\n",tr.cdp);    
	
      } 
      
    }
    
  }
  
  /* SIRMUR operator */
  free2complex(oper);
  
  if (tests==1){
    free3complex(LSoper);
    free3complex(Loper);
    free1float(mutevector);
    free1float(stackvector);
  }
  
  /*******************/
  free1int(cdpt);
  free1float(cdp);
  free2float(ovv);
  free1float(Wd);
  free1float(Wm);
  free2float(model);
  free1float(q);
  free2float(velatcsp);
  free1float(velint);
  free1float(hec);
  free1float(tc);  
  free1float(t);
  free1float(he);
  free1float(h);
  free1float(x);
  free1float(trace);
  free1float(trace2);
  free3((void ***) cspfold);
  
  free3float(csp);
  free1float(d);
  efclose(modelfilep);
  efclose(cspfilep);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);


  return EXIT_SUCCESS;
}


/* for graceful interrupt termination */
static void closefiles(void)
{
  //	efclose(stdout);
	exit(EXIT_FAILURE);
}


































