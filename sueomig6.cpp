/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUEOMIG:  $Date: March 1999  */

#include "eomradontd.h"
#include "header.h"
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
  " output=2        =0 test - Output=input                              ",
  "                 =1 Migrated zero offset section                     ",
  "                 =2 CSP gather                                       ",
  "                 =3 Least square migrated traces			",
  "                 =4 Least square CSP gather                          ",
  "                 =5 Radon EOM domain                                 ",
  "                 =6 Velocity Correction                              ",
  " nh=100       Number of traces for csp                            ",
  " ncsp=10         Number of CSP to compute for every file  reading    ",
  " hemax=1000      Maximum and minimum equivalent offset               ",
  " hemin=1000      Minimum and minimum equivalent offset               ",
  " pervmin = 10    Minimum perturbation on velocity                    ",
  " dperv=0.1       Defines the rate of incresing velocity spacing      ",  
  "                 Velocity space to use in Radon                      ",
  " nq=20           Size of Radon EOM                                   ",
  " testhe=0        =1 Computes equiv offset at every time              ",
  "                 =2 Uses one side CSP                                ", 
  " testadj=0       =1 Test adjoint with random numbers                 ",
  " smooth=0        =1 Pass triangular filter to smooth  output         ",
  " itercg=5        Internal iterations for CG                          ",
  " iter_end=1      External iterations for CG                          ",
  " eps1=1e-3       Numerator hyperparameter for Wm                     ",
  " eps2=1e-3       Denominator  hyperparameter for Wm                  ",
  " beta=90         Maximum angle to migrate                            ",
  " buffer=5        Size of filter along time for velocity corrections  ", 
  " aper=(MAX(fabs(hemax),fabs(hemin))/cdpspace):Number of CMPs for CSP ",
  " mute=0          =1 mute nmute traces at the sides of Velocity space ",
  " nmute=4         number of Velocity traces to mute if mute=1         ",
  " verbose=0       =1  Extra Information                               ",
  " t0=0            First useful time                                   ",
  " keepwm=1        =1 The value for Wm from previous CSP is used       ",
  " sinc=0          =1 uses sinc interpolator                 		",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset
 **************** end self doc ***********************************/

/* Globals (so can trap signal) defining temporary disk files */
FILE *modelfilep;       /* fp for model file output  		*/
FILE *cspfilep;         /* fp for csp file output  		*/
///////////////////////////////////////////////////////////////////
// Global Variables 
int FLAG=1; // The first time (computes model_adj, then it just use previous csp
int KEEP_WM; // If keepwm==1 the value for Wm from previous CSP is used


int main(int argc, char **argv)
{
  segy tr,trf; 
  cwp_String modelfile=""; /* output sufile for the model */
  cwp_String cspfile="";   /* output sufile for the csp   */
  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;

  int nt;   // number of time samples 
  int nx;    // number of midpoints
  int nh; // number of equivalent offsets.
  float dt; // Time sampling
  float dh; // offset interval
  float dhe; // equivalent offset interval
  float dx;  // midpoint interval
  float eps; // small tolerance number
  int testadj; // =1 test adjoint
  int smooth; // =1 smoothing filter
  /// GLOBALS 
  extern int FLAG; // The first time (computes model_adj, then it just use previous csp
  extern int KEEP_WM; // If keepwm==1 the value for Wm from previous CSP is used
  
  ///////////////////////////////////////////////////////////////
  int i, ih, ix, iq; // General counters 
  register int it;
  float *d=0;      /* single trace */
  float ***csp=0;   /* Common Scattering gather */
  unsigned short ***cspfold=0;
  float  *trace=0;    /* Final migrated trace   */
  float  *m=0;    /* Final migrated trace   */
  float *x=0;      /* axis for CDPs */
  float *he=0;     /* axis for equivalent offset */
  float  *t=0;     // time axis for input and output
  float  *h=0;      // offset
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
  
  unsigned long int ntr;     // Total number of input traces
  const double  pi=acos(-1.); 
  float  cdpmin;          //limit for the output
  float  cdpmax;          //limit for the output
  float  dxcdp;           // output cdp interval 
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
  float *q;
  float dq;
  float *model;            // Temporal array for the model

  //  LS mig  
  float step;
  float eps1;
  float eps2;
  int itercg;
  int itervel;
  int iter_end;
  int buffer;
  float dperv;
  float pervmin;
  int centralq; // defines the location for the central velocity law.
  int norm;
  float *Wm;   // Factorization of inv model covariance matrix 
  float *Wd;   // Factorization of inv residual covariance matrix 
  int mute1;
  int mute2;
  float parmute1;
  float parmute2;
  float scalefold;
  float scale;
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
  if (!(ntr=(unsigned long int) tr.ntr)) err("***ntr must be set\n");
  
  fprintf(stderr,"nt=%d,dt=%f\n",nt,dt);  
  //fprintf(stderr,"ntr=%d\n",ntr);

  // Get parameters 
  
  if (!getparfloat("cdpmin", &cdpmin))  cdpmin = 0;
  if (!getparfloat("cdpmax", &cdpmax))  cdpmax = 1000;
  if (!getparfloat("dxcdp", &dxcdp))  dxcdp = 10;
  if (!getparfloat("hmin", &hmin))  hmin = 0;hmin=fabs(hmin);
  if (!getparfloat("hmax", &hmax))  hmax = 1000;hmax=fabs(hmax);
  if (!getparint("nh", &nh))  nh = 50;
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparint("output", &output)) output = 2;
  if (!getparfloat("eps", &eps))  eps = 1e-7;
  if (!getparfloat("eps1", &eps1))  eps1 = 1e-3;
  if (!getparfloat("eps2", &eps2))  eps2 = 1e-3;
  if (!getparfloat("step", &step))  step = 0.9;
  if (!getparint("itercg", &itercg))  itercg = 5;
  if (!getparint("itervel", &itervel))  itervel = 3;
  if (!getparint("iter_end", &iter_end))  iter_end = 1;
  if (!getparint("testadj", &testadj))  testadj = 0;
  if (!getparint("smooth", &smooth))  smooth = 0;
  if (!getparfloat("qmin", &qmin))  qmin = -200;
  if (!getparfloat("qmax", &qmax))  qmax = 200;
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
  if (!getparfloat("pervmin",&pervmin)) pervmin = 10;
  if (!getparfloat("dperv",&dperv)) dperv =0.1;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparint("keepwm",&KEEP_WM)) KEEP_WM=0;
  if (!getparint("centralq",&centralq)) centralq=0;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
  if (!getparstring("cspfile",&cspfile)) cspfile="cspfile.su";
  if (!getparfloat("scalefold",&scalefold)) scalefold=0; //Default nofold scaling
  if (!getparfloat("scale",&scale)) scale=1;
  if (!getparfloat("fmax",&fmax)) fmax=0.8/(2*dt);

  modelfilep=fopen(modelfile,"w");
  cspfilep=efopen(cspfile,"w");
  nx=(int) floor((cdpmax-cdpmin)/dxcdp + 1);
  // Because this implementation computes ncsp csps together
  // we need to modify nx such that mod(nx/ncsp) =0 
  int rest=(int) fmod(nx,ncsp);
  int nxold=nx;
  if (rest) nx+=(ncsp-rest);
  fprintf(stderr,"rest=%d,ncsp=%d,nx=%d,nxold=%d,scalefold=%f\n",rest,
	  ncsp,nx, nxold,scalefold); 
  fprintf(stderr,"Aperture in cdp number = %f\n",aper);

  if ( output > 2){
    if ( (nt*nh) < USHRT_MAX && (nt*nq) < USHRT_MAX ) 
      fprintf(stderr,"Using unsigned SHORT for index\n");
    else
      fprintf(stderr,"Using unsigned INT for index\n");
  }

  if (verbose){
    fprintf(stderr,"***Equivalent Offset Migration***\n");
    fprintf(stderr,"output=%d\n",output);
    fprintf(stderr,"number of ouput cdps = %d\n",nx);
    if (smooth) fprintf(stderr,"Smoothed Output \n");
    if (precise) fprintf(stderr,"precise computation of v[t0] \n");
    if (testhe) fprintf(stderr,"slow method: he at every t[it]\n");
    fprintf(stderr,"beta=%f\n",beta);
  }

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
  if ((cspfold=(unsigned short ***) alloc3(nt,nh,ncsp,sizeof(unsigned short)))==NULL); 

  x=ealloc1float(nx);
  he=ealloc1float(nh);
  h=ealloc1float(nh);
  t=ealloc1float(nt);
  tc=ealloc1float(nh);
  hec=ealloc1float(nh);
  q=ealloc1float(nq+1);
  model=ealloc1float(nt*nq);
  Wm=ealloc1float(nt*nq);
  Wd=ealloc1float(nt*nh);
  m=ealloc1float(nt);
  trace=ealloc1float(nt);
  // Velocity axis for Radon space

  if (nq>1){
    dq = (qmax-qmin)/(nq-1);
    for (i=0;i<nq;i++) q[i] = qmin+i*dq;
  }
  else  q[0]=0;


  /* Create axis for output cpds, equivalent offset and time */
  //if (testhe==2) hemin=0;
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

  ///////////////////////////////////////////////////////////////////////
  if (0) plotgather(ovv,ncdp,nt,0.004);
  //Test for velocity interpolation

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
      memcpy((void *) velatcsp[icsp],(const void *) velint,nt*sizeof(float));
      
      for (it=0;it<nt;it++) velatcsp[icsp][it]*=scale;
    }

    gettr(&tr);
    ntr=1;

    /* Loop for mapping one trace to the CSP gather */
    do{
      ntr++;

      // Compute the contribution of the read trace to every csp

      for (icsp=0;(icsp<ncsp)&&(x[ix+icsp]<=cdpmax);icsp++){
	//if (tr.cdp==1000) fprintf(stderr,"tr.offset=%d\n",tr.offset);
	/* compute new square slowness and anis function */
	interpovv(nt,ncdp,cdp,ovv,x[ix+icsp],velint);
	memcpy((void *) velatcsp[icsp],
	       (const void *) velint,nt*sizeof(float));
	/* Map this trace to the current csp 
	Three options:
        equiv_offset_test does all computations, (precise and slow)
	equiv_offset_1side precomputes times, and uses 1 side of offset
        equiv_offset precomputes times, and uses + and - offsets */
	if (fabs(tr.cdp-x[ix+icsp])<aper){// && tr.cdp!=1000){
	  if (testhe==1) 
	    equiv_offset_test(tr,csp[icsp],cspfold[icsp],t,x[ix+icsp],
			      velatcsp[icsp],dt,he,dhe,nt,nh,beta,cdpspace,scalefold);
	  else if (testhe==2)
	    equiv_offset_1side(tr.data,tr.offset,tr.cdp,csp[icsp],cspfold[icsp],t,
			       x[ix+icsp],velatcsp[icsp],dt,he,dhe,nt,nh,tc,hec,
			       precise,beta,cdpspace,scalefold);
	  else
	    equiv_offset(tr,csp[icsp],cspfold[icsp],t,x[ix+icsp],
			 velatcsp[icsp],dt,he,dhe,nt,nh,tc,hec,precise,
			 beta,cdpspace,scalefold);
	}
      }
      //fprintf(stderr,"tr.tracl=%d\n",tr.tracl);
    }while(gettr(&tr) && (ntr < ntrmax));

    fprintf(stderr,"Rewind(stdin) ....\n");
    erewind(stdin);
    
    for (icsp=0;(icsp<ncsp)&&(x[ix+icsp]<=cdpmax);icsp++){
      /* compute new square slowness and anis function */
      // interpovv(nt,ncdp,cdp,ovv,oa1,oa2,x[ix+icsp],velint,oa1t,oa2t); 
      if (fold){
	for (ih=0;ih<nh;ih++) 
	  for (it=0;it<nt;it++)
	    if (cspfold[icsp][ih][it]>0.5) 
	      csp[icsp][ih][it]/=(float) cspfold[icsp][ih][it];
      }
      if (((output==1)&&(x[ix+icsp]==csptest))||(1)){
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
      if (1){
	save_gather(csp[icsp],nh,h,nt,dt,"csp.su"); 
	system("cat csp.su >> csps.su");      
	system("suxwigb key=offset < csp.su title=csp perc=99 xbox=600 & \n");
      } 

      if (0) plotgather_pipe(csp[icsp],nh,nt,"before_rt");        
      if (output==3 || output==4 || output==5){
	if ( (nt*nh) < USHRT_MAX && (nt*nq) < USHRT_MAX ) 
	  eomigls_sparse(m,model,t,he,q,velatcsp[icsp],csp[icsp],nt,nh,nq,fmax,step,eps,eps1,eps2,itercg,iter_end,smooth,testadj,buffer,pervmin,dperv,norm,parmute1,mute1,t0,Wm,Wd,centralq);
	else
	  eomigls_sparse_large(m,model,t,he,q,velatcsp[icsp],csp[icsp],nt,nh,nq,fmax,step,eps,eps1,eps2,itercg,iter_end,smooth,testadj,buffer,pervmin,dperv,norm,parmute1,mute1,t0,Wm,Wd,centralq);
      }
      if (1){
	save_gather(model,nq,q,nt,dt,"rt.su"); 
	system("cat rt.su >> rts.su");      
	system("suxwigb key=f2 < rt.su title=rt perc=99 xbox=0  \n");
      }
      if (0) plotgather_pipe(model,nq,nt,"rt");        
      if (1){
	save_gather(csp[icsp],nh,h,nt,dt,"csppred.su"); 
	system("cat csppred.su >> cspspred.su");      
	system("suxwigb key=offset < csppred.su title=csppred perc=99 xbox=100  \n");
      } 



      if (output==6 || output ==7 ){ 
	if ( (nt*nh) < USHRT_MAX && (nt*nq) < USHRT_MAX ) 
	  eomgetvel(m,model,t,he,q,velatcsp[icsp],csp[icsp],nt,nh,nq,fmax,step,eps,eps1,eps2,itercg,iter_end,smooth,testadj,buffer,pervmin,dperv,norm,itervel,t0,Wm,Wd);    
	else  
	  eomgetvel_large(m,model,t,he,q,velatcsp[icsp],csp[icsp],nt,nh,nq,fmax,step,eps,eps1,eps2,itercg,iter_end,smooth,testadj,buffer,pervmin,dperv,norm,itervel,t0,Wm,Wd);        
      }

      // Output trace
      //if (removemult) 
      if (output==2 || output==4){
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
	  tr.offset=(int) (h[ih]);
	  tr.sx=(int) x[ix];
	  tr.gx=(int) x[ix];
	  puttr(&tr);    
	}
      }
      fprintf(stderr,"output=%d\n",output);
      if (output==5 || output ==7 || output ==4){
	for (iq=0;iq<nq;iq++){
	  for (it=0;it<nt;it++) tr.data[it]=model[iq*nt+it];
	  tr.cdp=(int) x[ix+icsp]; // front of the CSP
	  tr.dt=(int) (dt*1e6);       
	  tr.ntr=nq;
	  tr.ns=nt;
	  tr.tracl=iq+1;
	  tr.tracr=iq+1;
	  tr.f2=q[iq];
	  tr.sx=(int) x[ix];
	  tr.gx=(int) x[ix];
	  fprintf(stderr,"writing to file\n");
	  fputtr(modelfilep,&tr);    
	}
      }
    
      if (((output==1)||(output==3)||(output==6))&&(x[ix+icsp]<=cdpmax)){

	eomig1(csp[icsp],he,trace,t,velatcsp[icsp],t0,nt,nh);
	for (it=0;it<nt;it++)
	  trace[it]/=(velatcsp[icsp][0]*sqrt(pi*MAX(t[it],1e-2)));

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

  free1int(cdpt);
  free1float(cdp);
  free2float(ovv);
  free1float(Wd);
  free1float(Wm);
  free1float(model);
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
  free1float(m);

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































