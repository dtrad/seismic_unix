/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUEOMIG:  $Date: March 1999  */

#include "header.h"
#include "eomig.h"
#include <signal.h>
#include <time.h>

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUCSP - CSPs from Equivalent Offset migration                       ",
  " 	   								",
  " sucsp < stdin > stdout [optional parameters]          		",
  " 									",
  " cdpmin=0        Fisrt CDP in meters                                 ",
  " cdpmax=100      Last CDP 						",
  " par=            file with stacking velocities and time as obtained  ",
  "                 from Velan                                          ",
  "                                                                     ",
  " output=2        =0 test - Output=input                              ",
  "                 =1 Migrated zero offset section                     ",
  "                 =2 CSP gather                                       ",
  " nh=100          Number of traces for csp                            ",
  " ncsp=10         Number of CSP to compute for every file  reading    ",
  " hemax=1000      Maximum and minimum equivalent offset               ",
  " hemin=1000      Minimum and minimum equivalent offset               ",
  " pervmin = 10    Minimum perturbation on velocity                    ",
  " dperv=0.1       Defines the rate of incresing velocity spacing      ",  
  "                 Velocity space to use in Radon                      ",
  " nq=20           Size of Radon EOM                                   ",
  " testhe=0        =1 Computes equiv offset at every time              ",
  "                 =2 Uses one side CSP                                ", 
  " beta=90         Maximum angle to migrate                            ",
  " buffer=5        Size of filter along time for velocity corrections  ", 
  " aper=(MAX(fabs(hemax),fabs(hemin))/cdpspace):Number of CMPs for CSP ",
  " verbose=0       =1  Extra Information                               ",
  " t0=0            First useful time                                   ",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset
 **************** end self doc ***********************************/

int main(int argc, char **argv)
{
  segy tr,trf; 
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
                        /* =2 output are CSP gathers,           */ 
  int *cdpt;         /* Vector with cdps number */
  int precise;       /* =1 Use velocity at to instead of t */
  int testhe;        /* =1 Compute he at every t           */
  float beta;        /* max angle (degrees) to include in migration */
  float cdpspace;    /* midpoint sampling interval; */
  float aper;

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
  if (!getparint("precise", &precise))  precise = 0;
  if (!getparint("testhe", &testhe))  testhe = 0;
  if (!getparfloat("beta", &beta))  beta = 45;
  if (!getparfloat("cdpspace",&cdpspace)) cdpspace = 1;
  if (!getparuint("ntrmax",&ntrmax)) ntrmax = 1000000;
  if (!getparuint("ncsp",&ncsp)) ncsp = 20 ;
  if (!getparushort("fold",&fold)) fold=1;
  if (!getparfloat("aper",&aper)) aper=(MAX(fabs(hmax/2),fabs(hmin/2))/cdpspace);
  if (!getparfloat("scalefold",&scalefold)) scalefold=0; //Default nofold scaling
  if (!getparfloat("scale",&scale)) scale=1;


  nx=(int) floor((cdpmax-cdpmin)/dxcdp + 1);
  // Because this implementation computes ncsp csps together
  // we need to modify nx such that mod(nx/ncsp) =0 
  int rest=(int) fmod(nx,ncsp);
  int nxold=nx;
  if (rest) nx+=(ncsp-rest);
  fprintf(stderr,"rest=%d,ncsp=%d,nx=%d,nxold=%d,scalefold=%f\n",rest,
	  ncsp,nx, nxold,scalefold); 
  fprintf(stderr,"Aperture in cdp number = %f\n",aper);


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
  trace=ealloc1float(nt);

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

      if (1){
	save_gather(csp[icsp],nh,h,nt,dt,"csp.su"); 
	system("cat csp.su >> csps.su");      
	system("suxwigb key=offset < csp.su title=csp perc=99 xbox=600 & \n");
      } 

      // Output trace
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
  free2float(velatcsp);
  free1float(velint);
  free1float(hec);
  free1float(tc);  
  free1float(t);
  free1float(he);
  free1float(h);
  free1float(x);
  free1float(trace);
  free3((void ***) cspfold);
  free3float(csp);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);


  return EXIT_SUCCESS;
}































