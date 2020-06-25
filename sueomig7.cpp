/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUEOMIG:  $Date: March 1999  */

#include "su.h"
#include "segy.h"
#include "clibrarytd.h"
#include "header.h"
#include <signal.h>
#include <math.h>
#include <time.h>
#include "eomig.h"

/* Maximum value an `unsigned short int' can hold.  (Minimum is 0).  */
#ifndef USHRT_MAX
#define USHRT_MAX 65535
#endif

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
  " nhcsp=100       Number of traces for csp                            ",
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
/**************** end self doc ***********************************/
static void closefiles(void);

/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
FILE *modelfilep;       /* fp for model file output  		*/
FILE *cspfilep;         /* fp for csp file output  		*/
///////////////////////////////////////////////////////////////////
// Global Variables 

int nh;   // number of offset traces
int nt;   // number of time samples 
int nx;    // number of midpoints
int nhcsp; // number of equivalent offsets.
float dt; // Time sampling
float dh; // offset interval
float dhe; // equivalent offset interval
float dx;  // midpoint interval
float eps; // small tolerance number
int testadj; // =1 test adjoint
int smooth; // =1 smoothing filter
int FLAG=1; // The first time (computes mm_adj, then it just use previous csp
int KEEP_WM; // If keepwm==1 the value for Wm from previous CSP is used

int main(int argc, char **argv)
{
  segy tr,trf; 
  cwp_String modelfile=""; /* output sufile for the model */
  cwp_String cspfile="";   /* output sufile for the csp   */
  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;

  //////////// Global///////////////////////////////////////////
  extern int  nhcsp;      //Maximum number of traces for the CSP
  extern int  nh;
  extern int  nt;
  extern int  nx;         // number of output cdps
  extern float dt,dh,dx,dhe,eps;
  extern int testadj;
  extern int smooth; 	
  extern int KEEP_WM;
  ///////////////////////////////////////////////////////////////
  int j, i, k, ih, ix; // General counters 
  register int it;
  float *d;      /* single trace */
  float ***csp;   /* Common Scattering gather */
  unsigned short ***cspfold;
  float  *m;    /* Final migrated trace   */
  float *x;      /* axis for CDPs */
  float *he;     /* axis for equivalent offset */
  float  *t;     // time axis for input and output
  float  h;      // halfoffset
  float t0;      // First useful time
  unsigned int ntrmax;  // Numebr of traces to process from data file
  unsigned int ncsp;  // Number of CSPs to process per reading loop
  unsigned int icsp;
  unsigned short fold;
  /* Velocity */
  int ncdp;	/* number of cdps specified */
  float *cdp;	/* array[ncdp] of cdps */
  int icdp;	/* index into cdp array */
  int jcdp;	/* index into cdp array */
  int nvnmo;	/* number of vnmos specified */
  float *vnmo;	/* array[nvnmo] of vnmos */
  int ntnmo;	/* number of tnmos specified */
  float *tnmo;	/* array[ntnmo] of tnmos */
  float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  float *velint;/* array[nt] of vel for a particular trace */
  float **velatcsp; /* array[nt] of vel for a particular trace */ 
  int nanis1;	/* number of anis1's specified */
  int nanis2;	/* number of anis2's specified */
  float *anis1;	/* array[nanis1] of anis1's */
  float *anis2;	/* array[nanis2] of anis2's */
  float **oa1;	/* array[ncdp][nt] of anis1 functions */
  float **oa2;	/* array[ncdp][nt] of anis2 functions */
  float *oa1t;	/* array[nt] of anis1 for a particular trace */
  float *oa2t;	/* array[nt] of anis2 for a particular trace */
  float tn;     /* temporary time */
  float acdp;	/* temporary used to sort cdp array */
  float *aovv;	/* temporary used to sort ovv array */
  float *aoa1;	/* temporary used to sort oa1 array */
  float *aoa2;	/* temporary used to sort oa2 array */
  char *tmpdir;		/* directory path for tmp files		*/
  cwp_Bool istmpdir=cwp_false;/* true for user-given path		*/ 

  unsigned long int ntr;     // Total number of input traces
  unsigned long int itr;
  const double  pi=acos(-1.); 
  float xm;    // midpoint distance from trace to scp = cdp_tr - x_scp

  float  cdpmin;          //limit for the output
  float  cdpmax;          //limit for the output
  float  dxcdp;           // output cdp interval 
  float  hemin;          //limit for the equivalent offset
  float  hemax;          //limit for the equivalent offset     
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
  float *mm;            // Temporal array for the model

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
  unsigned short mute;
  float parmute;
  float scalefold;
  int sinc; //=1 uses sinc interpolation

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
  if (!(ntr=(unsigned long int) tr.ntr)) err("***ntr must be set\n");
  fmax=0.8/(2*dt);
  
  fprintf(stderr,"nt=%d,dt=%f\n",nt,dt);  
  fprintf(stderr,"ntr=%d\n",ntr);

  // Get parameters 
  
  if (!getparfloat("cdpmin", &cdpmin))  cdpmin = 0;
  if (!getparfloat("cdpmax", &cdpmax))  cdpmax = 1000;
  if (!getparfloat("dxcdp", &dxcdp))  dxcdp = 10;
  if (!getparfloat("hemin", &hemin))  hemin = 0;
  if (!getparfloat("hemax", &hemax))  hemax = 1000;
  if (!getparint("nhcsp", &nhcsp))  nhcsp = 50;
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparint("output", &output)) output = 1;
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
  if (!getparfloat("parmute",&parmute)) parmute=0;
  if (!getparushort("mute",&mute)) mute=0;
  if (!getparushort("fold",&fold)) fold=1;
  if (!getparfloat("aper",&aper)) aper=(MAX(fabs(hemax),fabs(hemin))/cdpspace);
  if (!getparfloat("pervmin",&pervmin)) pervmin = 10;
  if (!getparfloat("dperv",&dperv)) dperv =0.1;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparint("keepwm",&KEEP_WM)) KEEP_WM=0;
  if (!getparint("sinc",&sinc)) sinc=0;
  if (!getparint("centralq",&centralq)) centralq=0;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
  if (!getparstring("cspfile",&cspfile)) cspfile="cspfile.su";
  if (!getparfloat("scalefold",&scalefold)) scalefold=1;
  if ((modelfilep=fopen(modelfile,"w"))==NULL)
    err("cannot open file=%s\n",modelfile);
  cspfilep=efopen(cspfile,"w");
  nx=(int) floor((cdpmax-cdpmin)/dxcdp + 1);
  // Because this implementation computes ncsp csps together
  // we need to modify nx such that mod(nx/ncsp) =0 
  int rest=(int) fmod(nx,ncsp);
  int nxold=nx;
  if (rest) nx+=(ncsp-rest);
  if (1) fprintf(stderr,"rest=%d,ncsp=%d,nx=%d,nxold=%d,scalefold=%f\n",rest,ncsp,nx,
		 nxold,scalefold); 
  fprintf(stderr,"Aperture in cdp number = %f\n",aper);

  if ( output > 2){
    if ( (nt*nhcsp) < USHRT_MAX && (nt*nq) < USHRT_MAX ) 
      fprintf(stderr,"Using unsigned SHORT for index\n");
    else
      fprintf(stderr,"Using unsigned INT for index\n");
  }

  if (verbose) fprintf(stderr,"***Equivalent Offset Migration***\n");
  if (verbose) fprintf(stderr,"output=%d\n",output);
  if (verbose) fprintf(stderr,"number of ouput cdps = %d\n",nx);
  if (verbose && smooth) fprintf(stderr,"Smoothed Output \n");
  if (verbose && precise) fprintf(stderr,"precise computation of v[t0] \n");
  if (verbose && testhe) fprintf(stderr,"slow method: he at every t[it]\n");
  if (verbose)  fprintf(stderr,"beta=%f\n",beta);

  beta=beta*pi/180.;
  /* get velocity functions, linearly interpolated in time */
  ncdp = countparval("cdp");
  if (ncdp>0) {
    if (countparname("vnmo")!=ncdp)
      err("a vnmo array must be specified for each cdp");
    if (countparname("tnmo")!=ncdp)
      err("a tnmo array must be specified for each cdp");
    if (countparname("anis1")!=ncdp &&
	countparname("anis1")!=0)
      err("an anis1 array must be specified for each cdp, "
	  "or omitted at all");
    if (countparname("anis2")!=ncdp &&
	countparname("anis2")!=0)
      err("an anis2 array must be specified for each cdp, "
	  "or omitted at all");
  } else {
    ncdp = 1;
    if (countparname("vnmo")>1)
      err("only one (or no) vnmo array must be specified");
    if (countparname("tnmo")>1)
      err("only one (or no) tnmo array must be specified");
    if (countparname("anis1")>1)
      err("only one (or no) anis1 array must be specified");
    if (countparname("anis2")>1)
      err("only one (or no) anis2 array must be specified");    
  }

  cdp = ealloc1float(ncdp);
  if (!getparfloat("cdp",cdp)) cdp[0] = tr.cdp;
  ovv = ealloc2float(nt,ncdp);
  oa1 = ealloc2float(nt,ncdp);
  oa2 = ealloc2float(nt,ncdp);
  for (icdp=0; icdp<ncdp; ++icdp) {
    nvnmo = countnparval(icdp+1,"vnmo");
    ntnmo = countnparval(icdp+1,"tnmo");
    nanis1 = countnparval(icdp+1,"anis1");
    nanis2 = countnparval(icdp+1,"anis2");
    if (nvnmo!=ntnmo && !(ncdp==1 && nvnmo==1 && ntnmo==0))
      err("number of vnmo and tnmo values must be equal");
    if (nanis1!=nvnmo && nanis1 != 0)
      err("number of vnmo and anis1 values must be equal");
    if (nanis2!=nvnmo && nanis2 != 0)
      err("number of vnmo and anis2 values must be equal");
    if (nvnmo==0) nvnmo = 1;
    if (ntnmo==0) ntnmo = nvnmo;
    if (nanis1==0) nanis1 = nvnmo;
    if (nanis2==0) nanis2 = nvnmo;
    /* equal numbers of parameters vnmo, tnmo, anis1, anis2 */
    vnmo = ealloc1float(nvnmo);
    tnmo = ealloc1float(nvnmo);
    anis1 = ealloc1float(nvnmo);
    anis2 = ealloc1float(nvnmo);
    if (!getnparfloat(icdp+1,"vnmo",vnmo)) vnmo[0] = 2000.0;
    if (!getnparfloat(icdp+1,"tnmo",tnmo)) tnmo[0] = 0.0;
    if (!getnparfloat(icdp+1,"anis1",anis1)) 
      for (i=0; i<nvnmo; i++) anis1[i] = 0.0;
    if (!getnparfloat(icdp+1,"anis2",anis2))
      for (i=0; i<nvnmo; i++) anis2[i] = 0.0;
    for (it=1; it<ntnmo; ++it)
      if (tnmo[it]<=tnmo[it-1]){
	fprintf(stderr,"Error for #cdp  %d\n",icdp);
	err("tnmo values must increase monotonically");
      }
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,vnmo,vnmo[0],vnmo[nvnmo-1],1,&tn,&ovv[icdp][it]);
    
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,anis1,anis1[0],anis1[nanis1-1],1,&tn,&oa1[icdp][it]);
    
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,anis2,anis2[0],anis2[nanis2-1],1,&tn,&oa2[icdp][it]);
    
    free1float(vnmo);
    free1float(tnmo);
    free1float(anis1);
    free1float(anis2);
  }

  /* sort (by insertion) sloth and anis functions by increasing cdp */
  for (jcdp=1; jcdp<ncdp; ++jcdp) {
    acdp = cdp[jcdp];
    aovv = ovv[jcdp];
    aoa1 = oa1[jcdp];
    aoa2 = oa2[jcdp];
    for (icdp=jcdp-1; icdp>=0 && cdp[icdp]>acdp; --icdp) {
      cdp[icdp+1] = cdp[icdp];
      ovv[icdp+1] = ovv[icdp];
      oa1[icdp+1] = oa1[icdp];
      oa2[icdp+1] = oa2[icdp];
    }
    cdp[icdp+1] = acdp;
    ovv[icdp+1] = aovv;
    oa1[icdp+1] = aoa1;
    oa2[icdp+1] = aoa2;
  } 

  /* allocate workspace */
 
  oa1t = ealloc1float(nt);
  oa2t = ealloc1float(nt);

  if ((cdpt=ealloc1int(ntr+1))==NULL)
    fprintf(stderr,"***Sorry, space for cdpt could not be allocated\n");

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
  
  ntr = 0;
  do {
    ntr++;
    cdpt[ntr]=(int) tr.cdp;
    // fprintf(stderr,"sx=%d,tracl=%d,tracf=%d\n",tr.sx,tr.tracl,tr.tracf);
    efwrite(&tr,HDRBYTES,1,headerfp);
    efwrite(tr.data,FSIZE, nt, tracefp);    
  } while (gettr(&tr) && (ntr < ntrmax) );
  erewind(tracefp);
  erewind(headerfp);  
  fprintf(stderr,"ntr=%d\n",ntr);
 
  
  // Allocate memory for data and model
  
  if ((d=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");

  if ((csp=ealloc3float(nt,nhcsp,ncsp))==NULL)
    fprintf(stderr,"***Sorry, space for csp could not be allocated\n");
 
  if ((cspfold=(unsigned short ***) alloc3(nt,nhcsp,ncsp,sizeof(unsigned short)))==NULL) 
    fprintf(stderr,"***Sorry, space for csp could not be allocated\n");
  
  if ((m=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for m could not be allocated\n");
    
  if ((x=ealloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for x could not be allocated\n");

  if ((he=ealloc1float(nhcsp))==NULL)
    fprintf(stderr,"***Sorry, space for he could not be allocated\n"); 
  
  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n");

  if ((tc=ealloc1float(nhcsp))==NULL)
    fprintf(stderr,"***Sorry, space for tc could not be allocated\n");

  if ((hec=ealloc1float(nhcsp))==NULL)
    fprintf(stderr,"***Sorry, space for hec could not be allocated\n");
 
  if ((q=ealloc1float(nq+1))==NULL)
    fprintf(stderr,"***Sorry, space for q could not be allocated\n");

  if ((mm=alloc1float(nt*nq))==NULL)
    fprintf(stderr,"***Sorry, space for mm could not be allocated\n");

  if ((Wm=alloc1float(nt*nq))==NULL)
    fprintf(stderr,"***Sorry, space for Wm could not be allocated\n");

  if ((Wd=alloc1float(nt*nhcsp))==NULL)
    fprintf(stderr,"***Sorry, space for Wd could not be allocated\n");

  // Velocity axis for Radon space

  if (nq>1){
    dq = (qmax-qmin)/(nq-1);
    for (i=0;i<nq;i++) q[i] = qmin+i*dq;
  }
  else  q[0]=0;

  // arrays for nmo model
  // Velint will contain velocity law for one csp
  if ((velint=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Space for velint could not be allocated\n");
  // velatcsp will contain velocity law for ncsps 
  if ((velatcsp=ealloc2float(nt,ncsp))==NULL)
    fprintf(stderr,"Cannot allocate velatcsp\n");

  /* Create axis for output cpds, equivalent offset and time */
  //if (testhe==2) hemin=0;
  dx=dxcdp;
  dhe=(hemax-hemin)/(nhcsp-1);      
  for(ix=0;ix<nx;ix++) x[ix]=cdpmin+dx*ix;
  for(it=0;it<nt;it++) t[it]=0+it*dt;
  for(ih=0;ih<nhcsp;ih++) he[ih]=hemin+ih*dhe;
  if (verbose) fprintf(stderr,"hemin=%f,hemax=%f,dhe=%f\n",he[0],he[nhcsp-1],dhe);
  ///////////////////////////////////////////////////////////////////////
  //if (verbose) for (icdp=0;icdp<ncdp;icdp++) fprintf(stderr,"ovv[%d][100]=%f\n",icdp,ovv[icdp][100]);   

  fprintf(stderr,"Starting cdp loop ****\n");  

  /* Loop to compute ncsp CSP gathers from x[ix] to x[ix+ncsp] */
  for (ix=0;ix<nx;ix+=ncsp){
    fprintf(stderr,"cdp=%f to cdp=%f\n",x[ix],x[ix+ncsp-1]);   

    /* Initialize to zero the ouput trace and the csp gather corresp to ix */  
    for (icsp=0;icsp<ncsp;icsp++)
      for (ih=0;ih<nhcsp;ih++) 
	for (it=0;it<nt;it++){
	  csp[icsp][ih][it]=0; // CSP trace
	  cspfold[icsp][ih][it]=0;
	}


    /* Loop for mapping one trace to the CSP gather */
    for (itr=0;itr<ntr;itr++){     
      efread(tr.data, FSIZE, nt, tracefp);
      efread(&tr,HDRBYTES,1,headerfp);
      // Compute the contribution of the read trace to every csp
      for (icsp=0;(icsp<ncsp)&&(x[ix+icsp]<=cdpmax);icsp++){
	/* compute new square slowness and anis function */
	interpovv(nt,ncdp,cdp,ovv,oa1,oa2,x[ix+icsp],velint,oa1t,oa2t);
	memcpy((void *) velatcsp[icsp],
	       (const void *) velint,nt*sizeof(float));
	/* Map this trace to the current csp 
	Three options:
        equiv_offset_test does all computations, (precise and slow)
	equiv_offset_1side precomputes times, and uses 1 side of offset
        equiv_offset precomputes times, and uses + and - offsets */
	if (fabs(tr.cdp-x[ix+icsp])<aper){
	  if (testhe==1) 
	    equiv_offset_test(tr,csp[icsp],cspfold[icsp],t,x[ix+icsp],
			      velatcsp[icsp],dt,he,dhe,nt,nhcsp,beta,cdpspace,scalefold);
	  else if (testhe==2)
	    equiv_offset_1side(tr,csp[icsp],cspfold[icsp],t,x[ix+icsp],
			       velatcsp[icsp],dt,he,dhe,nt,nhcsp,tc,hec,
			       precise,beta,cdpspace,scalefold);  
	  else
	    equiv_offset(tr,csp[icsp],cspfold[icsp],t,x[ix+icsp],
			 velatcsp[icsp],dt,he,dhe,nt,nhcsp,tc,hec,precise,
			 beta,cdpspace,scalefold);
	}
      }
    }
    erewind(tracefp);
    erewind(headerfp);
    
    for (icsp=0;(icsp<ncsp)&&(x[ix+icsp]<=cdpmax);icsp++){
      /* compute new square slowness and anis function */
      // interpovv(nt,ncdp,cdp,ovv,oa1,oa2,x[ix+icsp],velint,oa1t,oa2t); 
      if (fold)
	for (ih=0;ih<nhcsp;ih++) 
	  for (it=0;it<nt;it++)
	    if (cspfold[icsp][ih][it]>0.5) 
	      csp[icsp][ih][it]/=(float) cspfold[icsp][ih][it];

      if (output==4){
	for (itr=0;itr<nhcsp;itr++){
	  memcpy((void *) tr.data,
		 (const void *) csp[icsp][itr],nt*sizeof(float));
	  // for (it=0;it<nt;it++) tr.data[it]=csp[icsp][itr][it];
	  tr.cdp=(int) x[ix+icsp]; // front of the CSP
	  tr.dt=(int) (dt*1e6);       
	  tr.ntr=nhcsp;
	  tr.ns=nt;
	  tr.tracl=itr+1;
	  tr.tracr=itr+1;
	  tr.offset=(int) (2*he[itr]);
	  tr.sx=(int) x[ix];
	  tr.gx=(int) x[ix];
	  fputtr(cspfilep,&tr);    
	}
      }
      
      if (output==3 || output==4 || output==5){
	if (sinc){
	  eomigls_sparse_sinc(m,mm,t,he,q,velatcsp[icsp],csp[icsp],nt,nhcsp,nq,fmax,step,eps,eps1,eps2,itercg,iter_end,smooth,testadj,buffer,pervmin,dperv,norm,parmute,mute,t0,Wm,Wd,centralq);
	}
	else{
	  if ( (nt*nhcsp) < USHRT_MAX && (nt*nq) < USHRT_MAX ) 
	    eomigls_sparse(m,mm,t,he,q,velatcsp[icsp],csp[icsp],nt,nhcsp,nq,fmax,step,eps,eps1,eps2,itercg,iter_end,smooth,testadj,buffer,pervmin,dperv,norm,parmute,mute,t0,Wm,Wd,centralq);
	  else
	    eomigls_sparse_large(m,mm,t,he,q,velatcsp[icsp],csp[icsp],nt,nhcsp,nq,fmax,step,eps,eps1,eps2,itercg,iter_end,smooth,testadj,buffer,pervmin,dperv,norm,parmute,mute,t0,Wm,Wd,centralq);
	}
      }

      if (output==6 || output ==7 ){ 
	if ( (nt*nhcsp) < USHRT_MAX && (nt*nq) < USHRT_MAX ) 
	  eomgetvel(m,mm,t,he,q,velatcsp[icsp],csp[icsp],nt,nhcsp,nq,fmax,step,eps,eps1,eps2,itercg,iter_end,smooth,testadj,buffer,pervmin,dperv,norm,itervel,t0,Wm,Wd);    
	else  
	  eomgetvel_large(m,mm,t,he,q,velatcsp[icsp],csp[icsp],nt,nhcsp,nq,fmax,step,eps,eps1,eps2,itercg,iter_end,smooth,testadj,buffer,pervmin,dperv,norm,itervel,t0,Wm,Wd);        
      }
      // Output trace
     
      if (output==2 || output==4){
	for (itr=0;itr<nhcsp;itr++){
	  memcpy((void *) tr.data,
		 (const void *) csp[icsp][itr],nt*sizeof(float));
	  // for (it=0;it<nt;it++) tr.data[it]=csp[icsp][itr][it];
	  tr.cdp=(int) x[ix+icsp]; // front of the CSP
	  tr.dt=(int) (dt*1e6);       
	  tr.ntr=nhcsp;
	  tr.ns=nt;
	  tr.tracl=itr+1;
	  tr.tracr=itr+1;
	  tr.offset=(int) (2*he[itr]);
	  tr.sx=(int) x[ix];
	  tr.gx=(int) x[ix];
	  puttr(&tr);    
	}
      }
      fprintf(stderr,"output=%d\n",output);
      if (output==5 || output ==7 || output ==4){
	for (itr=0;itr<nq;itr++){
	  for (it=0;it<nt;it++) tr.data[it]=mm[itr*nt+it];
	  tr.cdp=(int) x[ix+icsp]; // front of the CSP
	  tr.dt=(int) (dt*1e6);       
	  tr.ntr=nq;
	  tr.ns=nt;
	  tr.tracl=itr+1;
	  tr.tracr=itr+1;
	  tr.f2=q[itr];
	  tr.sx=(int) x[ix];
	  tr.gx=(int) x[ix];
	  fprintf(stderr,"writing to file\n");
	  fputtr(modelfilep,&tr);    
	}
      }
    
      if (((output==1)||(output==3)||(output==6))&&(x[ix+icsp]<=cdpmax)){
	eomig1(csp[icsp],he,m,t,velatcsp[icsp],t0);
	for (it=0;it<nt;it++)
	  m[it]/=(velatcsp[icsp][0]*sqrt(pi*MAX(t[it],1e-2)));

	memcpy(tr.data, m, nt*sizeof(float));

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
  free2float(oa1);
  free2float(oa2);
  free1float(oa1t);
  free1float(oa2t);
  free1float(Wd);
  free1float(Wm);
  free1float(mm);
  free1float(q);
  free2float(velatcsp);
  free1float(velint);
  free1float(hec);
  free1float(tc);  
  free1float(t);
  free1float(he);
  free1float(x);
  free1float(m);

  free3((void ***) cspfold);

  free3float(csp);
  free1float(d);

  efclose(headerfp);
  efclose(tracefp);
  efclose(modelfilep);
  efclose(cspfilep);

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
































