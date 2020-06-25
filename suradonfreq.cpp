/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUEOMIG:  $Date: March 1999  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "radonfrequency.h"
#include <signal.h>
#include <math.h>
#include <time.h>


/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADONFREQ - Frequency domain Radon Transform                      ",
  "                                                                     ",
  " 	   								",
  " suradonfreq < stdin > stdout [optional parameters]       		",
  "                                                                     ",
  " IMPORTANT:****************************************                  ",
  " Input must be sort by cdp and offset, for example                   ",
  "                   susort cdp offset < input | .....                 ",
  "                                                                     ",
  " It removes This program is intended to input a line of data sort    ",
  " by cdps and, every cdp by offset.                                   ",
  "                                                                     ",
  " Notes:                                                              ",
  " Most parameters have reasonable default values. The most important  ",
  " parameters to set (for which there is no default for every data set ",
  " are: cdpmin, cdpmax, par and qmin.                                  ",
  " solver definers whether to have high resolution RT (cgfft)          ",
  " or standard RT (toep).                                              ",
  " rtmethod defines if working with linear (1) , parabolic (2) o       ",
  " pseudohyperbolic RT (3). In this last case need to set also depth   ",
  " (see Foster and Mosher,)                                            ",
  "                                                                     ",
  "                                                                     ",
  "                                                                     ",
  " cdpmin=0        Fisrt CDP in meters                                 ",
  " cdpmax=100      Last CDP 						",
  " dxcdp=1         Interval of cdps to process                         ",
  " par=            file with cdp numeber, stacking velocities and time ",
  "                 as obtained  from Velan (standard PARFILE in SU)    ",
  " rtmethod=2      shape of integration surface                        ",
  "                 =1 linear                                           ",
  "                 =2 parabolic                                        ",
  "                 =3 Pseudohyperbolic (Foster and Mosher)             ",
  " solver=cgfft    Method to solve the system of equations             ",
  "                 =cgfft very fast cg with fft                        ",
  "                 =toep  Toeplitz                                     ",
  "                 =adj Simple adjoint                                 ",
 " nhmax=200       Maximum number of traces for cdp                     ",
  " qmin= -5e-9     Minimum q parameter for the RT space                ", 
  " nq=100          number of traces in the RT space                    ",
  " itercg=50       Internal iterations for CG                          ",
  " iter_end=3      External iterations for IRLS                        ",
  " eps1=5e-1       Numerator hyperparameter for Wd (quantil of data)   ",
  " eps2=5e-1       Denominator  hyperparameter for Wm (quantil of model)",
  " fmax=0.8/(2*dt) Maximum frequency to preserve in the transform      ",
  "                 =0 Nyquist (1/2*dt)                                 ", 
  " norm=0          model weights derivated from Cauchy norm            ",
  "                 =1 model weights derivated from L1 norm             ", 
  " mute=0          =1 mute multiples                                   ",
  " verbose=0       =1  verbose output                                  ",
  " t0=0            First useful time                                   ",
  " smute=2         stretch greater than smutex100 % is muted		",
  " nmofactor=1     nmofactor * offset is used for NMO	                ",
  " quantil=1       filter for large data outliers                      ",
  "                                                                     ",
  " mute parameters :                                                   ",
  " tmin_m=0        min time for pass mute window                       ",
  " tmax_m=nt*dt    max time for pass mute window                       ",
  " ihmin_m=nq/2+2  first trace for pass mute window                    ",
  " ihmax_m=nq      last trace number for pass mute window              ",
  " thres_m=0.2    values less than threshold are removed inside window ",
  " slope_m=3      slope in pass mute window (ih top / ih bottom)       ",
  "                                                                    ",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset
 **************** end self doc ***********************************/
static void closefiles(void);
/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
FILE *modelfilep;       /* fp for model file output  		*/
///////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  segy tr,tr2;
  cwp_String modelfile=""; /* output sufile for the model */ 
  cwp_String cmpfile="";   /* output sufile for the cmp   */ 
  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;
  
  // Data space
  int nt;   // number of time samples 
  int nx;    // number of midpoints
  int nhmax; // max number of offsets.
  float dt; // Time sampling
  float dx;  // midpoint interval
  float eps; // small tolerance number
  
  ///////////////////////////////////////////////////////////////
  int ih, ix; // General counters 
  register int it;
  float *d;      /* single trace */
  float **cmp;   /* Common midpoint gather */
  float *x;      /* axis for CDPs */
  float *h;      /* axis for input offset */
  float *t;     // time axis for input and output
  float t0;      // First useful time

  unsigned int ntrmax;  // Numebr of traces to process from data file
  /* Velocity */
  int ncdp;	/* number of cdps specified */
  float *cdp;	/* array[ncdp] of cdps */
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
  float fmax;
  // Radon space  
  int nq;  
  float qmin;
  float *q;
  float **mm;            // Temporal array for the model
  
  //  LS mig  
  float step;
  float eps1;
  float eps2;
  int itercg;
  int iter_end;
  int norm;
  float *Wm;   // Factorization of inv model covariance matrix 
  float factor;
  float smute;
  float nmofactor;
  int *nhcmp;
  float quantil;
  int ngettr; // Number of bytes read by gettr (0 after)
  float depth;   // depth for pseudohyperbolic
  int rtmethod;
  char *solver;
  int keepwm; // If keepwm==1 the value for Wm from previous CMP is used  
  
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
  fprintf(stderr,"nt=%d,dt=%f\n",nt,dt);  
  
  // Get parameters 
  
  if (!getparfloat("cdpmin", &cdpmin))  cdpmin = 0;
  if (!getparfloat("cdpmax", &cdpmax))  cdpmax = 1000;
  if (!getparfloat("dxcdp", &dxcdp))  dxcdp = 1;
  if (!getparint("nhmax", &nhmax))  nhmax = 300;
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparfloat("eps", &eps))  eps = 1e-7;
  if (!getparfloat("eps1", &eps1))  eps1 = 5e-1;
  if (!getparfloat("eps2", &eps2))  eps2 = 5e-1;
  if (!getparfloat("step", &step))  step = 0.98;
  if (!getparint("itercg", &itercg))  itercg = 50;
  if (!getparint("iter_end", &iter_end))  iter_end = 3;
  if (!getparfloat("qmin", &qmin))  qmin = -5e-9;
  if (!getparint("nq", &nq))  nq = 100;
  if (!getparint("norm",&norm)) norm=0;
  if (!getparuint("ntrmax",&ntrmax)) ntrmax = 1000000;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparint("keepwm",&keepwm)) keepwm=0;
  if (!getparfloat("factor",&factor)) factor=0.8;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=1;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
  if (!getparstring("cmpfile",&cmpfile)) cmpfile="cmpfile.su";
  if (!getparfloat("quantil",&quantil)) quantil=1;
  if (!getparint("rtmethod",&rtmethod)) rtmethod=2;
  if (!getparfloat("depth",&depth)) depth=2000;
  if (!getparfloat("fmax",&fmax)) fmax=0.8/(2*dt);
  if (!getparstring("solver",&solver)) solver="cgfft";
  if (!getparfloat("quantil",&quantil)) quantil=1;
  /* Format the string solver to have the same length */ 
  if (STREQ(solver,"toep")) solver="toep__";
  else if (STREQ(solver,"cgfft")) solver="cgfft_";
  else if (STREQ(solver,"adj")) solver="adj___";

  if (fmax==0) fmax=1./(2*dt);
  



  /***********************************************************************/  
  /* for mute in the migrated space we need to set the geometry of the mask */
  int mute; 
  mutemask_par par;

  if (!getparint("mute",&mute)) mute = 0;
  if (!getparfloat("tmin_m",&par.tmin)) par.tmin = 0.; 
  if (!getparfloat("tmax_m",&par.tmax)) par.tmax = nt*dt;
  if (!getparint("ihmin_m",&par.ihmin)) par.ihmin = (int) (nq/2+2);
  if (!getparint("ihmax_m",&par.ihmax)) par.ihmax = (int) (nq);
  if (!getparfloat("slope_m",&par.slope)) par.slope = 3;     
  if (!getparfloat("thres_m",&par.threshold)) par.threshold = 0.2;     
  /* Write the new mute parameters to a file for plots */
  write_curve_mute(par);
  float **M=ealloc2float(nt,nq);
  if (mute) mutemask2(M,nq,nt,dt,par); // data independent mute
  /************************************************************************/  
  modelfilep=efopen(modelfile,"w");
  nx=(int) floor((cdpmax-cdpmin)/dxcdp + 1);

  if (verbose) fprintf(stderr,"number of ouput cdps = %d\n",nx);

  /* get velocity functions, linearly interpolated in time */
  // Velint will contain velocity law for one csp
  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  getvelocities(dt,nt,ncdp,cdp,ovv);
  
  /* Look for user-supplied tmpdir */

  if (!getparstring("tmpdir",&tmpdir) &&
      !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
  if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
    err("you can't write in %s (or it doesn't exist)", tmpdir);  
 
  // Store CMP traces in the disk for later use
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

  mm=ealloc2float(nt,nq);
  cmp=ealloc2float(nt,nhmax);
  d=ealloc1float(nt);
  x=ealloc1float(nx);
  h=ealloc1float(nhmax);
  t=ealloc1float(nt);
  q=ealloc1float(nq);
  Wm=ealloc1float(nt*nq);
  nhcmp=ealloc1int(nx);

  q[0]=qmin;

  /* Create axis for output cpds, equivalent offset and time */
  dx=dxcdp;
  for(ix=0;ix<nx;ix++) x[ix]=cdpmin+dx*ix;
  for(it=0;it<nt;it++) t[it]=t0+it*dt;
  
  ///////////////////////////////////////////////////////////////////////

  fprintf(stderr,"Starting cdp loop ****\n");

  /*************************************************************/
  /* Loop on the CMP gathers  */
  for (ix=0;ix<nx;ix++){
    if (verbose) fprintf(stderr,"cdp=%f\n",x[ix]);   
    /* put back last trace */
    tr=tr2;    
    // search for the first cdp  
    while (tr.cdp<x[ix]) if (!(ngettr=gettr(&tr))) break;
    
    if (verbose) fprintf(stderr,"tr.cdp=%d\n",tr.cdp);

    /* Initialize to zero the ouput trace and the cmp gather corresp to ix */  
    for (ih=0;ih<nhmax;ih++) memset(cmp[ih],(int) '\0',nt*FSIZE);

    erewind(tracefp);
    erewind(headerfp);   
    /**************************************************/    
    // Save the current cdp only (x[ix]) into a temporal file
    ntr = 0;
    
    do {
      ntr++;
      if (verbose) fprintf(stderr,"ntr=%d,tr.cdp=%d,tr.offset=%d\n",
			   ntr,(int) tr.cdp,(int) tr.offset);
      efwrite(&tr,HDRBYTES,1,headerfp);
      efwrite(tr.data,FSIZE, nt, tracefp);
    } while ( (ngettr=gettr(&tr)) && (tr.cdp == x[ix]));
    fprintf(stderr,"x=%f, ntr=%d\n",x[ix],ntr);
    /* Save last trace for next cdp */
    tr2=tr;    
    /**************************************************/
    // Read the current cdp
    erewind(tracefp);
    erewind(headerfp);  
    
    if (ntr<3) continue; // empty cdp
    for(ih=0;ih<ntr;ih++){
      efread(tr.data, FSIZE, nt, tracefp);
      efread(&tr,HDRBYTES,1,headerfp);
      h[ih]=tr.offset;
      memcpy(cmp[ih],tr.data,nt*FSIZE);
    }
    nhcmp[ix]=ntr;

    ///////////////////////////////////////////////////
    /* compute new square slowness  */
    interpovv(nt,ncdp,cdp,ovv,x[ix],velint);

    radonsolver0_subtract(cmp,h,nhcmp[ix],t,nt,dt,mm,q,nq,velint,itercg,iter_end,step,
		      eps2,eps1, quantil,norm,factor,smute,nmofactor,rtmethod,depth,
		      fmax,solver,M,mute); 

    // Output trace
    //erewind(tracefp);
    erewind(headerfp);  
    
    for (ih=0;ih<nhcmp[ix];ih++){ 
      efread(&tr,HDRBYTES,1,headerfp);
      memcpy((void *) tr.data,(const void *) cmp[ih],nt*sizeof(float));
      if (verbose) fprintf(stderr,"Original data ===> tr.offset=%d\n",tr.offset);
      puttr(&tr);
    }

    fprintf(stderr,"nhcmp[%d]=%d\n",ix,nhcmp[ix]);    

    for (itr=0;itr<nq;itr++){
      for (it=0;it<nt;it++) tr.data[it]=mm[itr][it];
      tr.cdp=(int) x[ix]; // front of the CMP
      tr.dt=(int) (dt*1e6);       
      tr.ntr=nq;
      tr.ns=nt;
      tr.tracl=itr+1;
      tr.tracr=itr+1;
      tr.f2=q[itr];
      tr.sx=(int) x[ix];
      tr.gx=(int) x[ix];
      fputtr(modelfilep,&tr);    
    }
  }

  free2float(M);
  free1float(cdp);
  free2float(ovv);
  free1float(Wm);
  free2float(mm);
  free1float(q);
  free1float(velint);
  free1float(t);
  free1float(h);
  free1float(x);
  free1int(nhcmp);
  free2float(cmp);
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



void write_curve_mute(mutemask_par par)
{
  
  FILE *fp;
  fp=efopen("curve1","w");
  
  fprintf(fp,"%f %f\n",par.tmin, par.ihmin);
  fprintf(fp,"%f %f\n",par.tmin, par.ihmax);
  fprintf(fp,"%f %f\n",par.tmax, par.ihmax);
  fprintf(fp,"%f %f\n",par.tmax, par.ihmin);
  fprintf(fp,"%f %f\n",par.tmin, par.ihmin);

  efclose(fp);
  
  return;
}



























