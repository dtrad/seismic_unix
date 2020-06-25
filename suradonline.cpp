/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUEOMIG:  $Date: March 1999  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "radonsolver.h"
#include <signal.h>
#include <math.h>
#include <time.h>


/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADONLINE - Frequency domain Parabolic Radon Transform            ",
  "                                                                     ",
  " 	   								",
  " suradonline < stdin > stdout [optional parameters]       		",
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
  " itercg=nq/5     Internal iterations for CG                          ",
  " iter_end=3      External iterations for CG                          ",
  " eps1=1e-3       Numerator hyperparameter for Wm                     ",
  " eps2=1e-3       Denominator  hyperparameter for Wm                  ",
  " mute=0          =1 mute multiples                                   ",
  " parnmute=1e-8   if q > parmute the radon trace is tapered           ",
  " verbose=0       =1  verbose output                                  ",
  " t0=0            First useful time                                   ",
  " quantil=1       Quantil to scale the model weight function          ",
  " smute=2         stretch greater than smutex100 % is muted		",
  " nmofactor=1     nmofactor * offset is used for nmo		        ",
  " pseudohyp=0     =1 Foster and Mosher                                ",
  "                                                                     ",
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
FILE *cmpfilep;         /* fp for cmp file output  		*/
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

  unsigned int ntrmax;  // Numebr of traces to process from data file
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
  float fmax;
  // Radon space  
  int nq;  
  float qmin;
  float *q;
  float dq;
  float **mm;            // Temporal array for the model
  
  //  LS mig  
  float step;
  float eps1;
  float eps2;
  int itercg;
  int iter_end;
  int norm;
  float *Wm;   // Factorization of inv model covariance matrix 
  float *Wd;   // Factorization of inv residual covariance matrix 
  float parmute;
  int mute;
  float t0mute;
  float factor;
  float smute;
  float nmofactor;
  int *nhcmp;
  float quantil;
  int ngettr; // Number of bytes read by gettr (0 after)
  int pseudohyp; // Apply pseudohyperbolic RT
  float depth;   // depth for pseudohyperbolic
  int rtmethod;
  char *solver;
  int keepwm; // If keepwm==1 the value for Wm from previous CMP is used  
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
  fprintf(stderr,"nt=%d,dt=%f\n",nt,dt);  
  
  // Get parameters 
  
  if (!getparfloat("cdpmin", &cdpmin))  cdpmin = 0;
  if (!getparfloat("cdpmax", &cdpmax))  cdpmax = 1000;
  if (!getparfloat("dxcdp", &dxcdp))  dxcdp = 1;
  if (!getparint("nhmax", &nhmax))  nhmax = 300;
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparfloat("eps", &eps))  eps = 1e-7;
  if (!getparfloat("eps1", &eps1))  eps1 = 5e-2;
  if (!getparfloat("eps2", &eps2))  eps2 = 5e-2;
  if (!getparfloat("step", &step))  step = 0.9;
  if (!getparint("itercg", &itercg))  itercg = 5;
  if (!getparint("iter_end", &iter_end))  iter_end = 1;
  if (!getparfloat("qmin", &qmin))  qmin = -5e-9;
  if (!getparint("nq", &nq))  nq = 100;
  if (!getparint("norm",&norm)) norm=0;
  if (!getparuint("ntrmax",&ntrmax)) ntrmax = 1000000;
  if (!getparfloat("parmute",&parmute)) parmute=1e-8;
  if (!getparint("mute",&mute)) mute=0;
  if (!getparfloat("t0mute",&t0mute)) t0mute=0;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparint("keepwm",&keepwm)) keepwm=0;
  if (!getparfloat("factor",&factor)) factor=0.8;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=1;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
  if (!getparstring("cmpfile",&cmpfile)) cmpfile="cmpfile.su";
  if (!getparfloat("quantil",&quantil)) quantil=1;
  if (!getparint("pseudohyp",&pseudohyp)) pseudohyp=0;
  if (!getparint("rtmethod",&rtmethod)) rtmethod=2;
  if (!getparfloat("depth",&depth)) depth=2000;
  if (!getparfloat("fmax",&fmax)) fmax=0.8/(2*dt);
  if (!getparstring("solver",&solver)) solver="toep__";
  if (!getparfloat("quantil",&quantil)) quantil=1;
  /* Format the string solver to have the same length */ 
  if (STREQ(solver,"toep")) solver="toep__";
  else if (STREQ(solver,"cgfft")) solver="cgfft_";
  else if (STREQ(solver,"adj")) solver="adj___";

  if (fmax==0) fmax=1./(2*dt);
  
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
  
  modelfilep=efopen(modelfile,"w");
  cmpfilep=efopen(cmpfile,"w");
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
			   ntr,tr.cdp,tr.offset);
      efwrite(&tr,HDRBYTES,1,headerfp);
      efwrite(tr.data,FSIZE, nt, tracefp);
    } while ( (ngettr=gettr(&tr)) && (tr.cdp == x[ix]));
    if (verbose) fprintf(stderr,"x=%f, ntr=%d\n",x[ix],ntr);
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

    // Check the cmp array
    if (0){
      erewind(headerfp);  
      for (itr=0;itr<nh;itr++){
	  efread(&tr,HDRBYTES,1,headerfp);
	  memcpy(tr.data,cmp[itr],nt*FSIZE);
	  fputtr(cmpfilep,&tr);
      }
    }

    radonsolver0_mute(cmp,h,nhcmp[ix],t,nt,dt,mm,q,nq,velint,itercg,iter_end,step,
		      eps2,eps1, quantil,norm,factor,smute,nmofactor,rtmethod,depth,
		      fmax,solver,0,parmute,0,mute,t0mute); 

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

    if (0){
      save_gather(mm,nq,nt,0.004,"model2");
      system("sufft  < model2 | suamp | suxwigb title=\"model2\" perc=99 & ");
    }


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

  free1float(amps);
  free1float(ffilter);
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
  efclose(cmpfilep);

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






























