/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUEOMIG:  $Date: March 1999  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "radontd_win2.h"
#include <signal.h>
#include <math.h>
#include <time.h>


/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADONLINETD - Time domain hyperbolic Radon Transform              ",
  "                                                                     ",
  " 	   								",
  " suradonlinetd < stdin > stdout [optional parameters]       		",
  "                                                                     ",
  " IMPORTANT:****************************************                  ",
  " Input must be sort by cdp and offset, for example                   ",
  "                   susort cdp offset < input | .....                 ",
  "                                                                     ",
  " This program removes multiples from a whole line of data.           ",
  "                                                                     ",
  " cdpmin=1        Fisrt CDP in meters                                 ",
  " cdpmax=1        Last CDP 						",
  " par=            file with stacking velocities and time as obtained  ",
  "                 from Velan (PARFILE)                                ",
  " modelfile=model.su  filename for the Radon transform                ", 
  " nhmax=200       Maximum number of traces for cmp                    ",
  " pervmin = 4e-9      minimum perturbation on velocity                ",
  " dperv=0             defines the rate of incresing velocity spacing  ",
  " dq=2e-8             defines the rate of constant slowness spacing   ",
  " nq=100          Size of Radon EOM                                   ",
  " itercg=10       Internal iterations for CG                          ",
  " iter_end=1      External iterations for CG                          ",
  " eps1=1e-3       Numerator hyperparameter for Wm                     ",
  " eps2=1e-3       Denominator  hyperparameter for Wm                  ",
  " mute=0          =1 mute multiples                                   ",
  " parnmute=1e-8   if q > parmute the radon trace is tapered           ",
  " verbose=0       =1  verbose output                                  ",
  " t0=0            First useful time                                   ",
  " quantil=1       Quantil to scale the model weight function          ",
  " tm=0            The time for the first  multiple to be subtracted   ",
  " plot=0          =1 =>      plot Radon space an primaries            ",
  "                 =2 =>      plots after mute                         ",
  "                                                                     ",
  " centralq=0      the q trace where the velocity law lays             ",
  " dataprec =0     event with energy higher than quantil 99 are        ",
  "                 downweighted by Q50/Q99                             ",
  " outputmodel=0   =1 output modelfile                                 ",
  " LI=0            =1 linear interpolation                             ",
  " smute=0         =1 stretching mute                                  ",
  " typewav=0       =1 include wavelet convolution and correlation      ",
  "                                                                     ",
  "                                                                     ",
  " Example: 		                                                ",
  "                                                                     ",
  " suradonlinetd < data.su par=stkvel.data.su itercg=25 iter_end=3     ",
  "                 modelfile=modelsp.su outputmodel=1 > dataout.su     ",
  "                                                                     ",
  "                                                                     ",
  "  suradonlinetd < filein.su  cdpmin=1    cdpmax=1000 dxcdp=1         ",
  "  nhcdp=200 par=stkvel.filein.su itercg=25 iter_end=1 eps2=1e-1      ",
  "  eps1=$eps1 nq=$nq verbose=$verbose ntrmax=$ntrmax tmpdir='./'      ",
  "  norm=0 t0=$t0 mute=$mute parmute=$parmute LI=$LI nw=$nw            ",
  "  centralq=$centralq dperv=$dperv pervmin=$pervmin filtout=0         ",
  "  restart=1  tm=$tm plot=1  outputmodel=1 smute=3 typewav=$typewav   ",
  "     	fpeak=$fpeak > $FILEREC                                 ", 
  "                                                                     ",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset */
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
  inv_par inv;
  segy tr,tr2;
  cwp_String modelfile=""; /* output sufile for the model */ 
  cwp_String cmpfile="";   /* output sufile for the cmp   */ 
  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;
  
  // Data space
  int nh;    // number of offset traces
  int nt;    // number of time samples 
  int nx;    // number of midpoints
  int *nhcmp;// number of offset traces for each CMP
  int nhmax; // max number of offsets.
  float dt;  // Time sampling
  float dx;  // midpoint interval
  int ngettr;// Number of bytes read by gettr (0 after)
  
  ///////////////////////////////////////////////////////////////
  int  k, ih, ix; // General counters 
  register int it;
  float *d;      // single trace 
  float **cmp;   // Common midpoint gather 
  float *x;      // axis for CDPs 
  float *h;      // axis for input offset 
  float *t;      // time axis for input and output
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
  // Radon space  
  int nq;  
  float qmin;
  float *q;
  float **mm;            // Temporal array for the model
  
  // CG  
  float *Wm;   // Factorization of inv model covariance matrix 
  //float *Wd;   // Factorization of inv residual covariance matrix 
  float parmute;
  int mute;
  float factor;
  float quantil;
  float smute;

  int rtmethod; // =3 HRT, =2 PRT, =1 LRT
  int keepwm; // If keepwm==1 the value for Wm from previous CMP is used  
  float pervmin; // minimum perturbation on velocity        
  float dperv;   // defines the rate of incresing velocity spacing
  int centralq;  // q trace with the velocity trend 
  int dataprec;  // dataprec =0      event with energy higher than quantil 99 are  
                 // downweighted by Q50/Q99                           
  // Wavelet for the RT operator
  int nw;        // number of point for the wavelet
  float fpeak;   // peak frequency for the wavelet
  int typewav;   // type of wavelet
  int LI;        // Linear interpolation 
  int nreg;      // defines how many traces besides the velocity trend have reg sampling
  int taperflag; 
  int smooth;  
  int plot;
  int outputmodel;
  float tm;
  int itm;
  // Filter
  float *ffilter=NULL;
  float *amps=NULL;
  
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
  
  if (!getparfloat("cdpmin", &cdpmin))  cdpmin = 1;
  if (!getparfloat("cdpmax", &cdpmax))  cdpmax = 1;
  if (!getparfloat("dxcdp", &dxcdp))  dxcdp = 1;
  if (!getparint("nhmax", &nhmax))  nhmax = 200;
  if (!getparint("verbose", &verbose)) verbose = 0; 
  // inv_par elements 
  if (!getparfloat("eps1", &inv.eps1))  inv.eps1 = 1e-1;
  if (!getparfloat("eps2", &inv.eps2))  inv.eps2 = 1e-1;
  if (!getparfloat("eps", &inv.eps))  inv.eps = 1e-7;
  if (!getparfloat("step", &inv.step))  inv.step =0.99; 
  if (!getparint("itercg", &inv.itercg))  inv.itercg = 10;
  if (!getparint("iter_end", &inv.iter_end))  inv.iter_end = 1;
  if (!getparint("restart",&inv.restart)) inv.restart = 1;
  if (!getparint("norm",&inv.norm)) inv.norm=0;
  if (!getparfloat("qmin", &qmin))  qmin = -5e-9;
  if (!getparint("nq", &nq))  nq = 100;

  if (!getparuint("ntrmax",&ntrmax)) ntrmax = 1000000;
  if (!getparfloat("parmute",&parmute)) parmute=1e-8;
  if (!getparint("mute",&mute)) mute=0;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparint("keepwm",&keepwm)) keepwm=0;
  if (!getparfloat("factor",&factor)) factor=0.8;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
  if (!getparstring("cmpfile",&cmpfile)) cmpfile="cmpfile.su";
  if (!getparfloat("quantil",&quantil)) quantil=1;
  if (!getparint("rtmethod",&rtmethod)) rtmethod=3; // HRT default
  if (!getparfloat("quantil",&quantil)) quantil=1;

  if (!getparfloat("pervmin",&pervmin)) pervmin = 40e-10;
  if (!getparfloat("dperv",&dperv)) dperv =0.;
  if (!getparfloat("factor", &factor))  factor =0.8;
  if (!getparint("taperflag", &taperflag))  taperflag =0;
  if (!getparint("smooth",&smooth)) smooth=0;
  if (!getparint("centralq",&centralq)) centralq=15;
  if (!getparint("dataprec",&dataprec)) dataprec =0;
  if (!getparint("nw",&nw)) nw =0;
  if (!getparfloat("fpeak",&fpeak)) fpeak =25;
  if (!getparint("typewav",&typewav)) typewav = 1;
  if (!getparint("LI",&LI)) LI = 0;
  if (!getparint("nreg",&nreg)) nreg = 5;
  if (!getparint("plot",&plot)) plot = 1;
  if (!getparfloat("tm",&tm)) tm=0;
  if (!getparfloat("smute",&smute)) smute=0;
  if (!getparint("outputmodel",&outputmodel)) outputmodel = 0;

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
  // Time for the first multiple in index units
  itm=(int) (tm/dt+0.5);
  
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


  // Initially define the time axis from zero but changed later to t[it0]=t0
  for(it=0;it<nt;it++) t[it]=0+it*dt;  
  
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
			   (int) ntr,tr.cdp,tr.offset);
      efwrite(&tr,HDRBYTES,1,headerfp);
      efwrite(tr.data,FSIZE, nt, tracefp);
    } while ( (ngettr=gettr(&tr)) && (tr.cdp == x[ix]));
    fprintf(stderr,"x=%f, ntr=%d\n",x[ix],(int) ntr);
    /* Save last trace for next cdp (it was not accepted in the while loop)*/
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
    TRACE;
    ///////////////////////////////////////////////////
    /* compute new square slowness  */
    interpovv(nt,ncdp,cdp,ovv,x[ix],velint);
    TRACE;
    // Check the cmp array
    if (0){
      erewind(headerfp);  
      for (itr=0;itr<nh;itr++){
	  efread(&tr,HDRBYTES,1,headerfp);
	  memcpy(tr.data,cmp[itr],nt*FSIZE);
	  fputtr(cmpfilep,&tr);
      }
    }
    TRACE;
    if (smute) smute_gather(cmp,nt,nhcmp[ix],t,h,velint,smute);
    TRACE;
    radontd_sparse(t,q,h,mm,cmp,nt,nhcmp[ix],nq,dt,velint,dperv,pervmin,t0,inv,
		   centralq,dataprec,nw,fpeak,typewav,LI,nreg,parmute,mute,itm,plot);
    // Output trace
    //erewind(tracefp);
    erewind(headerfp);  
    TRACE;
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
      if (outputmodel) fputtr(modelfilep,&tr);    
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






























