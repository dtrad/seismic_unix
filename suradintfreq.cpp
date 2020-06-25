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
#include "Complex.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADFREQ - Frequency domain Parabolic Radon Transform              ",
  "             (for interpolation)                                     ",
  " 	   								",
  " suradfreq < stdin > stdout [optional parameters]          		",
  "                                                                     ",
  " Input must be sort by cdp and offset, for example                   ",
  "                   susort cdp offset < input | .....                 ",
  "                                                                     ",
  " The output contains original traces at the given offsets if they   	",
  " are present or interpolated traces if the offsets are empty (gaps)  ",
  " This program is intended to input a whole line of data              ",
  " and produce a new line without gaps and regular offsets             ",
  "                                                                     ",
  " cdpmin=0        Fisrt CDP in meters                                 ",
  " cdpmax=100      Last CDP 						",
  " par=            file with stacking velocities and time as obtained  ",
  "                 from Velan (PARFILE)                                ",
  " offsetfile      ascii file that contains the desired output offset  ",
  "                 IMPORTANT: nhe has to be larger than the number of  ",
  "                 new offsets.                                        ",
  "                                                                     ", 
  " nhe=100         Maximum number of traces for cmp                    ",
  " ncmp=10         Number of CMP to compute for every file  reading    ",
  " nq=80           Size of Radon EOM                                   ",
  " testadj=0       =1 Test adjoint with random numbers                 ",
  " smooth=0        =1 Pass triangular filter to smooth  output         ",
  " itercg=5        Internal iterations for CG                          ",
  " iter_end=1      External iterations for CG                          ",
  " eps1=1e-3       Numerator hyperparameter for Wm                     ",
  " eps2=1e-3       Denominator  hyperparameter for Wm                  ",
  " mute=0          =1 mute nmute traces at the sides of Velocity space ",
  " parnmute=1      if q > parmute the radon trace is tapered           ",
  " verbose=0       =1  Extra Information                               ",
  " t0=0            First useful time                                   ",
  " smute=2         stretch greater than smutex100 % is muted		",
  " nmofactor=1.9   nmofactor * equiv offset is used for nmo		",
  " lsmethod=0      =0 Least Squares Radon with toeplitz                ",
  "                 =1 High resolution with WTLS                        ",
  " under=0.999     Tolerance to match offsets (from below)             ",
  " over =1.001     Tolerance to match offsets (from above)             ",
  "                                                                     ",
  "                 if under * h_0rig < h_new < over * h_Orig           ",
  "                 then original trace will be output                  ",
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
FILE *offsetfile; 
///////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  segy tr,trf,traux;
  cwp_String modelfile=""; /* output sufile for the model */ 
  cwp_String cmpfile="";   /* output sufile for the cmp   */ 
  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;

  int nh;   // number of offset traces
  int nt;   // number of time samples 
  int nx;    // number of midpoints
  int nhe; // number of offsets.
  int nh2;   // number of output offset traces
  float dt; // Time sampling
  float dh; // offset interval
  float dhe; // equivalent offset interval
  float dx;  // midpoint interval
  float eps; // small tolerance number
  int testadj; // =1 test adjoint
  int smooth; // =1 smoothing filter
  int KEEP_WM; // If keepwm==1 the value for Wm from previous CMP is used

  ///////////////////////////////////////////////////////////////
  int j, i, k, ih, ix; // General counters 
  register int it;
  float *d;      /* single trace */
  float **cmp;   /* Common Scattering gather */
  int oldcdp;
  float  *m;    /* Final migrated trace   */
  float *x;      /* axis for CDPs */
  float *he;     /* axis for input offset */
  float *h2;     /* axis for output offset if resample==1 */
  float  *t;     // time axis for input and output
  float  h;      // halfoffset
  float t0;      // First useful time
  unsigned int ntrmax;  // Numebr of traces to process from data file
  unsigned int ncmp;  // Number of CMPs to process per reading loop

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
  float *velatcmp; /* array[nt] of vel for a particular trace */ 
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
  float  h2min;          //limit for the equivalent offset
  float  h2max;          //limit for the equivalent offset     
  int verbose;		/* flag for echoing info		*/
  int output;           /* =2 output are original traces, used for test */
			/* =4 output are  Radon                 */
                        /* =2 output are CMP gathers,           */ 
  float fmax;
  int *cdpt;         /* Vector with cdps number */
  
  // deviations from velocity law
  int nq;  
  float qmin;
  float qmax;
  float *q;
  float dq;
  float **mm;            // Temporal array for the model

  //  LS mig  
  float step;
  float eps1;
  float eps2;
  int itercg;
  int itervel;
  int iter_end;
  int buffer;
  int norm;
  float *Wm;   // Factorization of inv model covariance matrix 
  float *Wd;   // Factorization of inv residual covariance matrix 
  unsigned short mute;
  float parmute;
  int sinc; //=1 uses sinc interpolation
  int lsmethod;
  float factor;
  float smute;
  float nmofactor;
  int **foldx;
  int *nhcmp;
  float under;
  float over;

  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
  // Register starting time
  start=time(0);
   // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  traux=tr; // Keep traux for use later
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
  if (!getparfloat("h2min", &h2min))  h2min = 0;
  if (!getparfloat("h2max", &h2max))  h2max = 1000;
  if (!getparint("nhe", &nhe))  nhe = 50;
  if (!getparint("verbose", &verbose)) verbose = 0;  
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
  if (!getparint("norm",&norm)) norm=0;
  if (!getparint("buffer", &buffer))  buffer = 5;
  if (!getparuint("ntrmax",&ntrmax)) ntrmax = 1000000;
  if (!getparuint("ncmp",&ncmp)) ncmp = 1 ;
  if (!getparfloat("parmute",&parmute)) parmute=1;
  if (!getparushort("mute",&mute)) mute=0;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparint("keepwm",&KEEP_WM)) KEEP_WM=0;
  if (!getparint("sinc",&sinc)) sinc=0;
  if (!getparint("lsmethod",&lsmethod)) lsmethod=0;
  if (!getparfloat("under",&under)) under=0.999;
  if (!getparfloat("over",&over)) over=1.001;
  if (!getparfloat("factor",&factor)) factor=0.8;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=1.9;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
  if (!getparstring("cmpfile",&cmpfile)) cmpfile="cmpfile.su";

  if ((modelfilep=fopen(modelfile,"w"))==NULL)
    err("cannot open file=%s\n",modelfile);

  cmpfilep=efopen(cmpfile,"w");
  FILE *offsetfilep;
  cwp_String offsetfile=""; /* file containing positions */
  int nn;
  if (getparstring("offsetfile",&offsetfile)){
    h2=ealloc1float(nhe);
    fprintf(stderr,"New offset given by %s\n",offsetfile);  
    offsetfilep=efopen(offsetfile,"r");
    ih=0;
    do{
      nn=fscanf(offsetfilep,"%f",&h2[ih]); 
      ih++;
    }while(nn==1);
    nh2=ih-1;
    fprintf(stderr,"nh2= %d\n",nh2);
    efclose(offsetfilep);
  }
  else{
    warn("No offset file, using the original \n");
    h2=he;
  }
  if (nh2>nhe) err("new offset require more traces than available ==> Increment nhe. \n");
  nx=(int) floor((cdpmax-cdpmin)/dxcdp + 1);

  if (verbose) fprintf(stderr,"output=%d\n",output);
  if (verbose) fprintf(stderr,"number of ouput cdps = %d\n",nx);
  if (verbose && smooth) fprintf(stderr,"Smoothed Output \n");
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
  
  // Allocate memory for data and model
  
  if ((d=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");

  if ((cmp=ealloc2float(nt,nhe))==NULL)
    fprintf(stderr,"***Sorry, space for cmp could not be allocated\n");

  if ((m=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for m could not be allocated\n");
    
  if ((x=ealloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for x could not be allocated\n");

  if ((he=ealloc1float(nhe))==NULL)
    fprintf(stderr,"***Sorry, space for he could not be allocated\n"); 
  
  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n");

  if ((q=ealloc1float(nq+1))==NULL)
    fprintf(stderr,"***Sorry, space for q could not be allocated\n");

  if ((mm=alloc2float(nt,nq))==NULL)
    fprintf(stderr,"***Sorry, space for mm could not be allocated\n");

  if ((Wm=alloc1float(nt*nq))==NULL)
    fprintf(stderr,"***Sorry, space for Wm could not be allocated\n");

  if ((Wd=alloc1float(nt*nhe))==NULL)
    fprintf(stderr,"***Sorry, space for Wd could not be allocated\n");

  foldx=ealloc2int(nhe,nx);

  nhcmp=ealloc1int(nx);

  int *hflag=ealloc1int(nhe);
  // Velocity axis for Radon space

  if (nq>1){
    dq = (qmax-qmin)/(nq-1);
    for (i=0;i<nq;i++) q[i] = qmin+i*dq;
  }

  else  q[0]=0;

  // arrays for nmo model
  // Velint will contain velocity law for one cmp
  if ((velint=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Space for velint could not be allocated\n");
  // velatcmp will contain velocity law for ncmps 
  if ((velatcmp=ealloc1float(nt))==NULL)
    fprintf(stderr,"Cannot allocate velatcmp\n");

  /* Create axis for output cpds, equivalent offset and time */
  dx=dxcdp;
  for(ix=0;ix<nx;ix++) x[ix]=cdpmin+dx*ix;
  for(it=0;it<nt;it++) t[it]=0+it*dt;

  ///////////////////////////////////////////////////////////////////////

  fprintf(stderr,"Starting cdp loop ****\n");  

  /* Loop to compute ncmp CMP gathers from x[ix] to x[ix+ncmp] */
  for (ix=0;ix<nx;ix++){
    fprintf(stderr,"cdp=%f\n",x[ix]);   

    /* Initialize to zero the ouput trace and the cmp gather corresp to ix */  
    for (ih=0;ih<nhe;ih++) 
      for (it=0;it<nt;it++)
	cmp[ih][it]=0; // CMP trace
    foldx[ix][ih]=0;	

    erewind(tracefp);
    erewind(headerfp);   
    // Save the current cdp (x[ix]) into a temporal file
    ntr = 0;
    tr=traux; //Use the trace already read in last pass 
    do {
      ntr++;
      cdpt[ntr]=(int) tr.cdp;
      //      fprintf(stderr,"tr.cdp=%d\n",tr.cdp);
      efwrite(&tr,HDRBYTES,1,headerfp);
      efwrite(tr.data,FSIZE, nt, tracefp);
    } while ((tr.cdp <= x[ix]) && gettr(&tr));
    traux=tr;
    erewind(tracefp);
    erewind(headerfp);  
    fprintf(stderr,"ntr=%d\n",ntr);
    if (ntr<3) continue;
    int ntread=0;
    ih=0;
    efread(tr.data, FSIZE, nt, tracefp);
    efread(&tr,HDRBYTES,1,headerfp);
    ntread++;

    /* Loop for mapping one trace to the CMP gather */

    /* set old cdp and old offset for first trace */
    oldcdp = tr.cdp;
    ih=0;
    while(ntread<ntr){     
      //      fprintf(stderr,"tr.cdp=%d,tr.offset=%d\n",tr.cdp,tr.offset);
      // Add the read trace to every cdp
      for (it=0;it<nt;it++) cmp[ih][it]=tr.data[it]; // CMP trace
      foldx[ix][ih]=1;
      he[ih]=tr.offset;
      efread(tr.data, FSIZE, nt, tracefp);
      efread(&tr,HDRBYTES,1,headerfp);
      ntread++;
      ih++;
    };
    nhcmp[ix]=ih;
    /* compute new square slowness and anis function */
    interpovv(nt,ncdp,cdp,ovv,oa1,oa2,x[ix],velint,oa1t,oa2t);
    memcpy((void *) velatcmp,
	   (const void *) velint,nt*sizeof(float));
    // Check the cmp array
    
    if (0){
      erewind(headerfp);  
      for (itr=0;itr<nhe;itr++){
	if (foldx[ix][itr]){
	  efread(&tr,HDRBYTES,1,headerfp);
	  memcpy((void *) tr.data,
		 (const void *) cmp[itr],nt*sizeof(float));
	  fputtr(cmpfilep,&tr);
	}    
      }
    }
    
    radonfreqint(he,nhcmp[ix],cmp,t,nt,dt,velatcmp,eps1,eps2,eps,mm,q,nq,itercg,iter_end,step,lsmethod,factor,smute,nmofactor,mute,parmute,h2,nh2);
    
    memset( (void *) hflag, (int) '\0', nh2 * ISIZE);
    sethflag2(hflag,nhe,he,nhe,h2,nh2,under,over);

    if(0)
      for  (ih=0;ih<nh2;ih++)
	fprintf(stderr,"hflag[%d]=%d,he[%d]=%f,h2[%d]=%f\n",ih,hflag[ih],ih,
		he[ih],ih,h2[ih]);

    // Output trace
    erewind(tracefp);
    erewind(headerfp);  
    

    for (ih=0;ih<nh2;ih++){
      if ((hflag[ih]==1)&&(ih<nhe)){
	efread(&tr,HDRBYTES,1,headerfp);
	efread(tr.data,FSIZE, nt, tracefp);
	if (verbose) fprintf(stderr,"Original data ===> tr.offset=%d\n",tr.offset);  
      }
      else{
	tr.offset=(int) h2[ih];
	tr.f2=0;
	tr.cdp=(int) x[ix]; // front of the CMP
	tr.ntr=nh2;
	memcpy((void *) tr.data,
	       (const void *) cmp[ih],nt*sizeof(float));
	if (verbose) fprintf(stderr,"Interpolated  data ======> tr.offset=%d\n",tr.offset);  
      }
      puttr(&tr);
    }

    for (itr=0;itr<nq;itr++){
      for (it=0;it<nt;it++) tr.data[it]=mm[itr][it];
      tr.cdp=(int) x[ix]; // front of the CMP
      tr.dt=(int) (dt*1e6);       
      tr.ntr=nq;
      tr.ns=nt;
      tr.tracl=itr+1;
      tr.tracr=itr+1;
      tr.offset=(int) (q[itr]);
      tr.f2=q[itr];
      tr.sx=(int) x[ix];
      tr.gx=(int) x[ix];
      fputtr(modelfilep,&tr);    
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
  free2float(mm);
  free1float(q);
  free1float(velatcmp);
  free1float(velint);
  free1float(t);
  free1float(he);
  free1float(x);
  free1float(m);
  free2int(foldx);
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






























