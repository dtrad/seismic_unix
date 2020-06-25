/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUEOMIG:  $Date: March 1999  */
#include <time.h>
#include "kmig.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUkMIGLS                   - Kirchoff  migration                     ",
  "                                                                     ", 
  "	                                                          	",
  " 	   								",
  " sukmigls < stdin > stdout [optional parameters]          		",
  " 									",
  " cdpmin=0        Fisrt CDP in meters                                 ",
  " cdpmax=100      Last CDP 						",
  " vel=2000        Velocity of migration (in meters) 			",
  "                                                                     ",
  " output=1        =0 test - Output=input                              ",
  "                 =1 Migrated zero offset section                     ",
  "                 =2 CSP gather                                       ",
  "                 =3 Least square migrated traces			",
  "                 =4 Least square CSP gather                          ",
  " aper=10         Number of CMPs to include at every CDP              ",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset
/**************** end self doc ***********************************/

/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/

//////////////////////////////////////////////////////////////////////
void kmig(float *d, float cdp, float h,float **m, float *t, float *x, float **ivel2, int nt, int nx, float dt, float cdpspace, float aper, int adj); 

int main(int argc, char **argv)
{ 
  /* time counting */
  time_t start,finish; 
  double elapsed_time;

  segy tr,trf; 
  int nh;   // number of offset traces
  int nt;   // number of time samples 
  int nx;    // number of midpoints
  float dt; // Time sampling
  float dh; // offset interval
  float dx;  // midpoint interval
  float eps; // small tolerance number
  int testadj; // =1 test adjoint
  int smooth; // =1 smoothing filter
  float aper;
  ///////////////////////////////////////////////////////////////
  int j, i, k, ih, ix, itr; // General counters 
  register int it;
  float *d;      /* single trace */
  float **m;    /* Final migrated output m[ix][it]   */
  float *x;     /* axis for CDPs */
  float *t;     // time axis for input and output
  float  h;     // halfoffset
 
  /* Velocity */
  // Variables used to interpolate velocities (copied from sunmo.c) 
  int ncdp;	/* number of cdps specified */
  float *cdp;	/* array[ncdp] of cdps */
  float **ovv;
  float **velint;
  float **ivel2;
  int icdp;
  /////////////////////////////////////////////////////


  const double  pi=acos(-1.); 
  float xm;    // midpoint distance from trace to scp = cdp_tr - x_scp
  int ntr;     // Total number of input traces
  float  cdpmin;          //limit for the output
  float  cdpmax;          //limit for the output
  float  dxcdp;           // output cdp interval 
  int verbose;		/* flag for echoing info		*/
  char *tmpdir;		/* directory path for tmp files		*/
  int output;           /* =0 output are original tracces       */ 
			/* =1 output are migrated traces        */
  float *fold;          /* Fold for final imaged model */
  float dfold;       // the fold is increase by dfold for every trace in cdp
  float fmax;
  int *cdpt;         /* Vector with cdps number */
  float beta;        /* max angle (degrees) to include in migration */
  int ntrmax;
  //  LS mig
  // Not used in this algorithm for now  
  float step;
  float eps1;
  float eps2;
  int itercg;
  int iter_end;
  int buffer;
  //////////////////////////
  float cdpspace;
  float cdptr;
  int adj;

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
  if (!(ntr=(int) tr.ntr)) err("***ntr must be set\n");
  nh=ntr;
  fmax=0.8/(2*dt);
  
  fprintf(stderr,"nt=%d,dt=%f\n",nt,dt);  

  // Get parameters 
  
  if (!getparfloat("cdpmin", &cdpmin))  cdpmin = 0;
  if (!getparfloat("cdpmax", &cdpmax))  cdpmax = 1000;
  if (!getparfloat("dxcdp", &dxcdp))  dxcdp = 10;
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparint("output", &output)) output = 1;
  if (!getparfloat("eps", &eps))  eps = 1e-7;
  if (!getparfloat("eps1", &eps1))  eps1 = 1e-3;
  if (!getparfloat("eps2", &eps2))  eps2 = 1e-3;
  if (!getparfloat("step", &step))  step = 0.9;
  if (!getparint("itercg", &itercg))  itercg = 10;
  if (!getparint("iter_end", &iter_end))  iter_end = 3;
  if (!getparint("testadj", &testadj))  testadj = 0;
  if (!getparint("smooth", &smooth))  smooth = 0;
  if (!getparfloat("dfold", &dfold))  dfold = 0.01;
  if (!getparfloat("cdpspace", &cdpspace))  cdpspace = 1; 
  if (!getparint("ntrmax", &ntrmax))  ntrmax = 900000;
  if (!getparint("adj", &adj))  adj = 1;
  if (!getparfloat("aper",&aper)) aper=10000;

  
  nx=(int) floor((cdpmax-cdpmin)/dxcdp + 1 + 0.5 );
  
  if (verbose) fprintf(stderr,"***Kirchhoff Migration***\n");
  if (verbose) fprintf(stderr,"output=%d\n",output);
  if (verbose) fprintf(stderr,"number of ouput cdps = %d\n",nx);
  if (verbose && smooth) fprintf(stderr,"Smoothed Output \n");
  //if (verbose)  fprintf(stderr,"beta=%f\n",beta);
  if (verbose)  fprintf(stderr,"dfold=%f\n",dfold);

  beta=beta*pi/180.;
  /// This part is copied from sunmo.c ///////////////////////
  // It interpolates the original velocities- time table in file
  // obtained by velocity analisys. The interpolated velocity model 
  // has the dimensions of the final output, i.e., velint[ix][it]
  // When Kirchhoff summation is performed it is used to computed the 
  // time in the data space that corresponds to the output
  // Anisotropy is included for future work but not used at all
  // in the program.
  //////////////////////////////////////////////////////////// 
  /* get velocity functions, linearly interpolated in time */
  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  getvelocities(dt,nt,ncdp,cdp,ovv);
  for (ix=0;ix<ncdp;ix++) 
    fprintf(stderr,"cdp[%d]=%f\n",ix,cdp[ix]); 
  ///////////////////////////////////////////////////////////////////////////
  // cdpt is used to keep the cdp number for every trace
  if ((cdpt=ealloc1int(ntr+1))==NULL)
    fprintf(stderr,"***Sorry, space for cdpt could not be allocated\n");

  
  // Allocate memory for data and model
  // The model can be referenced trace by trace, i.e., m[ix] 
  // or point by point m[ix][it]
  // The same for velint

  d=ealloc1float(nt);
  m=ealloc2float(nt,nx);
  x=ealloc1float(nx);
  fold=ealloc1float(nx);
  t=ealloc1float(nt);
   // array for interpolated velocity
  velint=ealloc2float(nt,nx);
  ivel2=ealloc2float(nt,nx);
  /* Create axis for output cpds, equivalent offset and time */
  
  dx=dxcdp;      
  for(ix=0;ix<nx;ix++) x[ix]=cdpmin+dx*ix;
  for(it=0;it<nt;it++) t[it]=0+it*dt;

  for (ix=0;ix<nx;ix++) fold[ix]=1;

  for (ix=0;ix<nx;ix++) for (it=0;it<nt;it++) m[ix][it]=0;  


  //  else for (ix=0;ix<nx;ix++) for (it=0;it<nt;it++) m[ix][it]=0;
  // fold is increased for every trace that belongs to the cdp but in a amount
  // proportional to the amplitude or energy
  // More precision would be better but for now you can set dfold until the 
  // output looks OK. A first value could be the amplitude (example as seen in   // clip with suxwigb

  for (itr=0;itr<ntr;itr++){
    ix=(int) floor((cdpt[itr]-cdpmin)/dx+.5);
    if ((ix<nx)&&(ix>=0)) fold[ix]+=dfold;
  }
  
  if (verbose){ 
    for (ix=0;ix<nx;ix++)
      fprintf(stderr,"fold[%d]=%f\n",ix,fold[ix]);
    for (icdp=0;icdp<ncdp;icdp++) 
      fprintf(stderr,"ovv[%d][100]=%f\n",icdp,ovv[icdp][100]);   
    fprintf(stderr,"**************************\n");
  } 
 
  /* Loop to compute velocities at every output by interpolation  */
  for (ix=0;ix<nx;ix++){
    interpovv(nt,ncdp,cdp,ovv,x[ix],velint[ix]);
    for (it=0;it<nt;it++) ivel2[ix][it]=1./(velint[ix][it]*velint[ix][it]);
  }
  
  itr = 0;
  do {
    cdpt[itr]=(int) tr.cdp;
    //    cdptr=tr.sx/1000.+tr.offset/2.;
    cdptr=(float) tr.cdp;    
    if (verbose) 
      fprintf(stderr,"tr.sx=%d,tr.tracl=%d,tr.tracf=%d\n",tr.sx,
	      tr.tracl,tr.tracf);
    h=(float) tr.offset;
    h=h/2.;
    memcpy(d, tr.data, nt*sizeof(float));
    
    //kmig3(d,cdptr,h,m,t,x,velint,nt,nx,dt);
    kmig(d,cdptr,h,m,t,x,ivel2,nt,nx,dt,cdpspace,aper,adj);

    if ((0)&&(fmod(tr.tracl,100.)==0)) fprintf(stderr,"tr.tracl=%d\n",tr.tracl);
    itr++;
  } while (gettr(&tr) && (itr < ntrmax));
  ntr=itr;
  fprintf(stderr,"ntr=%d\n",ntr);
  


  if (output==1){
    for (ix=0;ix<nx;ix++){
      /*
	for (it=0;it<nt;it++)
	m[ix][it]/=(fold[ix]*velint[ix][0]*sqrt(pi*MAX(t[it],1e-2)));
      */
      memcpy(tr.data, m[ix], nt*sizeof(float));
      tr.cdp=(int) x[ix];
      tr.dt=(int) (dt*1e6);
      tr.ntr=nx;
      tr.ns=nt;
      tr.tracl=ix+1;
      tr.tracr=ix+1;
      tr.offset=0;
      tr.sx=(int) x[ix];
      tr.gx=(int) x[ix];
      //filt(tr.data,nt,dt,fmax,0,50,trf.data);
      puttr(&tr);
      fprintf(stderr,"tr.cdp=%d\n",tr.cdp);    
    } 
  }

    
  free1int(cdpt);
  free1float(cdp);
  free2float(ovv);
  free2float(ivel2);
  free2float(velint);
  free1float(t);
  free1float(fold);
  free1float(x);
  free2float(m);
  free1float(d);

  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);
  
  return EXIT_SUCCESS;
}




















