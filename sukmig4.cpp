/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUEOMIG:  $Date: March 1999  */

#include "su.h"
#include "segy.h"
#include "clibrarytd.h"
#include "header.h"
#include <signal.h>

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUkMIG3                   - Kirchoff  migration                     ",
  "                                                                     ", 
  "	                                                          	",
  " 	   								",
  " sueomig0 < stdin > stdout [optional parameters]          		",
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

///////////////////////////////////////////////////////////////////
// Global Variables 


//////////////////////////////////////////////////////////////////
// Function Prototypes

void rjwfilter(float **d,int nt,int nh, float dt);
void filt(float *trace,int nt,float dt,float fmax,int ls,int m,float *trf);
void interpovv(int nt, int ncdp, float *cdp, float **ovv, float **oa1, 
	       float **oa2, float cdpt, float *ovvt, float *oa1t, float *oa2t);
void kmig3(float *d, float cdp, float h,float **m, float *t, float *x, float **vel, int nt, int nx, float dt, float cdpspace, float aper);
void kmig3(float *d, float cdp, float h,float **m, float *t, float *x, float **vel, int nt, int nx, float dt); 
//////////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{
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
  int icdp;	/* index into cdp array */
  int jcdp;	/* index into cdp array */
  int nvnmo;	/* number of vnmos specified */
  float *vnmo;	/* array[nvnmo] of vnmos */
  int ntnmo;	/* number of tnmos specified */
  float *tnmo;	/* array[ntnmo] of tnmos */
  float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  float **velint;/* array[nx][nt] of vel for a particular trace */
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
 
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

   // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");

  //if (!tr.offset) err("offset header field must be set");
  
  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  if (!(ntr=(int) tr.ntr)) err("***ntr must be set\n");
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
  if (!getparfloat("aper",&aper)) aper=10;

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
      if (tnmo[it]<=tnmo[it-1])
	err("tnmo values must increase monotonically");

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

  ///////////////////////////////////////////////////////////////////////////
  // cdpt is used to keep the cdp number for every trace
  if ((cdpt=ealloc1int(ntr+1))==NULL)
    fprintf(stderr,"***Sorry, space for cdpt could not be allocated\n");

  
  // Allocate memory for data and model
  // The model can be referenced trace by trace, i.e., m[ix] 
  // or point by point m[ix][it]
  // The same for velint
  
  if ((d=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");

  if ((m=ealloc2float(nt,nx))==NULL)
    fprintf(stderr,"***Sorry, space for m could not be allocated\n");
  
  if ((x=ealloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for x could not be allocated\n");
   
  if ((fold=ealloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for fold could not be allocated\n");   

  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n");

   // array for interpolated velocity
  
  if ((velint=ealloc2float(nt,nx))==NULL)
    fprintf(stderr,"***Space for velint could not be allocated\n");

  /* Create axis for output cpds, equivalent offset and time */
  
  dx=dxcdp;      
  for(ix=0;ix<nx;ix++) x[ix]=cdpmin+dx*ix;
  for(it=0;it<nt;it++) t[it]=0+it*dt;

  for (ix=0;ix<nx;ix++) fold[ix]=1;
  for (ix=0;ix<nx;ix++) for (it=0;it<nt;it++) m[ix][it]=0;  

  // fold is increased for every trace that belongs to the cdp but in a amount
  // proportional to the amplitude or energy
  // More precision would be better but for now you can set dfold until the 
  // output looks OK. A first value could be the amplitude (example as seen in   // clip with suxwigb

  for (itr=0;itr<ntr;itr++){
    ix=(int) floor((cdpt[itr]-cdpmin)/dx+.5);
    if ((ix<nx)&&(ix>=0)) fold[ix]+=dfold;
  }

  if (verbose){ 
    for (ix=0;ix<nx;ix++) fprintf(stderr,"fold[%d]=%f\n",ix,fold[ix]);
    for (icdp=0;icdp<ncdp;icdp++) fprintf(stderr,"ovv[%d][100]=%f\n",icdp,ovv[icdp][100]);   
    fprintf(stderr,"**************************\n");
  } 
 
  /* Loop to compute velocities at every output by interpolation  */
  for (ix=0;ix<nx;ix++)
    interpovv(nt,ncdp,cdp,ovv,oa1,oa2,x[ix],velint[ix],oa1t,oa2t);
      
  ntr = 0;
  do {
    ntr++;
    cdpt[ntr]=(int) tr.cdp;
    //    cdptr=tr.sx/1000.+tr.offset/2.;
    cdptr=(float) tr.cdp;    
    //fprintf(stderr,"tr.sx=%d,tr.tracl=%d,tr.tracf=%d\n",tr.sx,tr.tracl,tr.tracf);
    h=(float) tr.offset;
    h=h/2.;
    memcpy(d, tr.data, nt*sizeof(float));
    //kmig3(d,cdptr,h,m,t,x,velint,nt,nx,dt);
    kmig3(d,cdptr,h,m,t,x,velint,nt,nx,dt,cdpspace,aper);
  } while (gettr(&tr) && (ntr < ntrmax));

  fprintf(stderr,"ntr=%d\n",ntr);
        
  if (output==1){
    for (ix=0;ix<nx;ix++){
      for (it=0;it<nt;it++)
	m[ix][it]/=(fold[ix]*velint[ix][0]*sqrt(pi*MAX(t[it],1e-2)));
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
      filt(tr.data,nt,dt,fmax,0,50,trf.data);
      puttr(&tr);
      fprintf(stderr,"tr.cdp=%d\n",tr.cdp);    
    } 
  }

  free1int(cdpt);
  free1float(cdp);
  free2float(ovv);
  free2float(oa1);
  free2float(oa2);
  free1float(oa1t);
  free1float(oa2t);
  free2float(velint);
  free1float(t);
  free1float(fold);
  free1float(x);
  free2float(m);
  free1float(d);

  return EXIT_SUCCESS;
}




















