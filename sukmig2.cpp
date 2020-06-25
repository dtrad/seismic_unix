/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUKMIG:  $Date: March 1999  */

#include "su.h"
#include "segy.h"
#include "clibrarytd.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUKMIG0 - First attempt for Kirchhoff migration                     ", 
  "	                                                          	",
  " 	   								",
  " sukmig0 < stdin > stdout [optional parameters]          		",
  " 									",
  " cdpmin=0        Fisrt CDP in meters                                 ",
  " cdpmax=100      Last CDP 						",
  " nx=0                                                                ",
  " vel=2000        Array with velocities of migration (in meters) 	",
  "                 at times time[i]                     		",
  " time=0          Array with times at which vel[i] is defined         ",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset
/**************** end self doc ***********************************/

segy tr,trf; // for malloc use *trr
void kmig2(float *d, float cdp, float h, float **m, float *t, float *x, float *vel) ;
void rjwfilter(float **d,int nt,int nh, float dt);
void filt(float *trace,int nt,float dt,float fmax,int ls,int m,float *trf);

int nt,nh,nx; // Size of volume x,h,t
float dt,dh,dx; // Sampling
int main(int argc, char **argv)
{
  int j,i, k, ix;
  register int it;
  float *d, **m;
  float *x, *t, h, cdp, cdpmin, cdpmax; 
  extern int nt, nh, nx;
  extern float dt,dh,dx;	
  const double  pi=acos(-1.);
  float fmax;

  /* Velocity */
  float *vel;    // Array for velocity function
  float *time;   // Array for time in velocity function 
  int   nvel;    // Number of velocities
  int   ntime;   // Number of time points for Velocities
  float *velint; // Array for velocity function after interpolation
  
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
   
  // Get parameters 
  if (!getparfloat("cdpmin", &cdpmin))  cdpmin = 0;
  if (!getparfloat("cdpmax", &cdpmax))  cdpmax = 1000;
  if (!getparint("nx", &nx))  nx = 100;

  nvel = countnparval(1,"vel");
  ntime = countnparval(1,"time");
  intprint(nvel);
  intprint(ntime);
       
  if (nvel!=ntime)
    err("number of vel and time values must be equal");  
  
  // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");

  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  fmax=0.8/(2*dt);
  fprintf(stderr,"nt=%d,dt=%f,nx=%d,fmax=%f\n",nt,dt,nx,fmax); 
  
  // Allocate memory for data and model
  
  if ((d=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");
  
  if ((m=ealloc2float(nx,nt))==NULL)
    fprintf(stderr,"***Sorry, space for m could not be allocated\n");
    
  if ((x=ealloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for x could not be allocated\n");
    
  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n"); 

  // arrays for nmo model
  
  if ((time=ealloc1float(ntime))==NULL)
    fprintf(stderr,"***Space for time could not be allocated\n");
  
  if ((vel=ealloc1float(nvel))==NULL)
    fprintf(stderr,"***Space for vel could not be allocated\n");

  if ((velint=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Space for velint could not be allocated\n");

  if (!getnparfloat(1,"vel",vel)) vel[0] = 2000.0;
  if (!getnparfloat(1,"time",time)) time[0] = 0.0;        
  for (it=0; it<ntime; ++it){
    fprintf(stderr,"vel[%d]=%e\n",it,vel[it]);
    fprintf(stderr,"time[%d]=%e\n",it,time[it]);
  }
  for (it=1; it<ntime; ++it)
    if (time[it]<=time[it-1])
      err("tnmo values must increase monotonically");

  /* Create axis for output cpds, equivalent offset and time */
  
  dx=(cdpmax-cdpmin)/(nx-1);      
  for(ix=0;ix<nx;ix++) x[ix]=cdpmin+dx*ix;
  for(it=0;it<nt;it++) t[it]=it*dt;

  /* Create axis for velocities */
  intlin(ntime,time,vel,vel[0],vel[nvel-1],nt,t,velint);

  fprintf(stderr,"**************************\n");  
  
  for (ix=0;ix<nx;ix++) for (it=0;it<nt;it++) m[it][ix]=0; // Output trace
  j=0;
  do {   // Loop over traces
    j++; 
    //fprintf(stderr,"j=%d,tr.cdp=%d\n",j,tr.cdp);    
    h=(float) tr.offset;
    h/=2;  // halfoffset
    cdp=(float) tr.cdp;
    for (it=0;it<nt;it++)
      d[it]=(float) tr.data[it];    
    kmig2(d,cdp,h,m,t,x,velint);
  }while (gettr(&tr));
  for (ix=0;ix<nx;ix++) 
    for (it=0;it<nt;it++)
      m[it][ix]/=(vel[0]*sqrt(pi*MAX(t[it],1e-2)));

  for (ix=0;ix<nx;ix++){  
    for (it=0;it<nt;it++)
      tr.data[it]=m[it][ix];
    tr.cdp=(int) x[ix];
    tr.dt=(int) (dt*1e6);
    tr.ntr=nx;
    tr.ns=nt;
    tr.tracl=ix;
    tr.tracr=ix;
    tr.offset=0;
    tr.sx=(int) x[ix];
    tr.gx=(int) x[ix];
    // SU phase shift filter and lowpass filter
    filt(tr.data,nt,dt,fmax,0,50,trf.data);
    puttr(&tr);
  }

  free1float(velint);
  free1float(time);
  free1float(vel);   
  free1float(t);
  free1float(x);
  free2float(m);
  free1float(d);
  
  return EXIT_SUCCESS;
}




















