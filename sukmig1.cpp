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
  " vel=2000        Velocity of migration (in meters) 			",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset
/**************** end self doc ***********************************/

segy tr,trf; // for malloc use *trr
void kmig1(float *d, float cdp, float h, float **m, float *t, float *x, float vel) ;
void rjwfilter(float **d,int nt,int nh, float dt);
void filt(float *trace,int nt,float dt,float fmax,int ls,int m,float *trf);

int nt,nh,nx; // Size of volume x,h,t
float dt,dh,dx; // Sampling
int main(int argc, char **argv)
{
  int j,i, k, ix;
  register int it;
  float *d, **m;
  float *x, *t, h, cdp, vel, cdpmin, cdpmax; 
  extern int nt, nh, nx;
  extern float dt,dh,dx;	
  const double  pi=acos(-1.);
  float fmax;  
  // Initialize 
  initargs(argc, argv);
  requestdoc(1);
   
  // Get parameters 
  if (!getparfloat("vel", &vel))  vel = 2000;
  if (!getparfloat("cdpmin", &cdpmin))  cdpmin = 0;
  if (!getparfloat("cdpmax", &cdpmax))  cdpmax = 1000;
  if (!getparint("nx", &nx))  nx = 100;
  
  // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");

  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  fmax=0.8/(2*dt);
  fprintf(stderr,"nt=%d,dt=%f,nx=%d,vel=%f,fmax=%f\n",nt,dt,nx,vel,fmax); 
  
  // Allocate memory for data and model
  
  if ((d=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");
  
  if ((m=ealloc2float(nx,nt))==NULL)
    fprintf(stderr,"***Sorry, space for m could not be allocated\n");
    
  if ((x=ealloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for x could not be allocated\n");
    
  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n"); 
  
  dx=(cdpmax-cdpmin)/(nx-1);      
  for(ix=0;ix<nx;ix++) x[ix]=cdpmin+dx*ix;
  for(it=0;it<nt;it++) t[it]=it*dt;
  fprintf(stderr,"**************************\n");  
  
  for (ix=0;ix<nx;ix++) for (it=0;it<nt;it++) m[it][ix]=0; // Output trace
  j=0;
  do {   // Loop over traces
    j++; 
    //fprintf(stderr,"j=%d,tr.cdp=%d\n",j,tr.cdp);    
    h=(float) tr.offset;
    h/=2;  // halfoffset
    cdp=(float) tr.cdp;
    // SU phase shift filter and lowpass filter
    filt(tr.data,nt,dt,fmax,0,50,trf.data);
    for (it=0;it<nt;it++)
      d[it]=(float) tr.data[it];    
    kmig1(d,cdp,h,m,t,x,vel);
  }while (gettr(&tr));
  //for (ix=0;ix<nx;ix++) 
  //  for (it=0;it<nt;it++)
  //    m[it][ix]/=(vel*sqrt(pi*MAX(t[it],1e-3)));
  //rjwfilter(m,nt,nx,dt); // My phase shift filter
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
    puttr(&tr);
  }   
  free1float(t);
  free1float(x);
  free2float(m);
  free1float(d);
  
  return EXIT_SUCCESS;
}




















