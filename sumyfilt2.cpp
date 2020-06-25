/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUKMIG:  $Date: March 1999  */
#include <iostream.h>
#include "su.h"
#include "segy.h"
#include "clibrarytd.h"
#include "inversion_par.hpp"
/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUMYFIT - Filter data with different filters.                       ", 
  "	      Useful mainly for testing filters.                  	",
  " 	   								",
  " sukmig0 < stdin > stdout [optional parameters]          		",
  " 									",
  " filt = 0    No filter                                               ",
  "                          						",
  "                                                   			",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt
/**************** end self doc ***********************************/

segy tr,trf; 

void rjwfilter(float *d,int nt,float dt);
void filt(float *trace,int nt,float dt,float fmax,int ls,int m,float *trf);
void rho_filter(int npoints, int nt, float dt, float *rho);

int main(int argc, char **argv)
{
  inversion inv;
  inv.SetItercg(20);
  fprintf(stderr,"Itercg=%d\n",inv.GetItercg());
  //  ~inversion inv;

  //  inv.itercg=20;	/* Number maximum of iterations in CG*/
  inv.iter_end=1;	/* Number maximum of External Iterations */
  inv.eps1=1e-3;	/* noise covariance    */
  inv.eps2=1e-3;	/* model covariance */
  inv.eps=1e-7;	/* small number for tolerance */
  inv.step=0.95;	/* step length factor  */
  inv.norm=0;     /* norm to use in model weights; */
  inv.restart=1;  /* always set to 1 for now */
  inv.taperflag=0; /* If set applies taper to 5 outer traces */
  inv.mute=0;      /* if set applies mute */
  //  cout << " Defining inversion\n";

  int j,i,k,ix;
  register int it;
  float *d;
  float *t; 
  int nt;
  float dt;	
  const double  pi=acos(-1.);
  float fmax;
  int filter;
  int lfilter; /* length (number of points ) of the filter */
  float *rho;
  // Initialize 

  initargs(argc, argv);
  requestdoc(1);
   
  // Get parameters 
  if (!getparint("filter", &filter))  filter = 0;
  if (!getparint("lfilter", &lfilter))  lfilter = 50;  
  // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");

  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;

  fmax=0.8/(2*dt);
  fprintf(stderr,"nt=%d,dt=%f,filter=%d,fmax=%f\n",nt,dt,filter,fmax); 
  
  // Allocate memory for data and model
  
  if ((d=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for d could not be allocated\n");
        
  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n"); 

 
  if ((rho=ealloc1float(lfilter))==NULL)
    fprintf(stderr,"***Sorry, space for rho could not be allocated\n"); 
  
  if (filter==3) rho_filter(lfilter,nt,dt,rho);
  
  int np2=lfilter/2;	        
  for(it=0;it<nt;it++) t[it]=it*dt;
  fprintf(stderr,"**************************\n");  
  
  j=0;
  do {   // Loop over traces
    j++; 
    /* SU phase shift filter and lowpass filter */
    if (filter==1) filt(tr.data,nt,dt,fmax,0,lfilter,trf.data);
    /* My phase shift filter */
    if (filter==2) rjwfilter(tr.data,nt,dt);
    /* Rho filter */
    if (filter==3){
      trf=tr; // Warning: Check this
      conv(nt,-np2,trf.data,lfilter,0,rho,nt,0,tr.data);
    }
     
    puttr(&tr);
  }while (gettr(&tr));

  free1float(t);
  free1float(d);
  free1float(rho);

  return EXIT_SUCCESS;
}




















