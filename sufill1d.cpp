/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SUFILL1D:  $Date: November 1999  */

#include "su.h"
#include "segy.h"
#include "clibrarytd.h"

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUFILL1D - Fill data in one dimension using Claerbout's CG.         ", 
  "	      Useful mainly for testing filters.                  	",
  " 	   								",
  " sufill1d < stdin > stdout [optional parameters]          		",
  " 									",
  " filt = 0   interpolate time using given filter                      ",
  " filt = 1   interpolate time using training data 			",
  "            to find the filter                     			",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, dt
/**************** end self doc ***********************************/
void miss1(int na,int np,float *data, int niter, float *known);
segy tr,trf; 
void  missfip(int ntt, int na,int np,float *data,float *known,int niter);
void copy(int n, float *xx, float *yy);
void shaper(int nf,int nx,float *data, int niter);
void iner(int nf,int nr,float *data, int niter, int lag, int gap1, int gapn);
  
int main(int argc, char **argv)
{
  int j,i,k,ix;
  register int it;
  float *known;
  float *t; 
  int nt;
  float dt;	
  const double  pi=acos(-1.);
  float fmax;
  int filter;
  int niter;
  int lfilter; /* length (number of points ) of the filter */
  int ntraining; // length of training data
  int lag;
  int gap1;
  int gapn;

  // Initialize 

  initargs(argc, argv);
  requestdoc(1);
   
  // Get parameters 
  if (!getparint("filter", &filter))  filter = 0;
  if (!getparint("ntraining", &ntraining))  ntraining = 10; 
  if (!getparint("lfilter", &lfilter))  lfilter = 3; 
  if (!getparint("niter",&niter)) niter = 5; 
  if (!getparint("lag",&lag)) lag = 0; 
  if (!getparint("gap1",&gap1)) gap1 = 0; 
  if (!getparint("gapn",&gapn)) gapn = 1; 

  // Get info from first trace 
  
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");

  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;

  fmax=0.8/(2*dt);
  fprintf(stderr,"nt=%d,dt=%f,filter=%d,fmax=%f\n",nt,dt,filter,fmax); 
  
  // Allocate memory for data and model
  
  if ((known=alloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for known could not be allocated\n");
        
  if ((t=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for t could not be allocated\n"); 

  for(it=0;it<nt;it++) t[it]=it*dt;
  fprintf(stderr,"**************************\n");  

  for(it=0;it<nt;it++) known[it]=1;
  //known[10]=0; 
  if ((filter==1)||(filter==2)) 
    for(it=0;it<nt-10;it+=4) known[it]=0;
  
  j=0;
  do {   // Loop over traces
    j++; 
    if (filter==0) for(it=0;it<nt-10;it+=2) tr.data[it]=0; 
    else if (filter==1) miss1(lfilter,nt,tr.data,niter,known); 
    else if (filter==2) missfip(ntraining,lfilter,nt,tr.data,known,niter);
    else if (filter==3) shaper(lfilter,nt,tr.data,niter);
    else if (filter==4) iner(lfilter,nt,tr.data,niter,lag,gap1,gapn);
    puttr(&tr);
  }while (gettr(&tr)); 
  free1float(t);
  free1float(known);

  return EXIT_SUCCESS;
}




















