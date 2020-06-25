/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */

/* SURADTD:  $Date: June 1999  */
#include "su.h"
#include "clibrarytd.h"

void smoothing(float *d,int nt,int nx,int nl,int nr,int flag)
{
  float *dtemp1;
  float *dtemp2;
  float *filter;
  int ix,it; 
    
  dtemp1=ealloc1float(nt);
  dtemp2=ealloc1float(nt);
  filter=alloc1float(nl+nr+1);


  fprintf(stderr,"smoothing......\n");  
  rwa_smoothing_filter (flag,nl,nr,filter);
  for (ix=0;ix<nx;ix++){
    for (it=0;it<nt;it++) dtemp1[it]=d[it+ix*nt];
    convolve_cwp (nl+nr+1,-nl,filter,nt,0,dtemp1,nt,0,dtemp2);
    for (it=0;it<nt;it++) d[it+ix*nt]=dtemp2[it];    
  }
    
  free1float(filter);
  free1float(dtemp1);
  free1float(dtemp2);
  
  return;
} 


void smoothing(float **d,int nt,int nx,int nl,int nr,int flag)
{
  float *dtemp1;
  float *dtemp2;
  float *filter;
  int ix,it; 
    
  dtemp1=ealloc1float(nt);
  dtemp2=ealloc1float(nt);
  filter=alloc1float(nl+nr+1);


  fprintf(stderr,"smoothing......\n");  
  rwa_smoothing_filter (flag,nl,nr,filter);
  for (ix=0;ix<nx;ix++){
    for (it=0;it<nt;it++) dtemp1[it] = d[ix][it];
    convolve_cwp (nl+nr+1,-nl,filter,nt,0,dtemp1,nt,0,dtemp2);
    for (it=0;it<nt;it++) d[ix][it] = dtemp2[it];    
  }
  int nlh = 5;
  int nrh = 5;
  rwa_smoothing_filter (flag,nlh,nrh,filter);
  for (it=0;it<nt;it++){
       for (ix=0;ix<nx;ix++) dtemp1[ix] = d[ix][it];
       convolve_cwp (nlh+nrh+1,-nlh,filter,nx,0,dtemp1,nx,0,dtemp2);
       for (ix=0;ix<nx;ix++) d[ix][it] = dtemp2[ix];    
  }
       
  free1float(filter);
  free1float(dtemp1);
  free1float(dtemp2);
  
  return;
} 



void smoothing(float *d,int nt, int nl,int nr,int flag)
{
  float *dtemp1;
  float *dtemp2;
  float *filter;
  int ix,it; 
    
  dtemp1=ealloc1float(nt);
  dtemp2=ealloc1float(nt);
  filter=alloc1float(nl+nr+1);


  fprintf(stderr,"smoothing......\n");  
  rwa_smoothing_filter (flag,nl,nr,filter);
  

  // first filter out outliers
  for (it=0;it<nt;it++) dtemp1[it]=d[it];
  convolve_cwp (nl+nr+1,-nl,filter,nt,0,dtemp1,nt,0,dtemp2);
  // take care of borders
  for (it=0;it<nl;it++) dtemp2[it]=d[it];
  for (it=nt-nr;it<nt;it++) dtemp2[it]=d[it];
  // filter outliers in input;
  float sum=0;
  for (it=nl;it<nt-nr;it++) sum+=(fabs(d[it]-dtemp2[it]));
  sum/=(nt-nr-nl);
  for (it=nl;it<nt-nr;it++){
    if (fabs(d[it]-dtemp2[it])<=fabs(sum)) dtemp1[it]=d[it];
    else dtemp1[it]=dtemp2[it];
  }
  // filter again after outliers are out
  convolve_cwp (nl+nr+1,-nl,filter,nt,0,dtemp1,nt,0,dtemp2);
  //copy to output preserving original on borders.
  for (it=nl;it<nt-nr;it++) d[it]=dtemp2[it];    
  
    
  free1float(filter);
  free1float(dtemp1);
  free1float(dtemp2);
  
  return;
} 






