#include "su.h"
#include "clibrary.h"
#include "radon.h"
#include "radonl1freq.h"
/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void radonl1freq(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, float *vel, int itercg, int iter_end, op_param op2, float smute, float nmofactor, float *ffilter, float *amps)
{
  int it, ih;
  float *dtemp=0;
  int nq=op2.nq1+op2.nq2;

  dtemp=alloc1float(nt);
   
  //////////////////////////////////////////////////////////////////////
  
  for (ih=0;ih<nh;ih++){
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
  }

  plotgather(data,nh,nt,dt,"suxwigb perc=90 title=\"plotgather\"");
  
  memset( (void *) model[0], (int) '\0', nq * nt *FSIZE);  

  radonl1freq_loop(data,model,t,q,h,nt,nq,nh,dt,iter_end,itercg,op2);

  /////////////////////////////

  for (ih=0;ih<nh;ih++){
    for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
  }    

  free1float(dtemp);

  return;

}
























