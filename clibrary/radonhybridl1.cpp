#include "su.h"
#include "clibrary.h"
#include "radon.h"

/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void radonhybrid(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, int itercg, int iter_end, float step, float eps2, float eps1, float quantil, int norm, float qmin1, float qmin2,  int nq1, int nq2, float factor1, float factor2, float smute, float nmofactor, int rtmethod1, int rtmethod2, float depth1, float depth2, float fmax1, float fmax2, float *ffilter, float *amps)
{
  int it, ih;
  float *dtemp;
  float fmax;
  float freq=1/(2*dt)*0.5;
  float *Wd;
  float eps=1e-7;
  int testadj=0;
  int dq=0;

  dtemp=alloc1float(nt);
  Wd=alloc1float(nh);
   
  if (freq==0) fmax=1/(2*dt)*0.5;
  else fmax=freq;

  for (ih=0;ih<nh;ih++) Wd[ih]=1;

  //////////////////////////////////////////////////////////////////////
  
  for (ih=0;ih<nh;ih++){
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
  }
  
  memset( (void *) model[0], (int) '\0', nq * nt *FSIZE);  

  radonwtcgls_beam2(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj,quantil,nq1,nq2,rtmethod1,rtmethod2,depth1,depth2,qmin1,qmin2,factor1,factor2,fmax1,fmax2);

  radoninv_beam(h,nh,data,t,nt,dt,model,q,nq,nq1,nq2,qmin1,qmin2,rtmethod1,
		     rtmethod2,depth1,depth2,factor1,factor2,fmax1,fmax2,ffilter,amps);   

  /////////////////////////////

  for (ih=0;ih<nh;ih++){
    for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
  }    

  free1float(Wd);
  free1float(dtemp);

  return;

}
























