#include "radoninv.h"




/*******************************************************************************
radoninv0

Input: is a Radon model 
Output: is the data gather
Inverse NMO is performed after the Radon transform

Daniel Trad - September 20- 2000
********************************************************************************/


void radoninv0(float *h, int nh, float **data ,float *t, int nt, float dt, float *vel, float **model, float *q, int nq, float smute,float nmofactor, float depth, float fmax,int rtmethod)
{
  int i, it, iq, nqt, ih;
  float *dtemp;
  float dq;
  float *posh;  // Pseudo offset for pseudohyperbolic for h[]
  float *htemp;
  
  dtemp=ealloc1float(nt);
  
  zero_array(data,nh,nt);

  hrrti(data,h,dt,model,q,fmax,nt,nh,nq,rtmethod,depth);
  
  if (nmofactor)
    for (ih=0;ih<nh;ih++){
      xequaly(dtemp,data[ih],nt);
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
    }

  free1float(dtemp);

  return;

}

























