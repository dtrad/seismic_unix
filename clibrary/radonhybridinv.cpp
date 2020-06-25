#include "su.h"
#include "radonhybridinv.h"

/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void radonhybridinv(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, int nq1, int nq2, float smute, float nmofactor, int rtmethod1, int rtmethod2, float depth1, float depth2, float fmax1, float fmax2, int filtering, int npoly, float *ffilter, float *amps)
{
  int it, ih, iq;
  float *dtemp;
  int testadj=0;
  float dq1=0;
  float dq2=0;
  float *q1;
  float *q2;
  float *ph1;
  float *ph2; 
  

  dtemp=alloc1float(nt);
  q1=ealloc1float(nq1);  
  q2=ealloc1float(nq2);
  ph1=ealloc1float(nh);
  ph2=ealloc1float(nh);

  fprintf(stderr,"fmax1=%f,fmax2=%f\n",fmax1,fmax2);
   
  //////////////////////////////////////////////////////////////////////
  
  memset( (void *) data[0], (int) '\0', nh * nt *FSIZE);  

  radon_moveout(h,ph1,nh,rtmethod1,depth1);
  radon_moveout(h,ph2,nh,rtmethod2,depth2);
  
  memcpy(q1,q,nq1*FSIZE);
  memcpy(q2,&q[nq1],nq2*FSIZE);
  
  if (filtering==1) for (iq=0;iq<nq1;iq++) for(it=0;it<nt;it++) model[iq][it]=0;
  if (filtering==2) for (iq=nq1;iq<nq;iq++) for(it=0;it<nt;it++) model[iq][it]=0;
  
  radoninv_2op(data,ph1,ph2,nh,t,nt,model,q,nq,q1,nq1,q2,nq2,fmax1,fmax2,filtering,npoly,ffilter,amps);
  TRACE;
  if (0) plotgather(data[0],nt,nh,"data_before_inmo_inv");
  /////////////////////////////
  TRACE;
  for (ih=0;ih<nh;ih++){
    for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
  }    
  TRACE;
  free1float(q1);
  free1float(q2);
  free1float(ph1);
  free1float(ph2);
  free1float(dtemp);

  return;

}
























