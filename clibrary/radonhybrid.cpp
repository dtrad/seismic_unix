#include "su.h"
#include "radonhybrid.h"

/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void radonhybrid(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, inv_par inv, float qmin1, float qmin2,  int nq1, int nq2, float factor1, float factor2, float smute, float nmofactor, int rtmethod1, int rtmethod2, float depth1, float depth2, float fmax1, float fmax2, int filtering, int npoly, float *ffilter, float *amps, int symmetric)
{
  int it, ih;
  float *dtemp;
  float *Wd;
  int testadj=0;
  float dq1=0;
  float dq2=0;
  float *q1;
  float *q2;
  float *ph1;
  float *ph2; 
  TRACE;
  dtemp=alloc1float(nt);
  Wd=alloc1float(nh);
  q1=ealloc1float(nq1);  
  q2=ealloc1float(nq2);
  ph1=ealloc1float(nh);
  ph2=ealloc1float(nh);
  
  fprintf(stderr,"fmax1=%f,fmax2=%f\n",fmax1,fmax2);
   
  for (ih=0;ih<nh;ih++) Wd[ih]=1;

  //////////////////////////////////////////////////////////////////////
  if (nmofactor>0){
    for (ih=0;ih<nh;ih++){
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
      for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
    }
  }

  memset( (void *) model[0], (int) '\0', nq * nt *FSIZE);  
  

  radon_moveout(h,ph1,nh,rtmethod1,depth1);
  radon_moveout(h,ph2,nh,rtmethod2,depth2);
  
  dataweigths(h,nh,Wd,TRUE);
  TRACE;
  radon_param_2op(fmax1,fmax2,h,nh,q, nq, qmin1,qmin2,q1,q2, nq1, nq2, depth1, depth2,
		   rtmethod1, rtmethod2, factor1, factor2,&dq1,&dq2, symmetric);
  TRACE;

  radon_wtcgls_2op(data,ph1,ph2,nh,t,nt,dt,model,q,nq,q1,nq1,q2,nq2,inv,Wd,testadj,fmax1,fmax2);
  TRACE;
  memset( (void *) data[0], (int) '\0', nh * nt *FSIZE);  

  radoninv_2op(data,ph1,ph2,nh,t,nt,model,q,nq,q1,nq1,q2,nq2,fmax1,fmax2,filtering,npoly,ffilter,amps);
  
  /////////////////////////////

  plotgather(data[0],nt,nh,"data_before_inmo");

  if (nmofactor>0){
    for (ih=0;ih<nh;ih++){
      for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
    }    
  }

  free1float(q1);
  free1float(q2);
  free1float(ph1);
  free1float(ph2);
  free1float(Wd);
  free1float(dtemp);

  return;

}
























