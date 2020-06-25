#include "radoncgfft_tfd.h"

/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void radoncgfft0_tfd_mute(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, float fmax,  float quantil, float factor, float smute, float nmofactor, int rtmethod, float depth, inv_par inv, float parmute1, float parmute2, int mute1, int mute2, float t0)
{
  int i, it, ih;
  float  qmax, qmaxt;
  float  dx_av, dx_max;
  float *dtemp;
  float dq;
  float qmin=q[0];
  float **Wd;
  int testadj=0;

  dtemp=alloc1float(nt);
  Wd=ealloc2float(nt,nh);   

  //////////////////////////////////////////////////////////////////////

  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
  radon_param(fmax,h,nh,dx_av,qmin,&qmaxt,&qmax,&dq,nq,rtmethod,factor);
  fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
  fprintf(stderr,"freq max=%f,dq=%e\n", fmax,dq);
  for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 
  
  //////////////////////////////////////////////////////////////////////
  
  for (ih=0;ih<nh;ih++){
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
  }
  
  memset( (void *) model[0], (int) '\0', nq * nt *FSIZE);  

  for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) Wd[ih][it]=1;

  radoncgfft_tfd(h,nh,data,t,nt,dt,model,q,nq,dq,fmax,Wd,testadj,quantil,rtmethod,depth,inv);

  //hrrti(data,h,dt,model,q,fmax,nt,nh,nq,rtmethod,depth);
  /* Any automatic filtering process can be done here on model[nq][nt] */
  doublemute(data,model,nq,nh,nt,q,h,t,parmute1,parmute2,mute1,mute2,fmax,rtmethod,depth, t0);
  if (0) plotgather_pipe(data,nh,nt,"After double mute");
  /////////////////////////////

  for (ih=0;ih<nh;ih++){
    for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
  }    


  free2float(Wd);
  free1float(dtemp);
  return;
}
























