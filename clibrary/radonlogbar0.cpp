#include "su.h"
#include "radonlogbar.h"

void radonlogbar0(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, float factor, float smute, float nmofactor, int rtmethod, float depth, float fmax, char *solver, logb_param par)
{
  int i, it, ih;
  float  qmax, qmaxt; 
  float  dx_av, dx_max;
  float *dtemp;
  float dq;
  float qmin=q[0];
  float *Wd;

  dtemp=alloc1float(nt);
  Wd=alloc1float(nh);
   
  for (ih=0;ih<nh;ih++) Wd[ih]=1;

  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
  radon_param(fmax,h,nh,dx_av,qmin,&qmaxt,&qmax,&dq,nq,rtmethod,factor);
  fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
  fprintf(stderr,"freq max=%f,dq=%e\n", fmax,dq);
  for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 
  
  
  for (ih=0;ih<nh;ih++){
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
  }
  
  memset( (void *) model[0], (int) '\0', nq * nt *FSIZE);  

  radonlogbar(h,nh,data,t,nt,dt,model,q,nq,dq,fmax,Wd,rtmethod,depth,solver,par);

  /* Any automatic filtering process can be done here on model[nq][nt] */

  hrrti(data,h,dt,model,q,fmax,nt,nh,nq,rtmethod,depth);

  for (ih=0;ih<nh;ih++){
    for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
  }    

  free1float(Wd);
  free1float(dtemp);

  return;

}
























