#include "su.h"
#include "radonwtcgls.h"

/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void radonwtcgls0(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, int itercg, int iter_end, float step, float eps2, float eps1, float quantil, int norm, float factor, float smute, float nmofactor, int rtmethod, float depth, float fmax)
{
  int i, it, ih;
  float  qmax, qmaxt;
  float  dx_av, dx_max;
  float *dtemp;
  float dq;
  float qmin=q[0];
  float *Wd;
  float *posh;  // Pseudo offset for pseudohyperbolic for h[]
  float *htemp;
  int pseudohyp=0;
  float eps=1e-7;
  int testadj=0;



  if (rtmethod==3) pseudohyp=1;

  dtemp=alloc1float(nt);
  posh=ealloc1float(nh);  // Input offset
  Wd=alloc1float(nh);
  htemp=ealloc1float(nh);
   
  for (ih=0;ih<nh;ih++) Wd[ih]=1;

  //////////////////////////////////////////////////////////////////////
  // The following lines are for pseudo hyperbolic RT
  // All methods are the same with the only difference that the offsets are
  // transformed to pseudo offsets.

  if (pseudohyp){
    fprintf(stderr,"Using Pseudo hyperbolic RT#######################\n");
    rtmethod=1;

    for (i=0;i<nh;i++) posh[i]=sqrt(h[i]*h[i]+depth*depth)-depth;
  }

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

  // For pseudo hyperbolic replace h with pos and h2 with posh2
  // keeping the original offset h and h2 for the final process
  if (pseudohyp){

    for (ih=0;ih<nh;ih++) htemp[ih]=h[ih];
    for (ih=0;ih<nh;ih++) h[ih]=posh[ih];
  }


  wtcgls0(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,
	       norm,step,testadj,rtmethod);
  TRACE;
  hrrti(data,h,dt,model,q,fmax,nt,nh,nq,rtmethod);
  TRACE;
  // Before applying inverse nmo recover the offset and free auxiliar arrays
  if (pseudohyp) for (ih=0;ih<nh;ih++) h[ih]=htemp[ih];


  /////////////////////////////

  for (ih=0;ih<nh;ih++){
    for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
  }    
  TRACE;
  free1float(htemp);
  free1float(posh);
  free1float(Wd);
  free1float(dtemp);

  return;

}
























