#include "su.h"
#include "radonsolver.h"

/*
radonfixL0

Input is a seismic gather
Output is the same gather after Radon multiple removal.
Similar to radonsolver0 but the dq is computed for maximum frequency equal 
to freq=1/(2pi) ==> w=1
Daniel Trad - July 1- 2001
*/


void radonsolver0(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, int itercg, int iter_end, float step, float eps2, float eps1, float quantil, int norm, float factor, float smute, float nmofactor, int rtmethod, float depth, float fmax, char *solver)
{
  int i, it, ih;
  float  qmax, qmaxt;
  float  dx_av, dx_max;
  float *dtemp;
  float dq;
  float qmin=q[0];
  float *Wd;
  float eps=1e-7;
  int testadj=0;
  float frequnit;

  dtemp=alloc1float(nt);
  Wd=alloc1float(nh);
   
  for (ih=0;ih<nh;ih++) Wd[ih]=1;

  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
  // Changed for fixed L
  // The dq parameter is computed for w=1 ==> fmax=1/(2*pi)
  frequnit=1/(2*PI); // Note that radon_param uses f not w
  radon_param(frequnit,h,nh,dx_av,qmin,&qmaxt,&qmax,&dq,nq,rtmethod,factor);
  fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
  fprintf(stderr,"freq max=%f,dq=%e\n", fmax,dq);
  // First q axis to be sent to radonfixL. However it changes at every freq
  for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 
  
  
  for (ih=0;ih<nh;ih++){
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
  }
  //  plotgather(data,nh,nt,dt,"suxwigb");
  if (0) plotgather_pipe(data,nh,nt,"After NMO");
  memset( (void *) model[0], (int) '\0', nq * nt *FSIZE);  

  radonsolver(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,
	       norm,step,testadj,rtmethod,depth,solver);

  /* Any automatic filtering process can be done here on model[nq][nt] */
  
  // The q axis contains the original dq0 for frequnit
  // It is created again because it changed inside radonsolver
  for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 
  hrrti(data,h,dt,model,q,fmax,nt,nh,nq,rtmethod,depth);
  // The q axis contains the original dq0 for frequnit
  // It is created again because it changed inside radonsolver
  for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 

  for (ih=0;ih<nh;ih++){
    for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
  }    

  free1float(Wd);
  free1float(dtemp);

  return;

}
























