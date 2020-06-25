#include "su.h"
#include "radonsolver.h"


/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void radonsolver0_mute(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, int itercg, int iter_end, float step, float eps2, float eps1, float quantil, int norm, float factor, float smute, float nmofactor, int rtmethod, float depth, float fmax, char *solver, float parmute1, float parmute2, int mute1, int mute2, float t0)
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


  dtemp=ealloc1float(nt);
  Wd=ealloc1float(nh);
   
  for (ih=0;ih<nh;ih++) Wd[ih]=1;

  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
  radon_param(fmax,h,nh,dx_av,qmin,&qmaxt,&qmax,&dq,nq,rtmethod,factor);
  fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
  fprintf(stderr,"freq max=%f,dq=%e, nh=%d, nt=%d\n", fmax,dq,nh,nt);
  for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 
  

  for (ih=0;ih<nh;ih++){
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
  }

  //  plotgather(data,nh,nt,dt,"suxwigb");
  if (0) plotgather_pipe(data,nh,nt,"After_NMO");
  memset( (void *) model[0], (int) '\0', nq * nt *FSIZE);  

  radonsolver(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,
	       norm,step,testadj,rtmethod,depth,solver);

  if (0) plotgather_pipe(model,nq,nt,"Radon_before_mute");
  if (1){
    save_gather(model,nq,q,nt,dt,"RTmodel");
    system("cat RTmodel >> radondata &");
    system("suxwigb key=f2 < RTmodel title=\"RT\" perc=99 & \n");
  }
  /* Any automatic filtering process can be done here on model[nq][nt] */
  if ((mute1)||(mute2))
    doublemute(data,model,nq,nh,nt,q,h,t,parmute1,parmute2,mute1,mute2,fmax,rtmethod,
	       depth,t0);
  else hrrti(data,h,dt,model,q,fmax,nt,nh,nq,rtmethod,depth);
  if (0) plotgather_pipe(data,nh,nt,"After_double_mute");



  for (ih=0;ih<nh;ih++){
    for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
  }    

  free1float(Wd);
  free1float(dtemp);

  return;

}
























