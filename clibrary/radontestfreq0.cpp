#include "su.h"
#include "clibrary.h"
//#include "eomig.h"
int taper(float **data, int nt, int nh, int ntaper,int flag);
void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,float smute);
void toepradon(float *pos, float **d, float dt,float **m, float *q, float dq,
	       float eps1, float eps2, float eps, float fmax, int nt, int nh,
	       int nq, int rtmethod);
void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int  nh, int nq, int rtmethod);

void radonwtcgls_tfd(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float **Wd, int itercg, int iter_end,int norm,float step, int testadj, float quantil);

void radon_mpcgnec(float *pos, int nh, float **d,float *t,int nt, float dt,float **m, float *q, int nq,float dq,float eps1, float eps2, float eps, float fmax,float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil);

void radoncg_levinson_alias(float *pos, int nh, float **data, float *t,int nt,float dt,float **model, float *q, int nq, float dq,float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil);

void radon_mpcgnec_alias(float *pos, int nh, float **d,float *t,int nt, float dt,float **m, float *q, int nq,float dq,float eps1, float eps2, float eps, float fmax,float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil);

void radoncg_levinson_beam(float *pos, int nh, float **data, float *t,int nt,float dt,float **model, float *q, int nq, float dq,float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil);

/*
radonfreq
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/

void radonfreq(float *h, int nh, float **data ,float *t, int nt, float dt, float *vel, float eps1,float eps2, float eps, float **model, float *q, int nq, int itercg, int iter_end, float step, int lsmethod, float factor, float smute,float nmofactor,unsigned short  mute, float parmute, float quantil)
{
  int i, it, iq, nqt, ih;
  float  qmax, t0, dqt, qmaxt;
  float  dx_av, dx_max, dh;
  float *dtemp;
  float fmax;
  float freq=1/(2*dt)*0.6;
  float dq;
  float qmin=q[0];
  int rtmethod=2;
  float *Wd;
  int norm=0;
  int testadj=0;
  float scalemute;
  int nmute=0;
  int pp=0;
  if ((dtemp=alloc1float(nt))==NULL)
     err("cannot allocate memory for dtemp\n");
  if ((Wd=alloc1float(nh))==NULL)
     err("cannot allocate memory for Wd\n"); 
  
  if (freq==0)
    fmax=1/(2*dt)*0.5;
  else
    fmax=freq;

  for (ih=0;ih<nh;ih++) Wd[ih]=1;

  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
  radon_param(fmax,h,nh,dx_av,qmin,&qmaxt,&qmax,&dq,nq,rtmethod,factor);
  fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
  fprintf(stderr,"freq max=%f,dq=%e\n", fmax,dq);
  for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 

  //////////////////////////////////////////////////////////////////////
  // The following lines are for pseudo hyperbolic RT

  float *posh=ealloc1float(nh);
  float depth=1000;

  if (lsmethod==8){
    for (i=0;i<nh;i++) posh[i]=sqrt(h[i]*h[i]+depth*depth)-depth;
    radon_param(fmax,h,nh,dx_av,qmin,&qmaxt,&qmax,&dq,nq,1,factor);
    qmin=-1e-4;
    for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 
  }
  //////////////////////////////////////////////////////////////////////

  fprintf(stderr,"pp=%d******************\n",++pp);
  for (ih=0;ih<nh;ih++){
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
  }
  fprintf(stderr,"pp=%d******************\n",++pp);  
  for (iq=0;iq<nq;iq++) 
    for(it=0;it<nt;it++)
      model[iq][it]=0;

  if (lsmethod==0) 
    toepradon(h,data,dt,model,q,dq,eps1,eps2,eps,fmax,nt,nh,nq,rtmethod);
  else if (lsmethod==1) 
    wtcgls0(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj);
  else if (lsmethod==2) 
    wtcgls0_td(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj);
  else if (lsmethod==3) 
    radoncg_levinson(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj,quantil);
  else if (lsmethod==4) 
    radon_mpcgnec(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj,quantil);
  else if (lsmethod==5){
    float **Wdt=ealloc2float(nt,nh);
    for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) Wdt[ih][it]=1;
    radonwtcgls_tfd(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wdt,itercg,iter_end,norm,step,testadj,quantil);
    free2float(Wdt);
  }
  else if (lsmethod==6) 
    radoncg_levinson_alias(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj,quantil);
  else if (lsmethod==7) 
    radon_mpcgnec_alias(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj,quantil);
  else if (lsmethod==8) 
    radoncg_levinson_beam(posh,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj,quantil);
  else{
    fprintf(stderr,"Method no implemented\n");
    return;
  }						
  fprintf(stderr,"pp=%d******************\n",++pp);
  if (mute){
    iq=0; while(q[iq]<parmute) iq++; nmute=iq;
    fprintf(stderr,"MUTING at nmute=%d************************\n",nmute);
    taper(model,nt,nq,nq-nmute,2);

  }
  if (1)
    for (ih=0;ih<nh;ih++)
      for (it=0;it<nt;it++)    
	data[ih][it]=0;
  fprintf(stderr,"pp=%d******************\n",++pp);
  hrrti(data,h,dt,model,q,fmax,nt,nh,nq,rtmethod);
    
  /////////////////////////////
  if (1)
    for (ih=0;ih<nh;ih++){
      for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
    }

  free1float(Wd);
  free1float(dtemp);
  free1float(posh);


  return;

}

