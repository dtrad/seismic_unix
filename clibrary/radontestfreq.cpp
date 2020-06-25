#include "su.h"
#include "clibrary.h"
//#include "eomig.h"
int taper(float **data, int nt, int nh, int ntaper,int flag);
void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,
	 float smute);

void toepradon(float *pos, float **d, float dt,float **m, float *q, float dq,
	       float eps1, float eps2, float eps, float fmax, int nt, int nh,
	       int nq, int rtmethod);

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int  nh, int nq, int rtmethod);

void wtcgls0_td(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, int rtmethod);

void wtcgls0(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, int rtmethod);

void radon_cholund(float *pos, int nh, float **d,float *t,int nt, float dt,float **m, float *q, int nq,float dq,float eps1, float eps2, float eps, float fmax,float *Wd, int itercg, int iter_end, int norm, float step, int testadj, int rtmethod);

void radonwtcgls_tfd(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float **Wd, int itercg, int iter_end,int norm,float step, int testadj, float quantil, int method);

void radon_mpcgnec(float *pos, int nh, float **d,float *t,int nt, float dt,float **m, float *q, int nq,float dq,float eps1, float eps2, float eps, float fmax,float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil, int rtmethod);

void radoncg_levinson(float *pos, int nh, float **data, float *t,int nt,float dt,float **model, float *q, int nq, float dq,float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil, int rtmethod);

void radoncg_levinson_alias(float *pos, int nh, float **data, float *t,int nt,float dt,float **model, float *q, int nq, float dq,float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil, int rtmethod);

void radoncg_levinson_alias(float *pos, int nh, float **data, float *t,int nt,float dt,float **model, float *q, int nq, float dq,float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil, int rtmethod);

void radon_mpcgnec_alias(float *pos, int nh, float **d,float *t,int nt, float dt,float **m, float *q, int nq,float dq,float eps1, float eps2, float eps, float fmax,float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil, int rtmethod);

/*
radonfreq
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void radonfreqint(float *h, int nh, float **data ,float *t, int nt, float dt, float *vel, float eps1,float eps2, float eps, float **model, float *q, int nq, int itercg, int iter_end, float step, int lsmethod, float factor, float smute,float nmofactor,unsigned short  mute, float parmute, float *h2, int nh2, float quantil, int pseudohyp, float depth)
{
  int i, it, iq, nqt, ih;
  float  qmax, t0, dqt, qmaxt;
  float  dx_av, dx_max, dh;
  float *dtemp;
  float fmax;
  float freq=1/(2*dt)*0.9;
  float dq;
  float qmin=q[0];
  int rtmethod=2;
  float *Wd;
  int norm=0;
  int testadj=0;
  float scalemute;
  int nmute=0;
  int pp=0;
  float *posh;  // Pseudo offset for pseudohyperbolic for h[]
  float *posh2; // Pseudo offset for pseudohyperbolic for h2[]
  float *htemp;
  float *h2temp;
  
  if ((dtemp=alloc1float(nt))==NULL)
     err("cannot allocate memory for dtemp\n");
  if ((Wd=alloc1float(nh))==NULL)
     err("cannot allocate memory for Wd\n"); 
  
  if (freq==0) fmax=1/(2*dt)*0.5;
  else fmax=freq;

  for (ih=0;ih<nh;ih++) Wd[ih]=1;

  //////////////////////////////////////////////////////////////////////
  // The following lines are for pseudo hyperbolic RT
  // All methods are the same with the only difference that the offsets are
  // transformed to pseudo offsets.

  if (pseudohyp){
    fprintf(stderr,"Using Pseudo hyperbolic RT#######################\n");
    rtmethod=1;
    posh=ealloc1float(nh);  // Input offset
    posh2=ealloc1float(nh2);// Output offset
    for (i=0;i<nh;i++) posh[i]=sqrt(h[i]*h[i]+depth*depth)-depth;
    for (i=0;i<nh;i++) posh2[i]=sqrt(h2[i]*h2[i]+depth*depth)-depth;
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
  
  for (iq=0;iq<nq;iq++) for(it=0;it<nt;it++) model[iq][it]=0;

  // For pseudo hyperbolic replace h with pos and h2 with posh2
  // keeping the original offset h and h2 for the final process
  if (pseudohyp){
    htemp=ealloc1float(nh);
    h2temp=ealloc1float(nh2);
    for (ih=0;ih<nh;ih++) htemp[ih]=h[ih];
    for (ih=0;ih<nh2;ih++) h2temp[ih]=h2[ih];
    for (ih=0;ih<nh;ih++) h[ih]=posh[ih];
    for (ih=0;ih<nh2;ih++) h2[ih]=posh2[ih];
  }

  if (lsmethod==0) 
    toepradon(h,data,dt,model,q,dq,eps1,eps2,eps,fmax,nt,nh,nq,rtmethod);
  else if (lsmethod==1) 
    wtcgls0(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,
	       norm,step,testadj,rtmethod);
  else if (lsmethod==2) 
    wtcgls0_td(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,
	       norm,step,testadj,rtmethod);
  else if (lsmethod==3) 
    radoncg_levinson(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj,quantil,rtmethod);
  else if (lsmethod==4) 
    radon_mpcgnec(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj,quantil,rtmethod);
  else if (lsmethod==5){
    float **Wdt=ealloc2float(nt,nh);
    for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) Wdt[ih][it]=1;
    radonwtcgls_tfd(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wdt,itercg,iter_end,norm,step,testadj,quantil,rtmethod);
    free2float(Wdt);
  }
  else if (lsmethod==6) 
    radoncg_levinson_alias(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj,quantil,rtmethod);
  else if (lsmethod==7) 
    radon_mpcgnec_alias(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj,quantil,rtmethod);
  else{
    fprintf(stderr,"Method no implemented\n");\
    return;
  }						
  if (mute){
    iq=0; while(q[iq]<parmute) iq++; nmute=iq;
    fprintf(stderr,"MUTING at nmute=%d************************\n",nmute);
    taper(model,nt,nq,nq-nmute,2);
  }
  
  if (1) for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) data[ih][it]=0;

  hrrti(data,h2,dt,model,q,fmax,nt,nh2,nq,rtmethod);
  
  // Before applying inverse nmo recover the offset and free auxiliar arrays
  if (pseudohyp){
    for (ih=0;ih<nh;ih++) h[ih]=htemp[ih];
    for (ih=0;ih<nh2;ih++) h2[ih]=h2temp[ih];
    free1float(htemp);
    free1float(h2temp);
    free1float(posh);
    free1float(posh2);
  }     

  /////////////////////////////
  if (1)
    for (ih=0;ih<nh2;ih++){
      for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
      nmo(data[ih],dtemp,t,nmofactor*h2[ih],vel,1,nt,dt,smute);  
    }

  free1float(Wd);
  free1float(dtemp);

  return;

}
























