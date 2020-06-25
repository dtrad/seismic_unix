#include "su.h"
#include "clibrary.h"
//#include "eomig.h"
int taper(float **data, int nt, int nh, int ntaper,int flag);
void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,
	 float smute);

void radoncg_levinson_beam(float *pos, int nh, float **data, float *t,int nt,float dt,float **model, float *q, int nq, float dq,float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil, int rtmethod);

void radonwtcgls_beam2(float *h, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, float quantil, int nq1, int nq2, int rtmethod1, int rtmethod2, float depth1, float depth2, float qmin1, float qmin2, float factor1, float factor2, float fmax1, float fmax2);

void radonwtcgls_beam(float *h, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, float quantil, int nq1, int nq2, int rtmethod1, int rtmethod2, float depth1, float depth2, float qmin1, float qmin2, float factor1, float factor2);


void radoninv_beam(float *h, int nh, float **data, float *t, int nt, float dt,
		   float **model, float *q,  int nq, int nq1, int nq2, float qmin1, 
		   float qmin2, float fmax, int rtmethod1, int rtmethod2, float depth1, 
		   float depth2, float factor1,float factor2);

void radoninv_beam(float *h, int nh, float **data, float *t, int nt, float dt,
		   float **model, float *q,  int nq, int nq1, int nq2, float qmin1, 
		   float qmin2, int rtmethod1, int rtmethod2, float depth1, 
		   float depth2, float factor1,float factor2, float fmax1, float fmax2,
		   float *f, float *amps);

/*
radonfreqint_beam
Radon transform of a gather using two different operators
d = L1 m1 + L2  m2

The shape of the two operators is given by rtmethod=1,2,3 (LRT,PRT, PseudoHyp from Foster)
Model m1 goes from q[0]=qmin1 to q[nq1]
Model m1 goes from q[nq1]=qmin2 to q[nq]

The dq values are computed using the corresponding alias condition (for LRT, PRT, PHRT)

Daniel Trad - August 9- 2000
*/


void radonfreqint_beam(float *h, int nh, float **data ,float *t, int nt, float dt, float *vel, float eps1,float eps2, float eps, float **model, float *q, int nq, int itercg, int iter_end, float step, int lsmethod, float factor, float smute, float nmofactor, unsigned short mute, float parmute, float *h2, int nh2, float quantil, float qmin1, float qmin2, int nq1, int nq2, int rtmethod1, int rtmethod2, float depth1, float depth2, float factor1, float factor2, int mute1, int mute2, float fmax1, float fmax2, float *ffilter, float *amps)
{
  int i, it, iq, nqt, ih;
  float  qmax, t0, dqt, qmaxt;
  float  dx_av, dx_max, dh;
  float *dtemp;
  float fmax;
  float freq=1/(2*dt)*0.5;
  float dq;
  float qmin=q[0];
  float *Wd;
  int norm=0;
  int testadj=0;
  float scalemute;
  int nmute=0;
  int pp=0;
  float *htemp;
  float *h2temp;

  nq2=nq-nq1; // nq1+nq2 must be = nq , because of this  nq2 given is ignored 

  // Combination method for L1 and L2
  // According to what combination method we use different parameters
  
  if ((dtemp=alloc1float(nt))==NULL)
     err("cannot allocate memory for dtemp\n");
  if ((Wd=alloc1float(nh))==NULL)
     err("cannot allocate memory for Wd\n"); 
  
  // fmax is always give for freq=fn/2 for now. 
  if (freq==0) fmax=1/(2*dt)*0.5;
  else fmax=20;

  // For now we are not really using Wd but is kept for future.
  for (ih=0;ih<nh;ih++) Wd[ih]=1; 
  
  //////////////////////////////////////////////////////////////////////
  
  for (ih=0;ih<nh;ih++){
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
  }
  
  for (iq=0;iq<nq;iq++) for(it=0;it<nt;it++) model[iq][it]=0;
  int onemaxfreq=0;
  if (lsmethod==9)
    if (onemaxfreq)
      radonwtcgls_beam(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj,quantil,nq1,nq2,rtmethod1,rtmethod2,depth1,depth2,qmin1,qmin2,factor1,factor2);
    else
      radonwtcgls_beam2(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj,quantil,nq1,nq2,rtmethod1,rtmethod2,depth1,depth2,qmin1,qmin2,factor1,factor2,fmax1,fmax2);
  else{
    fprintf(stderr,"Method no implemented\n");\
    return;
  }						

  fprintf(stderr,"nq=%d,nq1=%d,nq2=%d,mute1=%d,mute2=%d\n",nq,nq1,nq2,mute1,mute2);
  
  if(mute1) for (iq=0;iq<nq1;iq++)  for(it=0;it<nt;it++) model[iq][it]=0;



  if(mute2) for (iq=nq1;iq<nq;iq++) for(it=0;it<nt;it++) model[iq][it]=0;


  
  if (1) for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) data[ih][it]=0;

  if (onemaxfreq)
  radoninv_beam(h2,nh2,data,t,nt,dt,model,q,nq,nq1,nq2,qmin1,qmin2,fmax,rtmethod1,
		     rtmethod2,depth1,depth2,factor1,factor2); 
  else 
  radoninv_beam(h2,nh2,data,t,nt,dt,model,q,nq,nq1,nq2,qmin1,qmin2,rtmethod1,
		     rtmethod2,depth1,depth2,factor1,factor2,fmax1,fmax2,ffilter,amps); 

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
























