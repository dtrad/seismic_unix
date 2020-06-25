#include "su.h"
#include "clibrary.h"
//#include "eomig.h"
void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,float smute);
void toepradon(float *pos, float **d, float dt,float **m, float *q, float dq,
	       float eps1, float eps2, float eps, float fmax, int nt, int nh,
	       int nq, int rtmethod);
void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int  nh, int nq, int rtmethod);

/*
radonfreq
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/

void radonfreq(float *h, int nh, float **data ,float *t, int nt, float dt, float *vel, float eps1,float eps2, float eps, float *model, float *q, int nq, int itercg, int iter_end, float step, int lsmethod, float factor, float smute,float nmofactor,unsigned short  mute, unsigned short  nmute)
{
  int i, it, iq, nqt, ih;
  float  qmax, t0, dqt, qmaxt;
  float  dx_av, dx_max, dh;
  float *dtemp, *dtemp2;
  float **d;
  float **m;
  float fmax;
  float freq=1/(2*dt)*0.5;
  float dq;
  float qmin=q[0];
  int rtmethod=2;
  float *haux;
  float *Wd;
  int norm=0;
  int testadj=0;
  float scalemute;
  

  if ((haux=alloc1float(nh))==NULL)
     err("cannot allocate memory for haux\n");
  if ((dtemp=alloc1float(nt))==NULL)
     err("cannot allocate memory for dtemp\n");
  if ((d=alloc2float(nh,nt))==NULL)
     err("cannot allocate memory for d\n");
  if ((m=alloc2float(nq,nt))==NULL)
     err("cannot allocate memory for m\n");
  if ((dtemp2=alloc1float(nt))==NULL)
     err("cannot allocate memory for dtemp2\n"); 
  if ((Wd=alloc1float(nh))==NULL)
     err("cannot allocate memory for Wd\n"); 

  //for (ih=0;ih<nh;ih++) haux[ih]=h[ih];
  //interval(h,nh,&dx_max,&dx_av);
  //fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);

  // maximum and minimum q parameter to invert
  
  if (freq==0)
    fmax=1/(2*dt)*0.5;
  else
    fmax=freq;

  for (ih=0;ih<nh;ih++) Wd[ih]=1;

  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
  radon_param(fmax,h,nh,dx_av,qmin,&qmaxt,&qmax,&dq,nq,rtmethod,factor);
  //radon_param(fmax,h,nh,dx_av,qmin,&qmaxt,&qmax,&dq,nq,rtmethod,factor);
  fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
  fprintf(stderr,"freq max=%f,dq=%e\n", fmax,dq);
  for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 

  for (ih=0;ih<nh;ih++){
    //for (it=0;it<nt;it++) dtemp[it]=data[ih][it];
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) d[it][ih]=dtemp[it];
  }
  //for (ih=0;ih<nh;ih++) fprintf(stderr,"h[%d]=%f\n",ih,h[ih]);
  
  for (iq=0;iq<nq;iq++) 
    for(it=0;it<nt;it++)
      m[it][iq]=0;

  if (lsmethod==0) 
    toepradon(h,d,dt,m,q,dq,eps1,eps2,eps,fmax,nt,nh,nq,rtmethod);
  else if (lsmethod==1) 
    wtcgls0(h,nh,d,t,nt,dt,m,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,norm,step,testadj);
  else{
    fprintf(stderr,"Method no implemented\n");\
    return;
  }						

  if (mute){
    fprintf(stderr,"MUTING ********************************\n");
    for (iq=0;iq<nmute;iq++){
      scalemute=pow(((float) iq/nmute),2);
      //fprintf(stderr,"iq/nmute=%f\n",scalemute);
      for (it=0;it<nt;it++){
	m[it][iq]*=scalemute;
        m[it][nq-iq-1]*=scalemute;
      }
    }
  }
  if (1)
    for (ih=0;ih<nh;ih++)
      for (it=0;it<nt;it++)    
	d[it][ih]=0;

  hrrti(d,h,dt,m,q,fmax,nt,nh,nq,rtmethod);
  
  for (iq=0;iq<nq;iq++) 
    for(it=0;it<nt;it++)
      model[iq*nt+it]=m[it][iq];
  
  ////Test for NMO and INMO///
  if (0)
    for (ih=0;ih<nh;ih++)
      for (it=0;it<nt;it++)    
	data[ih][it]=d[it][ih];
  /////////////////////////////
 
  
  if (1)
    for (ih=0;ih<nh;ih++){
      for (it=0;it<nt;it++) dtemp[it]=d[it][ih];    
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
    }

  //for (ih=0;ih<nh;ih++) h[ih]=haux[ih];
  free1float(Wd);
  free1float(dtemp2);
  free1float(dtemp);
  free2float(d);
  free2float(m);
  free1float(haux);

  return;

}









