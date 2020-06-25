#include "su.h"
#include "clibrary.h"
#include <math.h>

void wtcgls(complex *b,complex **L, complex *x,complex *Wm,
	  float *Wd,int nh,int nq,complex *q,complex *q1,complex *s,
	  complex *x1,complex *z,complex *z1,complex *r,complex *Az,
	  float *eta,float *rho, float *gcv, float tol, float step, int itercg);

void matrix(complex **L, float *h, float *q,int nh,int nq,float w,float dq,int rtmethod);

/* This is an interface to the WTCGLS method in routine wtcgls.cpp

   This function solves the system of equations 
     (FH WdT Wd FH + WmT Wm ) m = FH WdT Wd d 
   Notice that LH=FH WdT and L= Wd F are computed with matrix_3.
      
   If we assumed noise and model are uncorrelated,
   i.e., data and model covariances matrices Cd and Cm  are diagonal 
   
   Wd is the factorization of the inverse data Covariance matrix
   inv(Cd)=WdT Wd
   
   Wm is the factorization of the inverse model Covariance matrix
   inv(Cm)=WmT Wm

   Wm is useful to increase resolution 
   Wd is useful to taper data at the boundaries and 
      to reduce effect of killed or bad traces

   Wm depends on the model so iteration are necessary
   Wd does not change. 

*/

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void radonwtcgls_beam(float *h, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, int rtmethod)
{
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,i,j,k,tempint,iter,laststep,flag;
  complex **m2;
  complex **d2;
  complex czero;
  complex **L;
  complex **L1;
  complex **L2;

  float w,wa,*dh,*dh2,df,*Jtot, Jdata, Jmod, normu;
  complex *d,*u;
  complex *dtemp, *dc, *Wm;
  complex *qcg, *q1cg, *scg, *x1cg, *zcg, *z1cg, *rcg, *Azcg;
  float  *eta, *rho, *gcv;
  float power, Jtotlast, Jtotprev, bb, tempfloat; 
  float *uaux;
  float dx_max;
  float dx_av;
  int rtmethod1=2;
  int rtmethod2=2;
  float dq1;
  float dq2;
  float qmin1;
  float qmin2;
  int nq1=nq/2;
  int nq2=nq-nq1;
  float *q1;
  float *q2;
  int costflag=0;
  int freqflag=0;  
  float factor=1;
  float qmaxt;
  float qmax;
  float *ph1;
  float *ph2;
  
  czero.r=czero.i=0;
  eps1/=100.;
  
  d2=ealloc2complex(nh,nt);
  
  m2=ealloc2complex(nq,nt);//==> m2(nt x nq) 
  
  L=ealloc2complex(nq,nh);   //==> l(nh x nq)
  
  L1=ealloc2complex(nq1,nh);   //==> l(nq xnh)
  
  L2=ealloc2complex(nq2,nh);   //==> l(nq xnh)

  q1=ealloc1float(nq1);  

  q2=ealloc1float(nq2);

  dh=ealloc1float(nh);
  
  dh2=ealloc1float(nh);
  
  Jtot=ealloc1float(20);
  
  qcg=ealloc1complex(nq); 

  q1cg=ealloc1complex(nq);

  scg=ealloc1complex(nq);
  
  x1cg=ealloc1complex(nq);
  
  zcg=ealloc1complex(nq);

  z1cg=ealloc1complex(nq);
  
  rcg=ealloc1complex(nh);
  
  Azcg=ealloc1complex(nh);
  
  eta=ealloc1float(nq);
  
  rho=ealloc1float(nq);

  gcv=ealloc1float(nq);

  Wm=ealloc1complex(nq);
  
  dc=ealloc1complex(nh);
  
  uaux=ealloc1float(nq);
  
  ph1=ealloc1float(nh);

  ph2=ealloc1float(nh);

  ////////////////////////////////////////////////////
  
  float depth1=1000;
  float depth2=1200;
  
  int pseudohyp=0;

  if (pseudohyp){  
    factor=4;
    qmin1=-1e-4;
    rtmethod1=1;
    rtmethod2=1;
  }
  else{
    qmin1=-1e-4;
    qmin2=-5e-8;
    factor=4;
    rtmethod1=1;
    rtmethod2=2;
  }

  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
    
  for (ih=0;ih<nh;ih++) ph1[ih]=sqrt(h[ih]*h[ih]+depth1*depth1)-depth1;
  for (ih=0;ih<nh;ih++) ph2[ih]=sqrt(h[ih]*h[ih]+depth2*depth2)-depth2;
 

  radon_param(fmax,h,nh,dx_av,qmin1,&qmaxt,&qmax,&dq1,nq1,rtmethod1,factor);
    
  fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
  fprintf(stderr,"freq max=%f,dq=%e\n", fmax,dq);
  
  for (iq=0;iq<nq1;iq++) q1[iq]=qmin1+iq*dq1;

  //qmin2=q1[nq1-1]+dq1;
  
  radon_param(fmax,h,nh,dx_av,qmin2,&qmaxt,&qmax,&dq2,nq1,rtmethod2,factor);

  for (iq=0;iq<nq2;iq++) q2[iq]=qmin2+iq*dq2;

  for (iq=0;iq<nq1;iq++) q[iq]=q1[iq];
  for (iq=0;iq<nq2;iq++) q[iq+nq1]=q2[iq];
  
  ///////////////////////////////////////////////////
  
  float *mtemp=ealloc1float(nq);
  
  for (ih=1;ih<nh;ih++) dh[ih]=h[ih]-h[ih-1];
  dh[0]=dh[1];

  for (i=0;i<=iter_end;i++) Jtot[i]=0.;
  fftgo0(-1,data,d2,nh,nt,dt,&nf0);
  nfreq=nf0/2;
  df=1/(nf0*dt);
  floatprint(df);
  maxfreq=(int) (fmax/df);
  if (freqflag==1) maxfreq=nfreq;
  fprintf(stderr,"maxfreq=%d, dt=%f, df=%f\n",maxfreq,dt,df);
  const double  pi=acos(-1.);
  float test; 
  for (freq=1;freq<maxfreq;freq++){
     w=2*pi*freq*df;
     wa=freqweight(freq,df,fmax-10,fmax);
     for (iq=0;iq<nq;iq++){
       Wm[iq].r=eps1;
       Wm[iq].i=0.0;
     }       
     if (pseudohyp){
       matrix(L1,ph1,q1,nh,nq1,w,dq1,rtmethod1);
       matrix(L2,ph2,q2,nh,nq2,w,dq2,rtmethod2);
     }
     else{
       matrix(L1,h,q1,nh,nq1,w,dq1,rtmethod1);
       matrix(L2,h,q2,nh,nq2,w,dq2,rtmethod2);
     }
     for (ih=0;ih<nh;ih++){
       for (iq=0;iq<nq1;iq++)
	 L[ih][iq]=L1[ih][iq];
       for (iq=0;iq<nq2;iq++) 
	 L[ih][iq+nq1]=L2[ih][iq];
     }	 
     
     //if (testadj) test=testadj_rad_f(L,LH);
     fprintf(stderr,"freq=%d\n",freq);
     power=rcdot(nh,d2[freq],d2[freq]);
     
     for (iter=1;iter<=iter_end;iter++){

       if (iter>1) for (iq=0;iq<nq;iq++) mtemp[iq]=abs(m2[freq][iq]); 
    
       for (iq=0;iq<nq;iq++)  {
	 if ((iter!=1)){
	   float qmodel=quest(.7,nq,mtemp);
	   //qmodel=1;
       	   normu=abs(m2[freq][iq])/qmodel+eps2;
	   Wm[iq].r=1/normu+eps1;
	 }
	 else if (0){	   	   
	   normu=abs(uaux[iq])+eps2;
	   Wm[iq].r=1/normu+eps1;
	 }
	 else{ 
	   Wm[iq].r=eps1;
	   Wm[iq].i=0.;
	 }
       }
       wtcgls(d2[freq],L,m2[freq],Wm,Wd,nh,nq,qcg,q1cg,scg,x1cg,zcg,
	      z1cg,rcg,Azcg,eta,rho,gcv,0,step,itercg);
     }
     for (iq=0;iq<nq;iq++) uaux[iq]=abs(m2[freq][iq]);            
     /////////////////////////////////////////////////
     for (iq=0;iq<nq;iq++) m2[freq][iq]/=nq;      
     if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++)  m2[freq][iq]*=wa;
     //for (iq=0;iq<nq;iq++) m2[freq][iq]=u[iq];          
   }
   for (iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++){
       for (iq=0;iq<nq;iq++) m2[freq][iq]=czero;
   }     
   fprintf(stderr,"w=%f\n",w);
   exit=fftback0(1,model,m2,nq,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");

   free1float(q1);
   free1float(q2);
   free1float(mtemp);
   free1float(uaux);
   free1complex(dc);
   free1complex(Wm);
   free1float(gcv);
   free1float(rho);
   free1float(eta);
   free1complex(Azcg);
   free1complex(rcg);
   free1complex(z1cg);
   free1complex(zcg);
   free1complex(x1cg);
   free1complex(scg);
   free1complex(q1cg);
   free1complex(qcg);
   free1float(Jtot);
   free1float(dh2);
   free1float(dh);
   free2complex(L2);
   free2complex(L1);
   free2complex(L);
   free2complex(m2);
   free2complex(d2);
   return;
}


void matrix(complex **L, float *h, float *q,int nh,int nq,float w, float dq, int rtmethod)
{

        register int ih;
	register int iq;  
        complex  arg;

        for (iq=0;iq<nq;iq++){
	  for (ih=0;ih<nh;ih++){
              arg.r=0;
              if (rtmethod==1) arg.i=w*h[ih]*q[iq];
              else if (rtmethod==2) arg.i=w*h[ih]*h[ih]*q[iq];
	      
	      L[ih][iq]=exp(-arg)/sqrt(nq*nh);
	  }
	}
	
        return;
}






















