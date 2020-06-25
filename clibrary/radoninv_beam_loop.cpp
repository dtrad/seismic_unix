#include "su.h"
#include "clibrary.h"
#include <math.h>

void matrix(complex **L, float *h, float *q,int nh,int nq,float w,float dq,int rtmethod);

void polygonalFilter(float *f, float *amps, int npoly,
				int nfft, float dt, float *filter);

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

void radoninv_beam_loop(float *h, int nh, float **data, float *t, int nt, float dt,
			float **model, float *q,  int nq1, int nq2, int rtmethod1, 
			int rtmethod2, float depth1, float depth2, float fmax1, 
			float fmax2, float *f, float *amps)
{
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,i,j,k;
  complex **m2;
  complex **d2;
  complex czero;
  complex **L;
  complex **L1;
  complex **L2;
  float *q1;
  float *q2;
  float *ph1;
  float *ph2;
  float w,wa,df;
  float dq1, dq2;
  int nq=nq1+nq2;

  czero.r=czero.i=0;

  d2=ealloc2complex(nh,nt);
  
  m2=ealloc2complex(nq,nt);
  
  L=ealloc2complex(nq,nh);  
  
  L1=ealloc2complex(nq1,nh);
  
  L2=ealloc2complex(nq2,nh);

  ph1=ealloc1float(nh);

  ph2=ealloc1float(nh);

  q1=ealloc1float(nq1);  

  q2=ealloc1float(nq2);

  if (rtmethod1==3) for (ih=0;ih<nh;ih++) ph1[ih]=sqrt(h[ih]*h[ih]+depth1*depth1)-depth1;
  if (rtmethod2==3) for (ih=0;ih<nh;ih++) ph2[ih]=sqrt(h[ih]*h[ih]+depth2*depth2)-depth2;

  for (iq=0;iq<nq1;iq++) q1[iq]=q[iq];
  for (iq=0;iq<nq2;iq++) q2[iq]=q[nq1+iq];

  dq1=q1[2]-q1[1];
  dq2=q2[2]-q2[1];

  fprintf(stderr,"dq1=%f;dq2=%f\n",dq1,dq2);
  //for (iq=0;iq<nq;iq++) fprintf(stderr,"q[%d]=%e\n",iq,q[iq]);

  ////////////////////////////////////////////////////


  fftgo0(-1,model,m2,nq,nt,dt,&nf0);
  nfreq=nf0/2;
  df=1/(nf0*dt);

  int filtering=1;

  if (filtering){
    int npoly=4;
    int nfft=nf0-2;
    float *filter;
    // float *f;
    // float *amps;
    int pp;

    filter = ealloc1float(nfreq);
    // f = ealloc1float(npoly);    
    // amps = ealloc1float(npoly);   
    fprintf(stderr,"nfft=%d\n",nfft);
    // f[0]=1;
    // f[1]=15;
    // f[2]=20;
    // f[3]=70;

    int k;
    for (k=0;k<npoly;k++) fprintf(stderr,"f[%d]=%f\n",k,f[k]);
    for (k=0;k<npoly;k++) fprintf(stderr,"amps[%d]=%f\n",k,amps[k]);


    amps[0]*=nfft;
    amps[1]*=nfft;
    amps[2]*=nfft;
    amps[3]*=nfft;

    polygonalFilter(f,amps,npoly,nfft,dt,filter);
   

    //for(freq=0;freq<nfreq;freq++) fprintf(stderr,"filter[%d]=%f\n",freq,filter[freq]);
    
    for (iq=0;iq<nq1;iq++) for(freq=0;freq<nfreq;freq++) m2[freq][iq]*=filter[freq];

    free1float(filter);
    //free1float(f);
    //free1float(amps);

  }


  floatprint(df);
  float fmax=fmax1;
  maxfreq=(int) (fmax/df);
  fprintf(stderr,"maxfreq=%d, dt=%f, df=%f\n",maxfreq,dt,df);
  const double  pi=acos(-1.);
  float test; 

  for (freq=1;freq<maxfreq;freq++){
    w=2*pi*freq*df;
    wa=freqweight(freq,df,fmax-10,fmax);
    if (rtmethod1==3) matrix(L1,ph1,q1,nh,nq1,w,dq1,rtmethod1);
    else  matrix(L1,h,q1,nh,nq1,w,dq1,rtmethod1);
    if (rtmethod2==3) matrix(L2,ph2,q2,nh,nq2,w,dq2,rtmethod2);
    else matrix(L2,h,q2,nh,nq2,w,dq2,rtmethod2);
 
    for (ih=0;ih<nh;ih++){
      for (iq=0;iq<nq1;iq++)
	L[ih][iq]=L1[ih][iq];
      for (iq=0;iq<nq2;iq++) 
	L[ih][iq+nq1]=L2[ih][iq];
    }

    Atimesx(d2[freq],L,m2[freq],nh,nq,0);      
    fprintf(stderr,".");
    
    /////////////////////////////////////////////////
    //    for (ih=0;ih<nh;ih++) fprintf(stderr,"d2[%d][%d]=%f\n",freq,ih,d2[freq][ih]);
    // for (ih=0;ih<nq;ih++) fprintf(stderr,"d2[%d][%d]=%f\n",freq,ih,m2[freq][ih]);
    if ((wa<1)&&(wa>0)) for (ih=0;ih<nh;ih++) d2[freq][ih]*=wa;
    //fprintf(stderr,"wa=%f,freq=%d\n",wa,freq);
  }
  fprintf(stderr,"freq=%d\n",freq);      
  for (ih=0;ih<nh;ih++) d2[0][ih]=czero;
  fprintf(stderr,"freq=%d\n",freq);      
  for (freq=maxfreq;freq<nfreq;freq++){
    for (ih=0;ih<nh;ih++) d2[freq][ih]=czero;
  }     
  fprintf(stderr,"w=%f\n",w);
  
  exit=fftback0(1,data,d2,nh,nt,dt,nf0);
  if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
  

  free1float(q1);
  free1float(q2);
  free1float(ph1);
  free1float(ph2);
  free2complex(L2);
  free2complex(L1);
  free2complex(L);
  free2complex(m2);
  free2complex(d2);
  return;
}































