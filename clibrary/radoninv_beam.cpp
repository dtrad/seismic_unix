#include "su.h"
#include "clibrary.h"
#include <math.h>

void polygonalFilter(float *f, float *amps, int npoly,
				int nfft, float dt, float *filter);


void radon_param_beam(float fmax1, float fmax2, float *h, int nh, float *q, 
			   float nq, float qmin1, float qmin2, float *q1, 
			   float *q2, float *ph1, float *ph2, int nq1, int nq2, 
			   float depth1, float depth2, int  rtmethod1, int rtmethod2, 
			   float factor1, float factor2,float *pdq1, float *pdq2);


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
void radoninv_beam(float *h, int nh, float **data, float *t, int nt, float dt,
		   float **model, float *q,  int nq, int nq1, int nq2, float qmin1, 
		   float qmin2, float fmax, int rtmethod1, int rtmethod2, float depth1, 
		   float depth2, float factor1, float factor2)
{
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,i,j,k;
  complex **m2;
  complex **d2;
  complex czero; czero.r=czero.i=0;
  complex **L;
  complex **L1;
  complex **L2;
  float *q1;
  float *q2;
  float *ph1;
  float *ph2;
  float w,wa,df;
  float dq1;
  float dq2;


  d2=ealloc2complex(nh,nt);
  m2=ealloc2complex(nq,nt);
  L=ealloc2complex(nq,nh);  
  L1=ealloc2complex(nq1,nh);
  L2=ealloc2complex(nq2,nh);
  ph1=ealloc1float(nh);
  ph2=ealloc1float(nh);
  q1=ealloc1float(nq1);  
  q2=ealloc1float(nq2);

  radon_moveout(h,ph1,nh,rtmethod1,depth1);
  radon_moveout(h,ph2,nh,rtmethod2,depth2);

  radon_param_2op(fmax,fmax,h,nh,q,nq,qmin1,qmin2,q1,q2,nq1,nq2,depth1,
		   depth2,rtmethod1,rtmethod2,factor1,factor2,&dq1,&dq2);


  fft_parameters(nt,dt,&nfft,&nf,&df);
  fftgo_xt2fx(-1,model,m2,nq,nt,dt,nfft,nf);

  maxfreq1=(int) (fmax1/df);  if (maxfreq1==0) maxfreq1=nf;
  maxfreq2=(int) (fmax2/df);  if (maxfreq2==0) maxfreq2=nf;

  fprintf(stderr,"maxfreq1=%d, maxfreq2=%d dt=%f, df=%f\n",maxfreq1,maxfreq2,dt,df);

  for (freq=1;freq<maxfreq;freq++){
    w=2*PI*freq*df;
    matrix_irreq(L1,ph1,q1,nh,nq1,w,dq1);
    matrix(L2,ph2,q2,nh,nq2,w,dq2);
    radon_matrix_2op(L,L1,L2,nh,nq1,nq2);
 
    Atimesx(d2[freq],L,m2[freq],nh,nq,FALSE);
    fprintf(stderr,":");
    
    //for (ih=0;ih<nh;ih++) d2[freq][ih]/=nh;
  }
  if (nq2>0){ // Between maxfreq1 and freqmax2   we use only one operator
    fprintf(stderr,"From now on one operator only\n");
    for (freq=maxfreq1;freq<maxfreq2;freq++){
      w=2*PI*freq*df;
      radon_matrix(L2,ph2,q2,nh,nq2,w);
      Atimesx(d2[freq],L,m2[freq],nh,nq,FALSE);
      fprintf(stderr,".");      
      fprintf(stderr,"freq=%d\n",freq);
    }
   }

  fprintf(stderr,"freq=%d\n",freq*df);      

  freqweighting(d2,nf,nh,df,fmin,fmax2);
  fftback_fx2xt(1,data,d2,nh,nt,dt,nfft,nf); 


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


void radoninv_beam(float *h, int nh, float **data, float *t, int nt, float dt,
		   float **model, float *q,  int nq, int nq1, int nq2, float qmin1, 
		   float qmin2,  int rtmethod1, int rtmethod2, float depth1, 
		   float depth2, float factor1, float factor2, float fmax1, float fmax2,
  		   float *f, float *amps)
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
  float dq1;
  float dq2;

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

  radon_param_beam(fmax1,fmax2,h,nh,q, nq, qmin1,qmin2,q1,q2,ph1,ph2,nq1,nq2,depth1,
		   depth2,rtmethod1, rtmethod2,factor1,factor2,&dq1,&dq2);

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



  if (0){
    int temp=(int) (30/df+0.5);
    for (iq=0;iq<nq1;iq++) for(freq=0;freq<temp;freq++) m2[freq][iq]=czero;
  }

  if (0){
    int temp=(int) (30/df+0.5);
    for(freq=0;freq<temp;freq++){
      wa=freqweight(freq,df,30-10,30);
      if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++) m2[freq][iq]*=wa;
    }    
    for (iq=0;iq<nq1;iq++) for(freq=temp;freq<nfreq;freq++) m2[freq][iq]=czero;
    for (iq=nq1;iq<nq;iq++) for(freq=0;freq<nfreq;freq++) m2[freq][iq]=czero;
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
    //for (ih=0;ih<nh;ih++) d2[freq][ih]/=nh;
    if ((wa<1)&&(wa>0)){
      for (ih=0;ih<nh;ih++) d2[freq][ih]*=wa;
      fprintf(stderr,"wa=%f,freq=%d\n",wa,freq);
    }
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



void polygonalFilter(float *f, float *amps, int npoly,
				int nfft, float dt, float *filter)
/*************************************************************************
polygonalFilter -- polygonal filter with sin^2 tapering
**************************************************************************
Input:
f		array[npoly] of frequencies defining the filter
amps		array[npoly] of amplitude values
npoly		size of input f and amps arrays
dt		time sampling interval
nfft		number of points in the fft

Output:
filter		array[nfft] filter values
**************************************************************************
Notes: Filter is to be applied in the frequency domain
**************************************************************************
Author:  CWP: John Stockwell   1992
*************************************************************************/
#define PIBY2   1.57079632679490
{
        int *intfr;             /* .... integerizations of f		*/
        int icount,ifs;		/* loop counting variables              */
	int taper=0;		/* flag counter				*/
        int nf;                 /* number of frequencies (incl Nyq)     */
        int nfm1;               /* nf-1                                 */
        float onfft;            /* reciprocal of nfft                   */
        float df;               /* frequency spacing (from dt)          */

        
	intfr=alloc1int(npoly);

        nf = nfft/2 + 1;
        nfm1 = nf - 1;
        onfft = 1.0 / nfft;

        /* Compute array of integerized frequencies that define the filter*/
        df = onfft / dt;
        for(ifs=0; ifs < npoly ; ++ifs) {
                intfr[ifs] = NINT(f[ifs]/df);
                if (intfr[ifs] > nfm1) intfr[ifs] = nfm1;
        }

	/* Build filter, with scale, and taper specified by amps[] values*/
	/* Do low frequency end first*/
	for(icount=0; icount < intfr[0] ; ++icount) 
		filter[icount] = amps[0] * onfft;

	/* now do the middle frequencies */
	for(ifs=0 ; ifs<npoly-1 ; ++ifs){
	   if(amps[ifs] < amps[ifs+1]) {	
		++taper;
		for(icount=intfr[ifs]; icount<=intfr[ifs+1]; ++icount) {
		    float c = PIBY2 / (intfr[ifs+1] - intfr[ifs] + 2);
		    float s = sin(c*(icount - intfr[ifs] + 1));
		    float adiff = amps[ifs+1] - amps[ifs];
		    filter[icount] = (amps[ifs] + adiff*s*s) * onfft;
		}
	   } else if (amps[ifs] > amps[ifs+1]) {	
		++taper;
		for(icount=intfr[ifs]; icount<=intfr[ifs+1]; ++icount) {
			   float c = PIBY2 / (intfr[ifs+1] - intfr[ifs] + 2);
                	   float s = sin(c*(intfr[ifs+1] - icount + 1));
			   float adiff = amps[ifs] - amps[ifs+1];
                	   filter[icount] = (amps[ifs+1] + adiff*s*s) * onfft;
		  }
	   } else 
		if(!(taper)){
		for(icount=intfr[ifs]; icount <= intfr[ifs+1]; ++icount)
		   	   filter[icount] = amps[ifs] * onfft;
		} else {
		for(icount=intfr[ifs]+1; icount <= intfr[ifs+1]; ++icount)
		   	   filter[icount] = amps[ifs] * onfft;
		}
	}

	/* finally do the high frequency end */
	for(icount=intfr[npoly-1]+1; icount<nf; ++icount){
		filter[icount] = amps[npoly-1] * onfft;
	}

}































