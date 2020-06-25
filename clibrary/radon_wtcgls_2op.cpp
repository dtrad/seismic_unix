#include "su.h"
#include "radonhybrid.h"


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


void radon_wtcgls_2op(float **data, float *ph1, float *ph2, int nh, float *t, int nt, float dt, float **model, float *q, int nq, float *q1, int nq1, float *q2, int nq2, inv_par inv, float *Wd, int testadj, float fmax1, float fmax2)
{
  int iq,iter;
  complex czero;  czero.r=czero.i=0;
  float w=0,df;
  int freq, nf,nfft;
  int maxfreq1, maxfreq2;
  complex **m2;
  complex **d2;
  complex **L;
  complex **L1;
  complex **L2;
  float *Wm;
  float  *Jtot;
  float fmin=0; // Minimum freq to compute in Hz;
  float sigmam;
  float sigmad;
  float quantil1=inv.eps1;
  float quantil2=inv.eps2;

  
  d2=ealloc2complex(nh,nt);
  m2=ealloc2complex(nq,nt);//==> m2(nt x nq) 
  L=ealloc2complex(nq,nh);   //==> l(nh x nq)
  L1=ealloc2complex(nq1,nh);   //==> l(nq xnh)
  L2=ealloc2complex(nq2,nh);   //==> l(nq xnh)
  Jtot=ealloc1float(20);
  Wm=ealloc1float(nq);

  //for (iq=0;iq<nq;iq++) fprintf(stderr,"q[%d]=%f\n",iq,q[iq]);

  zero_vector(Jtot,inv.iter_end);

  fft_parameters(nt,dt,&nfft,&nf,&df);fmin=df;
  fftgo_xt2fx(-1,data,d2,nh,nt,dt,nfft,nf);
  
  maxfreq1=(int) (fmax1/df);  if (maxfreq1==0) maxfreq1=nf;
  maxfreq2=(int) (fmax2/df);  if (maxfreq2==0) maxfreq2=nf;

  fprintf(stderr,"maxfreq1=%d, maxfreq2=%d dt=%f, df=%f\n",maxfreq1,maxfreq2,dt,df);
  

  for (freq=1;freq<maxfreq1;freq++){
    w=2*PI*freq*df;
    radon_matrix_irrq(L1,ph1,q1,nh,nq1,w);
    radon_matrix(L2,ph2,q2,nh,nq2,w);
    radon_matrix_2op(L,L1,L2,nh,nq1,nq2);
    
     //if (testadj) test=testadj_rad_f(L,LH);
    
    for (iter=1;iter<=inv.iter_end;iter++){
      weights_inv(m2[freq],nq,inv.norm,sigmam,Wm,iter);
      if (iter==1 && freq < 10 ){
	for (iq=0;iq<nq1;iq++) Wm[iq]*=2;
	for (iq=nq1;iq<nq;iq++) Wm[iq]*=0.5;
      }
      wtcgls(d2[freq],L,m2[freq],Wm,Wd,nh,nq,0,inv.step,inv.itercg);
    }
    
    fprintf(stderr,"index=%d,freq=%f,\n",freq,freq*df);
  }
  //freqweighting(m2,nf,nq,df,fmin,fmax1);

  if (nq2>0){
    fprintf(stderr,"From now on one operator only\n");
    //////////////////////////////////////////////////////////////
    // For Ground roll,
    // Between maxfreq1 and freqmax2   we use only one operator
    complex **m2w=window(m2,0,nf-1,nq1,nq-1);
    for (freq=maxfreq1;freq<maxfreq2;freq++){
      w=2*PI*freq*df;
      radon_matrix(L2,ph2,q2,nh,nq2,w);
      
      for (iter=1;iter<=inv.iter_end;iter++){
	if (iter==2) 
	  deviations(m2w[freq],nq2,d2[freq],nh,inv.norm,quantil1,quantil2,&sigmam,&sigmad);
	weights_inv(m2w[freq],nq2,inv.norm,sigmam,Wm,iter);   
	wtcgls(d2[freq],L2,m2w[freq],Wm,Wd,nh,nq2,0,inv.step,inv.itercg);
      }
      fprintf(stderr,"index=%d,freq=%f,\n",freq,freq*df);
    }
    freqweighting(m2w,nf,nq2,df,fmin,fmax2);
   }

   fprintf(stderr,"w=%f\n",w);
   
   
   for (iq=nq1;iq<nq;iq++) m2[freq][iq]/=nq;      
   fftback_fx2xt(1,model,m2,nq,nt,dt,nfft,nf);  



   free1float(Wm);
   free1float(Jtot);
   free2complex(L2);
   free2complex(L1);
   free2complex(L);
   free2complex(m2);
   free2complex(d2);
   return;
}








































