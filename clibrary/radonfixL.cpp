#include "radonsolver.h"
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

void radonsolver(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, int rtmethod, float depth, char *solver)
{
  int iq,ih,freq,iter;    // counters
  int nfft,maxfreq,nf;  // sizes
  complex **m2, **d2;     // Workable arrays 
  complex **L;            // operator
  float w=0,wa,*dh;
  float df;
  float  *Wm;
  float **J;
  // Foster and Mosher offset function
  float *g=0;  
  float quantil1=eps1;
  float quantil2=eps2;
  float sigmam;
  float sigmad;   
  int freqflag=0;  
  complex czero; czero.r=czero.i=0;
  float J0=0;
  int nf2=0;
  complex *RC=0;
  float *Cd;
  // Changed to apply fixed L
  // Original dq computed for the lowest frequency is saved in dq0
  float dq0=q[1]-q[0];
  float qmin=q[0];
  float qmax;
  int nqf; // nq at every frequency
  int nf2f; // 2nd dimension for RC at every freq
  float omegaunit=1;

   // Note that c allocates memory per column
   // Hence 2 D array with alloc2type requires reversing dimension
   // i.e., first columns, then rows.
   
   d2=ealloc2complex(nh,nt);
   m2=ealloc2complex(nq,nt);
   L=ealloc2complex(nq,nh); 
   dh=ealloc1float(nh);
   Wm=ealloc1float(nq);
   g=ealloc1float(nh);
   Cd=ealloc1float(nh);
   /****    cgfft   ****/
   nf2=npfa((int) 2*nq);     
   RC=ealloc1complex(nf2);
   /*********************/
   
   dataweigths(pos,nh,Wd,TRUE);

   //for (ih=0;ih<nh;ih++) Cd[ih]=Wd[ih]*nh/fabs(pos[nh-1]-pos[0]);
   for (ih=0;ih<nh;ih++) Cd[ih]=1;
   //for (ih=0;ih<nh;ih++) Wd[ih]=Wd[ih]*nh/fabs(pos[nh-1]-pos[0]);
   radon_moveout(pos,g,nh,rtmethod,depth);  
   fft_parameters(nt,dt,&nfft,&nf,&df);

   fftgo_xt2fx(-1,data,d2,nh,nt,dt,nfft,nf);
   maxfreq=(int) (fmax/df);
   if (freqflag==1) maxfreq=nf;

   fprintf(stderr,"maxfreq=%d, dt=%f, df=%f, nfft=%d nf=%d \n",maxfreq,dt,df,nfft,nf);
   J=ealloc2float(iter_end,maxfreq);

   radon_matrix(L,g,q,nh,nq,omegaunit);   
   
   // Compute the qmax at the maxfreq
   w=2*PI*fmax;
   dq=dq0/w; // New dq is the original dq0 for w=1 divided for dq
   qmax=qmin+(nq-1)*dq; 
   fprintf(stderr,"qmin=%e***dq=%e***qmax=%e\n",qmin,dq,qmax);

   for (freq=1;freq<maxfreq;freq++){
     w=2*PI*freq*df;
     dq=dq0/w; // New dq is the original dq0 for w=1 divided for dq
     q[0]=qmin; for (iq=1;iq<nq;iq++) q[iq]=q[iq-1]+dq;
     // The qmax is the same for all frequencies, hence we need to take 
     // part of L. 
     // Compute the nq at every frequency.
     nqf=4;while (q[nqf-1]<qmax) nqf++;
     if (0) fprintf(stderr,"w=%f,q[0]=%e***dq=%e***q[%d]=%e\n",w,q[0],dq,nqf-1,q[nqf-1]);
     wa=freqweight(freq,df,fmax-10,fmax);

     //fprintf(stderr,"freq=%d\n",freq); 
     // CAREFUL: the following line has to be change for fix L
     nf2f=npfa((int) 2*nqf);     
     if (STREQ(solver,"cgfft_")) radon_matrix_cgfft(L,RC,nh,nqf,nf2f,Cd);

     for (iter=1;iter<=iter_end;iter++){

       if (iter==2) 
	 deviations(m2[freq],nq,d2[freq],nh,norm,quantil1,quantil2,&sigmam,&sigmad);
       weights_inv(m2[freq],nq,norm,sigmam,Wm,iter);
       
       //taper(Wm,nq,nq/5);
       

       if (STREQ(solver,"adj___")) 
	 Atimesx(d2[freq],L,m2[freq],nh,nqf,TRUE);       
       else if (STREQ(solver,"cgfft_"))
	 /*I change for now step to sigmam to test Wm */
	 J0=radon_cgfft(d2[freq],L,RC,m2[freq],nh,nqf,nf2f,Wm,Cd,eps,itercg,step);
       else if (STREQ(solver,"wtcgls"))
	 J0=wtcgls(d2[freq],L,m2[freq],Wm,Wd,nh,nqf,0,step,itercg);
       else if (STREQ(solver,"toep__"))
	 J0=radon_toeplitz(d2[freq],L,m2[freq],eps1,nh,nqf);      
       else if (STREQ(solver,"cgtoep"))
	 J0=radon_cgtoep(d2[freq],L,m2[freq],eps1,nh,nqf);      

       J[freq][iter-1]=J0; 

     }


     for (iq=nqf;iq<nq;iq++) m2[freq][iq]=czero;
     if ((STREQ(solver,"toep__"))||(STREQ(solver,"adj___"))) 
       for (iq=0;iq<nq;iq++) m2[freq][iq]/=nh;
     if ((STREQ(solver,"cgfft_"))) 
       for (iq=0;iq<nq;iq++) m2[freq][iq]/=nh;

      
     if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++)  m2[freq][iq]*=wa;
     normalize(iter_end,J[freq]);
   }
   
   for (iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nf;freq++)  for (iq=0;iq<nq;iq++) m2[freq][iq]=czero;
        
   fprintf(stderr,"w=%f\n",w);
   fftback_fx2xt(1,model,m2,nq,nt,dt,nfft,nf);  
   //plotcurves(J,iter_end,freq,10,50,"Cost");    

   free2float(J);
   free1complex(RC);
   free1float(Cd);
   free1float(g);
   free1float(Wm);
   free1float(dh);
   free2complex(L);
   free2complex(m2);
   free2complex(d2);

   return;
}
























