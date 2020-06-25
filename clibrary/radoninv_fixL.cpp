#include "radoninv.h"

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int nh, int nq, int rtmethod, float depth)
{
  /**************************************************************
       INVERSE RADON TRANSFORM IN FREQ DOMAIN
         Daniel Trad
         E-mail: dtrad@geop.ubc.ca 	  	  	          
  **************************************************************/
   int ih, iq, freq, nf, maxfreq, nfft;
   float df, w=0, dq;
   complex **L, **d2, **m2, czero;
   float *g;  
   float dq0=q[1]-q[0];
   float qmin=q[0];
   float qmax;
   int nqf; // nq at every frequency
   float omegaunit=1;
   
   czero.r=czero.i=0;

   d2=ealloc2complex(nh,nt);
   m2=ealloc2complex(nq,nt);
   L=ealloc2complex(nq,nh);
   g=ealloc1float(nh);

   fprintf(stderr,"In p_stack1 nh=%d, nt=%d nq=%d rtmethod=%d\n",nh,nt,nq,rtmethod);

   radon_moveout(pos,g,nh,rtmethod,depth);          

   fft_parameters(nt,dt,&nfft,&nf,&df);
   fftgo_xt2fx(-1,m,m2,nq,nt,dt,nfft,nf);

   radon_matrix(L,g,q,nh,nq,omegaunit);
  
   if (fmax==0) maxfreq=(int) (0.9*nf);
   else maxfreq=(int) (fmax/df);
 
   fprintf(stderr,"maxfreq=%d, pi=%f, dt=%f, df=%f\n",maxfreq,PI,dt,df);

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
     nqf=1;while (q[nqf-1]<qmax) nqf++;
     if (0) fprintf(stderr,"w=%f,q[0]=%e***dq=%e***q[%d]=%e\n",w,q[0],dq,nqf-1,q[nqf-1]);
     Atimesx(d2[freq],L,m2[freq],nh,nqf);
   }

   for (ih=0;ih<nh;ih++) d2[0][ih]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nf;freq++)
	 for (ih=0;ih<nh;ih++) d2[freq][ih]=czero;
   fprintf(stderr,"w=%f\n",w);

   fftback_fx2xt(1,d,d2,nh,nt,dt,nfft,nf); 

   free1float(g);   
   free2complex(L);
   free2complex(m2);
   free2complex(d2);     

   return;
} 
 




















