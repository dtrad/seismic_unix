#include "radoninv.h"

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int nh, int nq, int rtmethod, float depth)
{
  /**************************************************************
       INVERSE RADON TRANSFORM IN FREQ DOMAIN
         Daniel Trad
         E-mail: dtrad@geop.ubc.ca 	  	  	          
  **************************************************************/
   int ih, freq, nf, maxfreq, nfft;
   float df, w=0, dq;
   complex **L, **d2, **m2, czero;
   float *g;  
   
   czero.r=czero.i=0;

   d2=ealloc2complex(nh,nt);
   m2=ealloc2complex(nq,nt);
   L=ealloc2complex(nq,nh);
   g=ealloc1float(nh);

   fprintf(stderr,"In radoninv nh=%d, nt=%d nq=%d rtmethod=%d\n",nh,nt,nq,rtmethod);

   dq=q[1]-q[0];
     
   radon_moveout(pos,g,nh,rtmethod,depth);          

   fft_parameters(nt,dt,&nfft,&nf,&df);
   fftgo_xt2fx(-1,m,m2,nq,nt,dt,nfft,nf);
  
   if (fmax==0) maxfreq=(int) (0.9*nf);
   else maxfreq=(int) (fmax/df);
 
   fprintf(stderr,"maxfreq=%d, pi=%f, dt=%f, df=%f\n",maxfreq,PI,dt,df);

   for (freq=1;freq<maxfreq;freq++){
       w=2*PI*freq*df;
       radon_matrix(L,g,q,nh,nq,w);
       Atimesx(d2[freq],L,m2[freq],nh,nq);
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
 




















