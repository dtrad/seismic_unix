#include "radonlogbar.h"

/* 
  Frequency loop for calling logbarrier_method in each frequency.
  Call  logbarrier_interface 
  Daniel Trad - UBC March - 2001
*/

void radonlogbar(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float fmax, float *Wd,  int rtmethod, float depth, char *solver, logb_param par)

{
  int iq,ih,freq;    // counters
  int nfft,maxfreq,nf;  // sizes
  complex **m2=0, **d2=0;     // Workable arrays 
  complex **L=0;            // operator
  float w=0,wa,*dh;
  float df;   /* frequency sampling interval */ 
  float *g=0; /* offset function */ 
  int freqflag=0;  
  complex czero; czero.r=czero.i=0;
  float *Cd=0;  /* Data Covariance */
  int iter_done=0;

   // Note that c allocates memory per column
   // Hence 2 D array with alloc2type requires reversing dimension
   // i.e., first columns, then rows.
   
   d2=ealloc2complex(nh,nt);
   m2=ealloc2complex(nq,nt);
   L=ealloc2complex(nq,nh); 
   dh=ealloc1float(nh);
   g=ealloc1float(nh);
   Cd=ealloc1float(nh);
   /*********************/

   dataweigths(pos,nh,Wd,TRUE);

   for (ih=0;ih<nh;ih++) Cd[ih]=Wd[ih]*nh/fabs(pos[nh-1]-pos[0]);

   radon_moveout(pos,g,nh,rtmethod,depth);  
   fft_parameters(nt,dt,&nfft,&nf,&df);

   fftgo_xt2fx(-1,data,d2,nh,nt,dt,nfft,nf);
   maxfreq=(int) (fmax/df);
   if (freqflag==1) maxfreq=nf;

   fprintf(stderr,"maxfreq=%d, dt=%f, df=%f, nfft=%d nf=%d \n",maxfreq,dt,df,nfft,nf);
   
   for (freq=1;freq<maxfreq;freq++){
     w=2*PI*freq*df;
     wa=freqweight(freq,df,fmax-10,fmax);
     radon_matrix(L,g,q,nh,nq,w);
     
     if (STREQ(solver,"adj___")) 
       Atimesx(d2[freq],L,m2[freq],nh,nq,TRUE);       
     else if (STREQ(solver,"logbar"))
       iter_done=logbarrier_interface(d2[freq],L,m2[freq],nh,nq,par);      
     fprintf(stderr,"freq=%d,iter_done=%d\n",freq,iter_done);
     /* Check later to scale factor */
     if (!STREQ(solver,"logbar")) for (iq=0;iq<nq;iq++) m2[freq][iq]*=nq;
     if (STREQ(solver,"adj___")) for (iq=0;iq<nq;iq++) m2[freq][iq]/=nh;
      
     if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++)  m2[freq][iq]*=wa;
   }
   
   for (iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nf;freq++)  for (iq=0;iq<nq;iq++) m2[freq][iq]=czero;
        
   fprintf(stderr,"w=%f\n",w);
   fftback_fx2xt(1,model,m2,nq,nt,dt,nfft,nf);  

   free1float(Cd);
   free1float(g);
   free1float(dh);
   free2complex(L);
   free2complex(m2);
   free2complex(d2);

   return;
}
























