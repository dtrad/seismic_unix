#include "radonnlcg.h"
/* 
  Radon loop for calling nlcg
*/

void radonsolvernlcg(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, int rtmethod, float depth, char *solver)
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

   for (ih=0;ih<nh;ih++) Cd[ih]=Wd[ih]*nh/fabs(pos[nh-1]-pos[0]);
   //for (ih=0;ih<nh;ih++) Wd[ih]=Wd[ih]*nh/fabs(pos[nh-1]-pos[0]);
   radon_moveout(pos,g,nh,rtmethod,depth);  
   fft_parameters(nt,dt,&nfft,&nf,&df);

   fftgo_xt2fx(-1,data,d2,nh,nt,dt,nfft,nf);
   maxfreq=(int) (fmax/df);
   if (freqflag==1) maxfreq=nf;

   fprintf(stderr,"maxfreq=%d, dt=%f, df=%f, nfft=%d nf=%d \n",maxfreq,dt,df,nfft,nf);
   J=ealloc2float(iter_end,maxfreq);

   
   
   for (freq=1;freq<maxfreq;freq++){
     w=2*PI*freq*df;
     wa=freqweight(freq,df,fmax-10,fmax);
     radon_matrix(L,g,q,nh,nq,w);
     fprintf(stderr,"freq=%d\n",freq);
     if (STREQ(solver,"adj___")) 
       Atimesx(d2[freq],L,m2[freq],nh,nq,TRUE);       
     else if (STREQ(solver,"nlcg__"))
       J0=nlcg_interface(d2[freq],L,m2[freq],nh,nq,eps1,itercg,eps2);      
     J[freq][0]=J0; 
     if (!STREQ(solver,"wtcgls")) for (iq=0;iq<nq;iq++) m2[freq][iq]/=nq;
     if (STREQ(solver,"adj___")) for (iq=0;iq<nq;iq++) m2[freq][iq]/=nh;
      
     if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++)  m2[freq][iq]*=wa;
     normalize(iter_end,J[freq]);
   }
   
   for (iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nf;freq++)  for (iq=0;iq<nq;iq++) m2[freq][iq]=czero;
        
   fprintf(stderr,"w=%f\n",w);
   fftback_fx2xt(1,model,m2,nq,nt,dt,nfft,nf);  
   plotcurves(J,iter_end,freq,10,50,"Cost");    

   free1float(Cd);
   free2float(J);
   free1float(g);
   free1float(Wm);
   free1float(dh);
   free2complex(L);
   free2complex(m2);
   free2complex(d2);

   return;
}
























