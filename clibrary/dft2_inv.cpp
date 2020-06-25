#include "dft2.h"

/*

Daniel Trad - December - 2000

*/


void dft2_inv(float **d, complex **m2, float *t, float *pos, float *k,
	     int nt, int nh, int nk, float fmax, int nf0)
{
   int ih, freq, maxfreq, nfreq,  exit;
   complex **d2=0, czero;
   complex **F=0, **FH=0;
   float w=0,wa,*dh=0,df;
   complex *rtoep=0;
   float dt=t[1]-t[0];

   czero.r=czero.i=0;
   d2=alloc2complex(nh,2*nt);
   F=alloc2complex(nk,nh);
   FH=alloc2complex(nh,nk);
   dh=alloc1float(nh);
   rtoep=alloc1complex(nk);
   fprintf(stderr,"In hrfft2i nf0=%d\n",nf0);
   nfreq=nf0/2;
   df=1/(nf0*dt);
   floatprint(df);
   maxfreq=(int) (fmax/df);
   fprintf(stderr,"maxfreq=%d, pi=%f, dt=%f, df=%f, nt=%d nh=%d nk=%d \n",maxfreq,PI,dt,df,nt,nh,nk);
   const double  pi=acos(-1.);
   for (freq=1;freq<maxfreq;freq++){
       w=2*pi*freq*df;
       //fprintf(stderr,"w=%f freq=%d \n",w,freq);
       wa=freqweight(freq,df,fmax-10,fmax);
       dft_matrix(rtoep,F,FH,pos,k,nh,nk,dh);
          Atimesx(d2[freq],F,m2[freq],nh,nk);
/////////////////////////////////////////////////
       if ((wa<1)&&(wa>0)) for (ih=0;ih<nh;ih++)  d2[freq][ih]*=wa;
       
   }

   for (ih=0;ih<nh;ih++) d2[0][ih]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++) for (ih=0;ih<nh;ih++) d2[freq][ih]=czero;
        
   fprintf(stderr,"w=%f\n",w);
   exit=fftback0(1,d,d2,nh,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
   TRACE;

   free1float(dh);
   free2complex(FH);
   free2complex(F);
   free2complex(d2);

}


