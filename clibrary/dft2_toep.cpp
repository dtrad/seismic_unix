#include "dft2.h"

/*

Daniel Trad - December - 2000
*/

void dft2_toep(float **d, complex **m2, float *t, float *pos, float *k,
	    int nt, int nh, int nk, float eps1, float fmax, int *pnf0)
{
   int ih, ik, freq, maxfreq, nfreq,tempint, nf0;
   complex **d2=0, czero;
   complex **F=0, **FH=0;
   float w=0,wa,*dh=0,df;
   complex *maux=0;
   complex *rtoep=0, *ftoep=0, *gtoep=0;
   float dt=t[1]-t[0];
   czero.r=czero.i=0;
   d2=alloc2complex(nh,nt);

   F=alloc2complex(nk,nh);
   FH=alloc2complex(nh,nk);
   dh=alloc1float(nh);
   maux=alloc1complex(nk);
   rtoep=alloc1complex(nk);
   gtoep=alloc1complex(nk);
   ftoep=alloc1complex(nk);

   for (ih=1;ih<nh-1;ih++) dh[ih]=(pos[ih+1]-pos[ih-1])/2;   
   dh[0]=pos[1]-pos[0];
   dh[nh-1]=pos[nh-1]-pos[nh-2];
   


   fftgo0(-1,d,d2,nh,nt,dt,&nf0);
   fprintf(stderr,"In hrfft2 nh=%d, nt=%d eps1=%f nf0=%d\n",nh,nt,eps1,nf0);
   *pnf0=nf0;

   nfreq=nf0/2;
   df=1/(nf0*dt);
   floatprint(df);
   maxfreq=(int) (fmax/df);
   fprintf(stderr,"maxfreq=%d, pi=%f, dt=%f, df=%f\n",maxfreq,PI,dt,df);
   const double  pi=acos(-1.);
   for (freq=1;freq<maxfreq;freq++){

     w=2*pi*freq*df;
     //fprintf(stderr,"w=%f\n",w);
     wa=freqweight(freq,df,fmax-10,fmax);
     dft_matrix(rtoep,F,FH,pos,k,nh,nk,dh);

     Atimesx(maux,FH,d2[freq],nk,nh);
     rtoep[0].r*=(1.+eps1); 

     tempint=ctoeplitz(nk,rtoep,m2[freq],maux,ftoep,gtoep);
     fprintf(stderr,"Levinson, nstep=%d\n",tempint);
     //for (ik=0;ik<nk;ik++)  m2[freq][ik]=maux[ik];

     /////////////////////////////////////////////////
     if ((wa<1)&&(wa>0)) for (ik=0;ik<nk;ik++)  m2[freq][ik]*=wa;

   }

   for (ik=0;ik<nk;ik++) m2[0][ik]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++) for (ik=0;ik<nk;ik++) m2[freq][ik]=czero;
        
   fprintf(stderr,"w=%f\n",w);
   
   free1complex(gtoep);
   free1complex(ftoep);
   free1complex(rtoep);
   free1complex(maux);
   free1float(dh);
   free2complex(FH);
   free2complex(F);
   free2complex(d2);

}


























