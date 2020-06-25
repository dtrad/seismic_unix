#include "hrfft2.h"
#include "clibrary.h"

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	

/*

Daniel Trad - June 9- 2000
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
     fft_matrix(rtoep,F,FH,pos,k,nh,nk,dh);

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


void hrfft2i(float **d, complex **m2, float *t, float *pos, float *k,
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
       fft_matrix(rtoep,F,FH,pos,k,nh,nk,dh);
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


void fft_matrix(complex *R, complex **F,complex **FH,float *h,float *k,int nh,int nk, float *dh)
{

  /******************************************************************** 
     FFT Transformation matrices. F, FH and R= top row of LH*L
     Input parameters:
     Out parameter:
     Notes:
     Daniel Trad- 22-02-99
  *********************************************************************/
   
     int ih, ik;
     complex  co;
     complex  dco;
     complex phase, dphase;
     float dk=k[1]-k[0];
     float Aperture=0.0;

     for (ih=0;ih<nh;ih++) Aperture+=dh[ih];        

     for (ih=0;ih<nh;ih++){
       phase.r=dphase.r=0;
       phase.i=(h[ih]*(k[0]-dk));
       dphase.i=(h[ih]*dk);
       co=exp(phase);
       dco=exp(dphase);
       for (ik=0;ik<nk;ik++){
	 co*=dco;
	 F[ih][ik]=conjg(co);
	 FH[ik][ih]=(1./Aperture)*dh[ih]*co;
	 
       }
     }
  	      
     for (ik=0;ik<nk;ik++){
       R[ik].r=0;
       R[ik].i=0;
       for (ih=0;ih<nh;ih++)
	 R[ik]+=FH[0][ih]*F[ih][ik]; //Top row of LL=LH*L
       //fprintf(stderr,"R[%d].r=%f, R[%d].i=%f\n",iq,R[iq].r,iq,R[iq].i);	    
     }
     
     return;
}






















