#include "dft2.h"

/*

Daniel Trad - December - 2000
*/

void dft2_wtcgls(float **d, complex **m2, float *t, float *pos, float *k,
	    int nt, int nh, int nk, float **Wm, inv_par inv, float fmax, int *pnf0)
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
   /*Vectors required for wtcgls*/
   qcg=alloc1complex(nk);
   q1cg=alloc1complex(nk);
   scg=alloc1complex(nk);
   x1cg=alloc1complex(nk);
   zcg=alloc1complex(nk);
   z1cg=alloc1complex(nk);
   rcg=alloc1complex(nh);
   Azcg=alloc1complex(nh);
   eta=alloc1float(nk);
   rho=alloc1float(nk);
   gcv=alloc1float(nk);
   
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
   Wm=alloc2complex(maxfreq,nk);

   for (freq=1;freq<maxfreq;freq++){

     w=2*PI*freq*df;
     //fprintf(stderr,"w=%f\n",w);
     wa=freqweight(freq,df,fmax-10,fmax);
     dft_matrix(rtoep,F,FH,pos,k,nh,nk,dh);

     Atimesx(maux,FH,d2[freq],nk,nh);
     
     for (iter=1;iter<=iter_end;iter++){
       modelweight_inv(m2[freq],nk,inv.norm,inv.eps1,Wm[freq]);
       //     tempint=ctoeplitz(nk,rtoep,m2[freq],maux,ftoep,gtoep);
       wtcgls(d2[freq],F,m2[freq],Wm[freq],Wd[freq],nh,nk,qcg,q1cg,scg,x1cg,zcg,
	      z1cg,rcg,Azcg,eta,rho,gcv,0,inv.step,inv.itercg);
       fprintf(stderr,"Levinson, nstep=%d\n",tempint);
       //for (ik=0;ik<nk;ik++)  m2[freq][ik]=maux[ik];
     }
     /////////////////////////////////////////////////
     if ((wa<1)&&(wa>0)) for (ik=0;ik<nk;ik++)  m2[freq][ik]*=wa;

   }

   for (ik=0;ik<nk;ik++) m2[0][ik]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++) for (ik=0;ik<nk;ik++) m2[freq][ik]=czero;
        
   fprintf(stderr,"w=%f\n",w);

   free2complex(Wm);
   free1float(gcv);
   free1float(rho);
   free1float(eta);
   free1complex(Azcg);
   free1complex(rcg);
   free1complex(z1cg);
   free1complex(zcg);
   free1complex(x1cg);
   free1complex(scg);
   free1complex(q1cg);
   free1complex(qcg);
   
   free1complex(maux);
   free1float(dh);
   free2complex(FH);
   free2complex(F);
   free2complex(d2);

}


























