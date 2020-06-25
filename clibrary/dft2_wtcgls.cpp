#include "dft2.h"

/*

Daniel Trad - December - 2000
*/

#define LOOKFAC	1	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

void fftAlongTime(int sign,float **d,complex  **m, int nh, int nt, float dt, int *nf0){
   int it, ih;
   register float *rt;	/* real trace				*/
   register complex *ct;/* complex transformed trace		*/
   int nfft;		/* transform length			*/
   int nf;		/* number of frequencies		*/
   float d1;		/* output sample interval in Hz		*/   
   complex czero; czero.r=czero.i=0;
   nfft = npfar(nt); // returns nt <= nfft <= 2*nt
   //nfft = npfar0(nt, LOOKFAC * nt); // returns nt <= nfft <= 2*nt
   fprintf(stderr,"nt=%d nfft=%d\n",nt,nfft);
   if (nfft >= SU_NFLTS || nfft >= PFA_MAX) err("Padded nt=%d--too big", nfft);
   nf = nfft/2 + 1;
   d1 = 1.0/(nfft*dt);
   if (nf > nt) err("nf must <= nt\n");

   rt = ealloc1float(nfft);
   ct = ealloc1complex(nf);

   fprintf(stderr,"In fftgo nfft=%d, nf=%d nh=%d nt=%d\n",nfft,nf,nh,nt);   

   for (ih=0;ih<nh;ih++){
       for (it=0;it<nt;it++) rt[it]=d[ih][it];
       for (it=nt;it<nfft;it++) rt[it]=0;
       pfarc(sign, nfft, rt, ct);
       //for (it=0;it<nf;it++) ct[it]/=nfft;
       for (it=0;it<nf;it++) m[it][ih]=ct[it];
       for (it=nf;it<nt;it++) m[it][ih]=czero;       
   }
   *nf0=(2*nf); // number of frequencies 
   fprintf(stderr,"nf0=%d\n",*nf0);   
   
   free1complex(ct);
   free1float(rt);
   
   return;

}



void dft_matrix2(complex **F,float *h,float *k,int nh,int nk){
     int ih, ik;
     complex  co;
     complex  dco;
     complex phase, dphase;
     float dk=k[1]-k[0];
     
     for (ih=0;ih<nh;ih++){
       phase.r=dphase.r=0;
       phase.i=(h[ih]*(k[0]-dk));
       dphase.i=(h[ih]*dk);
       co=exp(phase);
       dco=exp(dphase);
       for (ik=0;ik<nk;ik++){
	 co*=dco;
	 F[ih][ik]=conjg(co);
       }
     }
  	      
     return;
}


void dft2_op(float **d, complex **m2, float* t, float* pos, float* k, int nt, int nh, int nk, float fmax, float& df){
  int ih, it, ik, freq, maxfreq, nfreq, nf0;
  float dt = t[1]-t[0];
  
  float w  = 0, wa = 0;
  complex **F= 0;
  complex **d2=0;
  //float kmin= -0.25;
  //float kmax= 0.25;
  //float dk = 1./(nk*20);
  //float kmin = (-nk/2+1)*dk;
  
  //float dk = (kmax-kmin)/(nk-1);
  
  //float dk = 0.0001;
  //for (ik=0;ik<nk;ik++) k[ik] = (kmin+ik*dk)*2*PI;
  fprintf(stderr,"npfa=%d\n",npfa(nk));
  F  = ealloc2complex(nh,nt);
  d2 = ealloc2complex(nh,nt);
  fftAlongTime(-1,d,d2,nh,nt,dt,&nf0);
  dft_matrix2(F,pos,k,nh,nk);  
  
  fprintf(stderr,"In dft2_op nh=%d, nt=%d nf0=%d\n",nh,nt,nf0);
  nfreq=nf0/2;
  df=1/(nf0*dt);

  floatprint(df);
  maxfreq=(int) (fmax/df);
  fprintf(stderr,"maxfreq=%d, pi=%f, dt=%f, df=%f\n",maxfreq,PI,dt,df);
  
  for (freq=1;freq<maxfreq;freq++){
    w=2*PI*freq*df;
    //fprintf(stderr,"w=%f\n",w);
    wa=freqweight(freq,df,fmax-10,fmax);
    Atimesx(d2[freq],F,m2[freq],nh,nk,1);
  }  

  return;

}



void dft2_wtcgls(float **d, complex **m2, float *t, float *pos, float *k,
	    int nt, int nh, int nk, float **Wd, inv_par inv, float fmax, int *pnf0)
{
   int ih, ik, freq, maxfreq, nfreq,nf0, iter;
   complex **d2=0, czero;
   complex **F=0, **FH=0;
   float w=0,wa,*dh=0,df;
   complex *maux=0;
   float dt=t[1]-t[0];
   complex *qcg, *q1cg, *scg, *x1cg, *zcg, *z1cg, *rcg, *Azcg;
   float  *eta, *rho, *gcv;
   float *Wdh;
   float **Wm;
   float *mabs;
   float quantil=inv.eps1;

   czero.r=czero.i=0;
   d2=ealloc2complex(nh,nt);

   F=ealloc2complex(nk,nh);
   FH=ealloc2complex(nh,nk);
   dh=ealloc1float(nh);
   maux=ealloc1complex(nk);
   mabs=ealloc1float(nk);
   Wdh=ealloc1float(nh);

   /*Vectors required for wtcgls*/
   qcg=ealloc1complex(nk);
   q1cg=ealloc1complex(nk);
   scg=ealloc1complex(nk);
   x1cg=ealloc1complex(nk);
   zcg=ealloc1complex(nk);
   z1cg=ealloc1complex(nk);
   rcg=ealloc1complex(nh);
   Azcg=ealloc1complex(nh);
   eta=ealloc1float(nk);
   rho=ealloc1float(nk);
   gcv=ealloc1float(nk);
   
   for (ih=1;ih<nh-1;ih++) dh[ih]=(pos[ih+1]-pos[ih-1])/2;   
   dh[0]=pos[1]-pos[0];
   dh[nh-1]=pos[nh-1]-pos[nh-2];
   float dhav=fabs(pos[nh-1]-pos[0])/(nh-1);
   for (ih=0;ih<nh;ih++){
     dh[ih]=MIN(dh[ih],dhav);
     fprintf(stderr," dh[%d]=%f\n",ih,dh[ih]);
   }
   

   fftgo0(-1,d,d2,nh,nt,dt,&nf0);
   fprintf(stderr,"In hrfft2 nh=%d, nt=%d eps1=%f nf0=%d\n",nh,nt,inv.eps1,nf0);
   *pnf0=nf0;

   nfreq=nf0/2;
   df=1/(nf0*dt);
   floatprint(df);
   maxfreq=(int) (fmax/df);
   fprintf(stderr,"maxfreq=%d, pi=%f, dt=%f, df=%f\n",maxfreq,PI,dt,df);
   Wm=ealloc2float(nk,maxfreq);

   for (freq=1;freq<maxfreq;freq++){

     w=2*PI*freq*df;
     //fprintf(stderr,"w=%f\n",w);
     wa=freqweight(freq,df,fmax-10,fmax);
     dft_matrix(F,pos,k,nh,nk);

     //Atimesx(maux,FH,d2[freq],nk,nh);
     //for (ik=0;ik<nk;ik++) m2[freq][ik]=maux[ik]/(nk*nh);    
     for (iter=1;iter<=inv.iter_end;iter++){
       
       if (iter>1) modelweight_inv(m2[freq],nk,inv.norm,inv.eps1,Wm[freq]);
       else modelweight_inv(m2[freq],nk,3,inv.eps1,Wm[freq]);
       
       //     tempint=ctoeplitz(nk,rtoep,m2[freq],maux,ftoep,gtoep);
       wtcgls(d2[freq],F,m2[freq],Wm[freq],dh,nh,nk,qcg,q1cg,scg,x1cg,zcg,
	      z1cg,rcg,Azcg,eta,rho,gcv,0,inv.step,inv.itercg);
       fprintf(stderr,"freq=%d iter=%d\n",freq,iter);
       //for (ik=0;ik<nk;ik++)  m2[freq][ik]=maux[ik];
     }
     /////////////////////////////////////////////////
     if ((wa<1)&&(wa>0)) for (ik=0;ik<nk;ik++)  m2[freq][ik]*=wa;
   }
   TRACE;
   for (ik=0;ik<nk;ik++) m2[0][ik]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++) for (ik=0;ik<nk;ik++) m2[freq][ik]=czero;
        
   fprintf(stderr,"w=%f\n",w);
   TRACE;

   free1float(mabs);
   free2float(Wm);
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
   free1float(Wdh);
   free1float(dh);
   free2complex(FH);
   free2complex(F);
   free2complex(d2);

}


























