#include "radontoep.h"

//dcomplex sin(dcomplex);
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void toepradon(float *pos, float **d, float dt,float **m, float *q, float dq,
	       float eps1, float eps2, float eps, float fmax, int nt, int nh,
	       int nq, int rtmethod, float depth)
{
   int ih, iq, freq, maxfreq, nfreq,tempint, exit, nf0, i;
   complex **m2, **d2, czero;
   complex **l, **lh;
   float w=0,wa,*dh,df,*gx,*Cd;
   complex *maux, *maux2, *daux, *daux2;
   complex *rtoep, *ftoep, *gtoep;
   float *g; /* Function of offset for moveout equation */

   czero.r=czero.i=0;

   d2=ealloc2complex(nh,nt); 
   m2=ealloc2complex(nq,nt);
   l=ealloc2complex(nq,nh);
   lh=ealloc2complex(nh,nq);
   dh=ealloc1float(nh);
   maux=ealloc1complex(nq);
   maux2=ealloc1complex(nq);
   daux=ealloc1complex(nh);
   daux2=ealloc1complex(nh);
   rtoep=ealloc1complex(nq);
   ftoep=ealloc1complex(nq);
   gtoep=ealloc1complex(nq);
   gx=ealloc1float(nh);
   Cd=ealloc1float(nh); 
   g=ealloc1float(nh);

   if (rtmethod==1) for (ih=0;ih<nh;ih++) g[ih]=pos[ih];
   else if (rtmethod==2) for (ih=0;ih<nh;ih++) g[ih]=pos[ih]*pos[ih];
   else if (rtmethod==3) for (ih=0;ih<nh;ih++) g[ih]=sqrt(pos[ih]*pos[ih]+depth*depth)-depth; 
   
   for (ih=1;ih<nh-1;ih++) dh[ih]=(pos[ih+1]-pos[ih-1])/2;   
   dh[0]=pos[1]-pos[0];
   dh[nh-1]=pos[nh-1]-pos[nh-2];
   
   for (ih=0;ih<nh;ih++) Cd[ih]=1;  // Optional 1/Covarianza diagonal matrix
   
   fprintf(stderr,"In p_stack1 nh=%d, nt=%d eps1=%f method=%d\n",
   nh,nt,eps1,rtmethod);
   fftgo0(-1,d,d2,nh,nt,dt,&nf0);
  
   nfreq=nf0/2;
   df=1/(nf0*dt);
   floatprint(df);
   maxfreq=(int) (fmax/df);
   fprintf(stderr,"maxfreq=%d, pi=%f, dt=%f, df=%f\n",maxfreq,PI,dt,df);

   for (freq=1;freq<maxfreq;freq++){
       w=2*PI*freq*df;
       wa=freqweight(freq,df,fmax-10,fmax);
       
       for (ih=0;ih<nh;ih++) daux[ih]=d2[freq][ih];
       //matrix_3(rtoep,l,lh,pos,q,nh,nq,w,dh,dq,rtmethod);
       radon_matrix(rtoep,l,lh,g,q,nh,nq,w,dh); 
       
       Atimesx(maux,lh,daux,nq,nh);
     
       for (i=0;i<nh;i++) gx[i]=pos[i]*pos[i]; 
       rtoep[0].r*=(1.+eps1); 
       tempint=ctoeplitz(nq,rtoep,maux2,maux,ftoep,gtoep);
       fprintf(stderr,"Levinson, nstep=%d\n",tempint);
           
/////////////////////////////////////////////////
       if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++)  maux2[iq]*=wa;
       for (iq=0;iq<nq;iq++) m2[freq][iq]=maux2[iq];        
   }
   for (iq=0;iq<nq;iq++)   m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++){
       for (iq=0;iq<nq;iq++)
          m2[freq][iq]=czero;
   }     
   fprintf(stderr,"w=%f\n",w);
   exit=fftback0(1,m,m2,nq,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
   fprintf(stderr,"After fftgo -1 nh=%d, nt=%d \n",nh,nt);

   free1float(g);
   free1float(Cd);
   free1float(gx);
   free1complex(gtoep);
   free1complex(ftoep);
   free1complex(rtoep);
   free1complex(daux2);
   free1complex(daux);
   free1complex(maux2);
   free1complex(maux);
   free1float(dh);
   free2complex(lh);
   free2complex(l);
   free2complex(m2);
   free2complex(d2);

}







