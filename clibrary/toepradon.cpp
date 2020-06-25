#include "su.h"
#include "clibrary.h"
#include <math.h>
dcomplex sin(dcomplex);
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void toepradon(float *pos, float **d, float dt,float **m, float *q, float dq,
	     float eps1, float eps2, float eps, float fmax)
{
   int ih, iq, freq, maxfreq, nfreq,tempint, exit, nf0, i;
   complex **m2, **d2, czero;
   complex **l, **lh;
   float w,wa,*dh,df,*gx,*Cd;
   complex *maux, *maux2, *daux, *daux2;
   complex *rtoep, *ftoep, *gtoep;
   extern int nt, nh, nq, method, rtmethod, freqflag;
   extern float depth;
   float *gFM;

   czero.r=czero.i=0;

   if ((d2=alloc2complex(nh,(nt+100)))==NULL)
     err("cannot allocate memory for d2\n"); 
   if ((m2=alloc2complex(nq,(nt+100)))==NULL)//==> m2(nt x nq) 
     err("cannot allocate memory for m2\n");
   if ((l=alloc2complex(nq,nh))==NULL)   //==> l(nh x nq)
     err("cannot allocate memory for l\n");
   if ((lh=alloc2complex(nh,nq))==NULL)  //==> l(nq xnh)
     err("cannot allocate memory for lh\n");
   if ((dh=alloc1float(nh))==NULL)
     err("cannot allocate memory for dh\n"); 
   if ((maux=alloc1complex(nq))==NULL)
     err("cannot allocate memory for maux\n");
   if ((maux2=alloc1complex(nq))==NULL)
     err("cannot allocate memory for maux2\n");
   if ((daux=alloc1complex(nh))==NULL)
     err("cannot allocate memory for daux\n");
   if ((daux2=alloc1complex(nh))==NULL)
     err("cannot allocate memory for daux\n");
   if ((rtoep=alloc1complex(nq))==NULL)
     err("cannot allocate memory for rtoep\n");
   if ((ftoep=alloc1complex(nq))==NULL)
     err("cannot allocate memory for ftoep\n");
   if ((gtoep=alloc1complex(nq))==NULL)
     err("cannot allocate memory for gtoep\n");
   if ((gx=alloc1float(nh))==NULL)
     err("cannot allocate memory for gx\n");
   if ((Cd=alloc1float(nh))==NULL)
     err("cannot allocate memory for gx\n");
   if (rtmethod==3){
     if ((gFM=alloc1float(nh))==NULL)
       err("cannot allocate memory for g\n");
     for (i=0;i<nh;i++) gFM[i]=sqrt(pos[i]*pos[i]+depth*depth)-depth;
   }


   //for (ih=1;ih<nh;ih++) dh[ih]=pos[ih]-pos[ih-1];
   for (ih=1;ih<nh;ih++) dh[ih]=1.;
   dh[0]=dh[1];
   for (ih=0;ih<nh;ih++) Cd[ih]=1;  // Optional 1/Covarianza diagonal matrix
   //Cd[nh-1]=0.1;
   //Cd[nh-2]=0.3;
   //Cd[nh-3]=0.5;
   //Cd[nh-4]=0.7;        
   
   fprintf(stderr,"In p_stack1 nh=%d, nt=%d eps1=%f method=%d\n",
   nh,nt,eps1,method);
   fftgo(-1,d,d2,nh,nt,dt,&nf0);
  
   nfreq=(nf0-2)/2+1;
   df=1/(nfreq*dt*2);
   floatprint(df);
   maxfreq=(int) ((fmax/df)-1);
   if (freqflag==1) maxfreq=nfreq;
   fprintf(stderr,"maxfreq=%d, pi=%f, dt=%f, df=%f\n",maxfreq,PI,dt,df);
   const double  pi=acos(-1.);
   for (freq=1;freq<maxfreq;freq++){
       w=2*pi*freq*df;
       wa=freqweight(freq,df,fmax-10,fmax);
       
       for (ih=0;ih<nh;ih++) daux[ih]=d2[freq][ih];
       matrix_4(rtoep,l,lh,pos,q,nh,nq,w,dh,dq,rtmethod,gFM);

       Atimesx(maux,lh,daux,nq,nh);
     
       for (i=0;i<nh;i++) gx[i]=pos[i]*pos[i]; 
       /*
       for (i=0;i<nq;i++){
          rtoep[i]=czero;
          for (j=0;j<nh;j++) 
	      rtoep[i]+=(lh[0][j]*l[j][i]);
	      }*/
       rtoep[0].r*=(1.+eps1); 
       tempint=ctoeplitz(nq,rtoep,maux2,maux,ftoep,gtoep);
       //fprintf(stderr,"Levinson, nstep=%d\n",tempint);
           
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
   exit=fftback(1,m,m2,nq,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
   fprintf(stderr,"After fftgo -1 nh=%d, nt=%d \n",nh,nt);


if (rtmethod==3) free1float(gFM);
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







