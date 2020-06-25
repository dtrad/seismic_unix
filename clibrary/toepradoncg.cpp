#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrary.h"
#include <math.h>
dcomplex sin(dcomplex);
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void toepradoncg(float *pos, float **d, float dt,float **m, float *q, float dq,
	     float eps1, float eps2, float eps, float fmax)
{
   int ih, iq, freq, maxfreq, nfreq, exit, nf0, i,j, tempint;
   complex **m2, **d2,czero;
   complex **l, **lh;
   float w,wa,*dh,df,*gx;
   complex *maux, *maux2, *daux, *daux2;
   complex *rtoep;
   extern int nt, nh, nq, method,rtmethod, itercg, freqflag;             
   complex *ss, *ss2, *gg, *rr; // for ctoephcg
 
   czero.r=czero.i=0;
   // Note that c allocates memory per column
   // Hence 2 D array with alloc2type requires reversing dimension
   // i.e., first columns, then rows.
   
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
   if ((gx=alloc1float(nh))==NULL)
     err("cannot allocate memory for gx\n");
   if ((ss=alloc1complex(nq))==NULL)
     err("cannot allocate memory for ss\n");
   if ((ss2=alloc1complex(nq))==NULL)
     err("cannot allocate memory for ss2\n");
   if ((gg=alloc1complex(nq))==NULL)
     err("cannot allocate memory for gg\n");
   if ((rr=alloc1complex(nq))==NULL)
     err("cannot allocate memory for rr\n");
   

   for (ih=1;ih<nh;ih++) dh[ih]=pos[ih]-pos[ih-1];
   dh[0]=dh[1];
  
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
       matrix_3(rtoep,l,lh,pos,q,nh,nq,w,dh,dq,rtmethod);

       Atimesx(maux,lh,daux,nq,nh);
       /*
       for (i=0;i<nq;i++){
          rtoep[i]=czero;
          for (j=0;j<nh;j++) 
              rtoep[i]+=lh[0][j]*l[j][i];
       }
       */
       //noise=((eps1/100)*rtoep[0].r);
       //rtoep[0].i=0;
       //fprintf(stderr,"rtoep[0]=%e, noise=%e\n",rtoep[0].r,noise.r);
       rtoep[0]*=(1+eps1);
       tempint=ctoephcg(itercg,nq,rtoep,maux2,maux,ss,ss2,gg,rr);
       fprintf(stderr,"Levinson, nstep=%d\n",tempint);
       
/////////////////////////////////////////////////
      if ((wa<1)&&(wa>0))
	for (iq=0;iq<nq;iq++)  maux2[iq]*=wa;
      for (iq=0;iq<nq;iq++)
          m2[freq][iq]=maux2[iq];

        
   }
   for (iq=0;iq<nq;iq++)     m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++){
       for (iq=0;iq<nq;iq++)
          m2[freq][iq]=czero;
   }     
   fprintf(stderr,"w=%f\n",w);
   exit=fftback(1,m,m2,nq,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
   fprintf(stderr,"After fftgo -1 nh=%d, nt=%d \n",nh,nt);
   

free1complex(rr);
free1complex(gg);
free1complex(ss2);
free1complex(ss);
free1float(gx);
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







