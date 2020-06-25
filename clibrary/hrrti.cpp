#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrary.h"

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax)
{
  //       INVERSE RADON TRANSFORM IN FREQ DOMAIN
  //       Daniel Trad
  //       E-mail: dtrad@geop.ubc.ca 	  	  	          

   int ih, iq, it, freq, nfreq, nf0, maxfreq, exit;
   float df, w, *dh, dq;
   complex **l, **lh, **d2, **m2, *daux, *maux, czero;
   extern int nt, nh, nq, rtmethod;

  // Foster and Mosher offset function
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
   if ((daux=alloc1complex(nh))==NULL)
     err("cannot allocate memory for daux\n");
   fprintf(stderr,"In p_stack1 nh=%d, nt=%d nq=%d rtmethod=%d\n"
   ,nh,nt,nq,rtmethod);
   dq=q[1]-q[0];
   dh[0]=0;
   for (ih=1;ih<nh;ih++)
       dh[ih]=pos[ih]-pos[ih-1];

   if ((gFM=alloc1float(nh))==NULL)
     err("cannot allocate memory for g\n");
   for (it=0;it<nh;it++) gFM[it]=sqrt(pos[it]*pos[it]+depth*depth)-depth;
            
   fftgo(-1,m,m2,nq,nt,dt,&nf0);

   //fprintf(stderr,"In p_stack2 nh=%d, nt=%d nf0=%d\n",nh,nt,nf0);   
   nfreq=(nf0-2)/2+1;
   df=1/(nfreq*dt*2);
   if (fmax==0) maxfreq=nfreq;
   else maxfreq=(int) ((fmax/df)-1);
   //maxfreq=nfreq;
   fprintf(stderr,"maxfreq=%d, pi=%f, dt=%f, df=%f\n",maxfreq,PI,dt,df);

   for (freq=1;freq<maxfreq;freq++){
       w=2*PI*freq*df;
       for (iq=0;iq<nq;iq++) maux[iq]=m2[freq][iq];
       matrix_3(l,pos,q,nh,nq,w,rtmethod,gFM);
       Atimesx(daux,l,maux,nh,nq);
       for (ih=0;ih<nh;ih++) d2[freq][ih]=daux[ih];
   }
   for (ih=0;ih<nh;ih++) d2[0][ih]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++)
	 for (ih=0;ih<nh;ih++)
	   d2[freq][ih]=czero;
   fprintf(stderr,"w=%f\n",w);
   exit=fftback(1,d,d2,nh,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
   fprintf(stderr,"After fftgo -1 nh=%d, nt=%d \n",nh,nt);

   free1float(gFM);   
   free1complex(daux);
   free1complex(maux);
   free1float(dh);
   free2complex(lh);
   free2complex(l);
   free2complex(m2);
   free2complex(d2);     

   return;
} 
 







