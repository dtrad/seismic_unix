#include "su.h"
#include "Complex.h"
#include "clibrary.h"

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int nh, int nq, int rtmethod)
{
  //       INVERSE RADON TRANSFORM IN FREQ DOMAIN
  //       Daniel Trad
  //       E-mail: dtrad@geop.ubc.ca 	  	  	          

   int ih, freq, nfreq, nf0, maxfreq, exit;
   float df, w=0, dq;
   complex **L, **d2, **m2, czero;
   // Foster and Mosher offset function
   float depth=10;
   float *g;  
   

   czero.r=czero.i=0;


   d2=ealloc2complex(nh,nt);
   m2=ealloc2complex(nq,nt);
   L=ealloc2complex(nq,nh);
   g=ealloc1float(nh);

   fprintf(stderr,"In p_stack1 nh=%d, nt=%d nq=%d rtmethod=%d\n",nh,nt,nq,rtmethod);

   dq=q[1]-q[0];

   if (rtmethod==1) for (ih=0;ih<nh;ih++) g[ih]=pos[ih];
   else if (rtmethod==2) for (ih=0;ih<nh;ih++) g[ih]=pos[ih]*pos[ih];
   else if (rtmethod==3) for (ih=0;ih<nh;ih++) g[ih]=sqrt(pos[ih]*pos[ih]+depth*depth)-depth;
     
            
   fftgo0(-1,m,m2,nq,nt,dt,&nf0);

   //fprintf(stderr,"In p_stack2 nh=%d, nt=%d nf0=%d\n",nh,nt,nf0);   
   nfreq=nf0/2;
   df=1/(nf0*dt);
   if (fmax==0) maxfreq=nfreq;
   else maxfreq=(int) (fmax/df);
 
   fprintf(stderr,"maxfreq=%d, pi=%f, dt=%f, df=%f\n",maxfreq,PI,dt,df);

   for (freq=1;freq<maxfreq;freq++){
       w=2*PI*freq*df;
       radon_matrix(L,g,q,nh,nq,w);
       //matrix_3(rtoep,l,lh,pos,q,nh,nq,w,dh,dq,rtmethod);
       Atimesx(d2[freq],L,m2[freq],nh,nq);
   }
   for (ih=0;ih<nh;ih++) d2[0][ih]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++)
	 for (ih=0;ih<nh;ih++)
	   d2[freq][ih]=czero;
   fprintf(stderr,"w=%f\n",w);
   exit=fftback0(1,d,d2,nh,nt,dt,nf0);

   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
   fprintf(stderr,"After fftgo -1 nh=%d, nt=%d \n",nh,nt);

   free1float(g);   
   free2complex(L);
   free2complex(m2);
   free2complex(d2);     

   return;
} 
 




















