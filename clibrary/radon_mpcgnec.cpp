#include "su.h"
#include "clibrary.h"
#include <math.h>
#include "nrutil.h"
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void Atimesx(complex *b,complex **A,complex *x,int n,int m, int adj);

float mpcgne(complex **A, int n, int m, complex *x, complex *b,float *MI, float tol, 
	     float step,int itercg,int restart, float eps1, float eps2, int itermin);

void radon_mpcgnec(float *pos, int nh, float **d,float *t,int nt, float dt,float **m, float *q, int nq,float dq,float eps1, float eps2, float eps, float fmax,float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil, int rtmethod)
{
   int ih, iq, freq, maxfreq, nfreq, exit, nf0,i,j,iter,nn;
   complex **m2, **d2, czero;
   complex **L, **LH;
   float w,wa,*dh,df,*Jtot,Jmod,Jdata,powd;
   float wmax, wmin;
   float *Qp;
   complex *maux, *maux2, *daux, *daux2;
   complex *rtoep;
   int costflag=0;
   int freqflag=0;
   const double  pi=acos(-1.);
   int pp=0;
   int itermin;
   FILE *myfile3;

   czero.r=czero.i=0;

   if ((d2=alloc2complex(nh,(nt)))==NULL)
     err("cannot allocate memory for d2\n"); 
   if ((m2=alloc2complex(nq,(nt)))==NULL)//==> m2(nt x nq) 
     err("cannot allocate memory for m2\n");
   if ((L=alloc2complex(nq,nh))==NULL)   //==> l(nh x nq)
     err("cannot allocate memory for l\n");
   if ((LH=alloc2complex(nh,nq))==NULL)  //==> l(nq xnh)
     err("cannot allocate memory for lh\n");
   if ((dh=alloc1float(nh))==NULL)
     err("cannot allocate memory for dh\n"); 
   if ((Jtot=alloc1float(20))==NULL)
     err("cannot allocate memory for Jtot\n");
   if ((rtoep=alloc1complex(nq))==NULL)
     err("cannot allocate memory for rtoep\n");
   
   Qp=ealloc1float(nq);
   float *mtemp=ealloc1float(nq);
   fprintf(stderr,"quantil=%f\n",quantil);

   for (ih=1;ih<nh;ih++) dh[ih]=pos[ih]-pos[ih-1];
   dh[0]=dh[1];
   for (ih=0;ih<nh;ih++) dh[ih]=1; 
   for (iter=0;iter<iter_end;iter++) Jtot[iter]=0;  
   fftgo0(-1,d,d2,nh,nt,dt,&nf0);
   nfreq=(nf0-2)/2+1;
   df=1/(nfreq*dt*2);
   floatprint(df);
   maxfreq=(int) ((fmax/df)-1);
   if (freqflag==1) maxfreq=nfreq;
   fprintf(stderr,"maxfreq=%d\n",maxfreq);
   int freq0=5;
   for (freq=freq0;freq<maxfreq;freq++){
       w=2*pi*freq*df;
       wa=freqweight(freq,df,fmax-10,fmax);
       matrix_3(rtoep,L,LH,pos,q,nh,nq,w,dh,dq,rtmethod);
       //fprintf(stderr,"MPCGNE........Index =%d\n",freq);
       iter=0;
       
       while (iter<iter_end){
         iter++;
	 if (iter==1){ 
	   for (iq=0;iq<nq;iq++) Qp[iq]=eps1; 
	   itermin=25;
	 }
	 else{
	   for (i=0;i<nq;i++) mtemp[i]=abs(m2[freq][i]);
	   float qmodel=quest(quantil,nq,mtemp);  
	   for (iq=0;iq<nq;iq++) Qp[iq]=pow(abs(m2[freq][iq]),2)/(qmodel*qmodel)+1e-5;
	   itermin=25;
	 }
	 mpcgne(L,nh,nq,m2[freq],d2[freq],Qp,eps,step,itercg,1,eps1,eps2,itermin);
	 //	 fprintf(stderr,"AFTERpp**************=%d\n",++pp);  
       }
       /////////////////////////////////////////////////
       if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++)  m2[freq][iq]*=wa;
   }
   for (freq=0; freq<freq0;freq++)
     for (iq=0;iq<nq;iq++)     
       m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++){
     for (iq=0;iq<nq;iq++)
       m2[freq][iq]=czero;
   }     
   fprintf(stderr,"w=%f\n",w);
   
   exit=fftback0(1,m,m2,nq,nt,dt,nf0);   
   if (1){ 
     save_gather(m,nq,q,nt,dt,"mtemp.su");
     system("suxwigb < mtemp.su perc=99 title=\"mtemp\" &");     
     system("sufft < mtemp.su | suamp | suxwigb perc=99 title=\"mtemp\" &");
   }
   if (0){
     exit=fftback0(1,d,d2,nh,nt,dt,nf0);   
     save_gather(d,nh,nt,dt,"dtemp.su");
     system("suxwigb < dtemp.su perc=99 title=\"dtemp\" &");     
     system("sufft < dtemp.su | suamp | suxwigb perc=99 title=\"dtemp\" &");
   }
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
   
   if((myfile3=fopen("Jtottemp","w"))==NULL)
     err("cannot open file=%s\n","Jtottemp");    
   
   fwrite(Jtot,sizeof(float),iter_end,myfile3);
   fclose(myfile3);
   free1float(mtemp);
   free1complex(rtoep);  
   free1float(Qp);
   free1float(Jtot);
   free1float(dh);
   free2complex(LH);
   free2complex(L);
   free2complex(m2);
   free2complex(d2);
   
}








