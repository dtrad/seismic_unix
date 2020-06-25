#include "su.h"
#include "clibrary.h"
#include <math.h>
#include "nrutil.h"
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

void Atimesx(float *b,float **A,float *x,int adj, int n,int m);

float dot(int n, float *a, float *b);

float mpcgne(float **A, int n, int m, float *x,float *b,float *MI, float tol, float step,int itercg,int restart);

void radon_mpcgne(float *pos, int nh, float **d,float *t,int nt, float dt,float **m, float *q, int nq,float dq,float eps1, float eps2, float eps, float fmax,float *Wd, int itercg, int iter_end, int norm, float step, int testadj)
{
   int ih, iq, freq, maxfreq, nfreq, exit, nf0,i,j,iter,nn;
   complex **m2, **d2, czero;
   complex **L, **LH;
   float w,wa,*dh,df,*Jtot,Jmod,Jdata,powd;
   float **U, *W,*RHS, *x, wmax, wmin;
   float *Qp;
   double  *diagll2;
   complex *maux, *maux2, *daux, *daux2;
   complex *rtoep;
   int rtmethod=2;
   int costflag=0;
   int freqflag=0;
   const double  pi=acos(-1.);
   int pp=0;
   FILE *myfile3;

   czero.r=czero.i=0;

   nn=2*nh;

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
   if ((diagll2=alloc1double(nh))==NULL)
     err("cannot allocate memory for diagll2\n");
   if ((rtoep=alloc1complex(nq))==NULL)
     err("cannot allocate memory for rtoep\n");
   
   Qp=ealloc1float(2*nq);
   U=ealloc2float(2*nh,2*nq);
   W=ealloc1float(2*nh);
   RHS=ealloc1float(2*nh);
   x=ealloc1float(2*nq);

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

   for (freq=1;freq<maxfreq;freq++){
       w=2*pi*freq*df;
       wa=freqweight(freq,df,fmax-10,fmax);

       matrix_3(rtoep,L,LH,pos,q,nh,nq,w,dh,dq,rtmethod);
       fprintf(stderr,"Cholesky.HR2.......Index =%d\n",freq);
       iter=0;
       while (iter<iter_end){
         iter++;
	 if (iter){
	   for (iq=0;iq<nq;iq++){
	     Qp[iq]=1;
	     Qp[nq+iq]=1;
	   }
	 }
	 else{
	   for (iq=0;iq<nq;iq++){
	     Qp[iq]=1./abs(x[iq])+1e-2;
	     Qp[nq+iq]=1./abs(x[nq+iq])+1e-2;
	   }
	 }
	 for (i=0;i<nh;i++) 
	   for (j=0;j<nq;j++){ 
	     U[i][j]=L[i][j].r;
	     U[nh+i][j]=L[i][j].i;
	     U[i][nq+j]=-L[i][j].i;
	     U[nh+i][nq+j]=L[i][j].r;
           }

         for (i=0;i<nh;i++){
	   RHS[i]=d2[freq][i].r;
           RHS[nh+i]=d2[freq][i].i;
         }

	 mpcgne(U,2*nh,2*nq,x,RHS,Qp,eps,step,itercg,1);
	 //	 fprintf(stderr,"AFTERpp**************=%d\n",++pp);  
       }
       /////////////////////////////////////////////////
       for (iq=0;iq<nq;iq++){ 
	 m2[freq][iq].r=x[iq];
	 m2[freq][iq].i=x[nq+iq];
       }
       if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++)  m2[freq][iq]*=wa;
       
   }
   for (iq=0;iq<nq;iq++)     m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++){
     for (iq=0;iq<nq;iq++)
       m2[freq][iq]=czero;
   }     
   fprintf(stderr,"w=%f\n",w);
   
   exit=fftback0(1,m,m2,nq,nt,dt,nf0);   
   if (0){ 
     save_gather(m,nq,q,nt,dt,"mtemp.su");
     system("suxwigb < mtemp.su perc=99 title=\"mtemp\" &");     
     system("sufft < mtemp.su | suamp | suxwigb perc=99 title=\"mtemp\" &");
   }
   exit=fftback0(1,d,d2,nh,nt,dt,nf0);   
   save_gather(d,nh,nt,dt,"dtemp.su");
   system("suxwigb < dtemp.su perc=99 title=\"dtemp\" &");     
   system("sufft < dtemp.su | suamp | suxwigb perc=99 title=\"dtemp\" &");
   
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
   
   if((myfile3=fopen("Jtottemp","w"))==NULL)
     err("cannot open file=%s\n","Jtottemp");    
   
   fwrite(Jtot,sizeof(float),iter_end,myfile3);
   fclose(myfile3);
   
   free1complex(rtoep);  
   free1double(diagll2);
   free1float(Qp);
   free1float(Jtot);
   free1float(dh);
   free2complex(LH);
   free2complex(L);
   free2float(U);
   free1float(x);
   free1float(RHS);
   free2complex(m2);
   free2complex(d2);
   
}








