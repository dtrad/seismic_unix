#include "su.h"
#include "clibrary.h"
#include <math.h>
#include "nrutil.h"
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

float misfit(complex *a,complex *b, float eps1,float *dh, int n, int norm);

void radon_cholund(float *pos, int nh, float **d,float *t,int nt, float dt,float **m, float *q, int nq,float dq,float eps1, float eps2, float eps, float fmax,float *Wd, int itercg, int iter_end, int norm, float step, int testadj, int rtmethod)
{
   int ih, iq, freq, maxfreq, nfreq, exit, nf0,i,j,iter,nn;
   complex **m2, **d2, czero;
   complex **l, **lh, *Qp, **ll2;
   float w,wa,*dh,df,*Jtot,Jmod,Jdata,powd;
   float **U, **V, *W,*RHS,*x, wmax, wmin;
   double  *diagll2;
   complex *maux, *maux2, *daux, *daux2;
   complex *rtoep;
   int costflag=0;
   int freqflag=0;
   FILE *myfile3;
   czero.r=czero.i=0;
   nn=2*nh;

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
   if ((Jtot=alloc1float(20))==NULL)
     err("cannot allocate memory for Jtot\n");
   if ((Qp=alloc1complex(nq))==NULL)
     err("cannot allocate memory for Qp\n");
   if ((ll2=alloc2complex(nh,nh))==NULL)
     err("cannot allocate memory for ll2\n");
   if ((diagll2=alloc1double(nh))==NULL)
     err("cannot allocate memory for diagll2\n");
   if ((rtoep=alloc1complex(nq))==NULL)
     err("cannot allocate memory for rtoep\n");
   
   U=matrix(1,nn,1,nn);
   V=matrix(1,nn,1,nn);
   W=vector(1,nn);
   RHS=vector(1,nn);
   x=vector(1,nn);

   //dh[0]=0;
   for (ih=1;ih<nh;ih++) dh[ih]=pos[ih]-pos[ih-1];
   dh[0]=dh[1];
   // Method 7 uses Cholesky and only can be used for regular offset
   // Reasonable results can be obtained just avoiding the offset weights
   for (ih=0;ih<nh;ih++) dh[ih]=1; 
   for (iter=0;iter<iter_end;iter++) Jtot[iter]=0;  
   fftgo0(-1,d,d2,nh,nt,dt,&nf0);
   nfreq=(nf0-2)/2+1;
   df=1/(nfreq*dt*2);
   floatprint(df);
   maxfreq=(int) ((fmax/df)-1);
   if (freqflag==1) maxfreq=nfreq;
   fprintf(stderr,"maxfreq=%d\n",maxfreq);
   const double  pi=acos(-1.);

   for (freq=1;freq<maxfreq;freq++){
       w=2*pi*freq*df;
       wa=freqweight(freq,df,fmax-10,fmax);
      
       for (ih=0;ih<nh;ih++) daux[ih]=d2[freq][ih]; 
       matrix_3(rtoep,l,lh,pos,q,nh,nq,w,dh,dq,rtmethod);
       powd=rcdot(nh,daux,daux);
       Atimesx(maux,lh,daux,nq,nh); // Initial model; the adjoint
       fprintf(stderr,"Cholesky.HR2.......Index =%d,powm=%f\n",freq,rcdot(nq,maux,maux));
       iter=0;
       while (iter<iter_end){
         iter++;

	 if (iter==1) eps2=modgrad(maux,nq,norm,1.,Qp,4);
         if (iter!=1) eps2=modgrad(maux2,nq,norm,1.,Qp,4); 
         
	 for (iq=0;iq<nq;iq++) {Qp[iq].r=1;Qp[iq].i=0;}  //dc can not be recovered  
	 
         //for (i=0;i<nh;i++) 
	 //  for (j=0;j<nh;j++)
         //    if (j>i) ll2[i][j]=rtoep[j-i];
	 //    else ll2[i][j]=conjg(rtoep[i-j]);
	 AtimesBm(ll2,l,lh,nh,nq,nh,Qp); // 't' ==> only top triang is comp

	 for (i=0;i<nh;i++) 
	   ll2[i][i].r+=(eps1/100*ll2[i][i].r+eps);
	 for (i=0;i<nh;i++) 
	   for (j=0;j<nh;j++){ 
	     U[i+1][j+1]=ll2[i][j].r;
	     U[nh+i+1][j+1]=ll2[i][j].i;
	     U[i+1][nh+j+1]=-ll2[i][j].i;
	     U[nh+i+1][nh+j+1]=ll2[i][j].r;
           }
	 
	 
         for (i=0;i<nh;i++){
	   RHS[i+1]=daux[i].r;
           RHS[nh+i+1]=daux[i].i;
         }
       	 choldc(U,nn,W);
	 cholsl(U,nn,W,RHS,x);
	 
         for (i=0;i<nh;i++){
	   daux2[i].r=x[i+1];
	   daux2[i].i=x[nh+i+1];           
         } 
	 
 	 Atimesx(maux2,lh,daux2,nq,nh);
         xtimesy(maux2,maux2,Qp,nq);       
         if (costflag==1){
	   //Atimesx(daux2,l,maux2,nh,nq);
	   //Jmod=modnorm(maux2,eps2,(int) 1,nq);
	   //Jdata=misfit(daux2,daux,powd,dh,nh,(int) 10); // keep 10=Cauchy 1=l1
           //Jtot[iter-1]+=log(Jmod+Jdata);
	   //fprintf(stderr,"Iter %d, Jd=%e, Jm=%e, J=%e, eps2=%e\n",
	   //	   iter,Jdata,Jmod,(Jdata+Jmod),eps2);
         }
       }
       /////////////////////////////////////////////////
       if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++)  maux2[iq]*=wa;
       for (iq=0;iq<nq;iq++) m2[freq][iq]=maux2[iq];
       
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
   
   free_vector(x,1,nn);
   free_vector(RHS,1,nn);
   free_vector(W,1,nn);
   free_matrix(V,1,nn,1,nn);
   free_matrix(U,1,nn,1,nn);  
   free1complex(rtoep);  
   free1double(diagll2);
   free2complex(ll2);
   free1complex(Qp);
   free1float(Jtot);
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








