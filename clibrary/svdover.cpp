#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrary.h"
#include <math.h>
#include "nrutil.h"
#include "nr.h"

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void svdover(float *pos, float **d, float dt,float **m, float *q, float dq,
	     float eps1, float eps2, float eps, float fmax)
{
   int ih, iq, freq, maxfreq, nfreq, exit, nf0, i,j, iter, nnq, nnh;
   complex **m2, **d2,czero, tempcomp;
   complex **l, **lh, **ll, *Qp, **ll0;
   float w, wa, *dh,df,*Jtot,Jmod,Jdata,powd,*Cd, **U, **V, *W,*RHS,*x, 
         wmin, wmax;
   double *diagll;
   complex *maux, *maux2, *daux, *daux2;
   extern int nt, nh, nq, method, iter_end, rtmethod, norm,freqflag,costflag;
   FILE *myfile3;                
    czero.r=czero.i=0;
   nnq=2*nq;
   nnh=2*nh;
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
   if ((diagll=alloc1double(nq))==NULL)
     err("cannot allocate memory for diagll\n");
   if ((daux=alloc1complex(nh))==NULL)
     err("cannot allocate memory for daux\n");
   if ((daux2=alloc1complex(nh))==NULL)
     err("cannot allocate memory for daux\n");
   if ((Jtot=alloc1float(20))==NULL)
     err("cannot allocate memory for Jtot\n");
   if ((Qp=alloc1complex(nq))==NULL)
     err("cannot allocate memory for Qp\n");
   if ((Cd=alloc1float(nh))==NULL)
     err("cannot allocate memory for Cd\n");
   //if ((U=alloc2float(nn,nn))==NULL)
   //  err("cannot allocate memory for U\n");
   U=matrix(1,nnh,1,nnq);
   V=matrix(1,nnq,1,nnq);
   W=vector(1,nnq);
   RHS=vector(1,nnh);
   x=vector(1,nnq);

   for (ih=0;ih<nh;ih++) Cd[ih]=1.;
   //for (ih=1;ih<nh;ih++) dh[ih]=pos[ih]-pos[ih-1];
   for (ih=1;ih<nh;ih++) dh[ih]=1;
   dh[0]=dh[1];
   for (iter=0;iter<iter_end;iter++) Jtot[iter]=0.;
   fprintf(stderr,"In p_stack1 nh=%d, nt=%d eps1=%f method=%d dh[0]=%f\n",
   nh,nt,eps1,method,dh[0]);
   fftgo(-1,d,d2,nh,nt,dt,&nf0);
   
   nfreq=(nf0-2)/2+1;
   df=1/(nfreq*dt*2);
   floatprint(df);
   maxfreq=(int) ((fmax/df)-1);
   if (freqflag==1) maxfreq=nfreq;
   fprintf(stderr,"maxfreq=%d, pi=%f, dt=%f, df=%f\n",maxfreq,PI,dt,df);
   const double  pi=acos(-1.);
   for (freq=2;freq<maxfreq;freq++){
       w=2*pi*freq*df;
       wa=freqweight(freq,df,fmax-10,fmax);
       
       for (ih=0;ih<nh;ih++) daux[ih]=d2[freq][ih];
       
       matrix_3(l,lh,pos,q,nh,nq,w,dh,dq,rtmethod);
       powd=rcdot(nh,daux,daux);
      
         for (i=0;i<nh;i++) 
	   for (j=0;j<nq;j++){ 
	     U[i+1][j+1]=l[i][j].r;
	     U[nh+i+1][j+1]=l[i][j].i;
	     U[i+1][nq+j+1]=-l[i][j].i;
	     U[nh+i+1][nq+j+1]=l[i][j].r;
           }

 
         for (i=0;i<nh;i++){
	   RHS[i+1]=daux[i].r;
           RHS[nh+i+1]=daux[i].i;
         }

         svdcmp2(U,nnh,nnq,W,V);

         /*Method 1
         wmax=0.0;for (j=1;j<=nnq;j++) if (W[j]>wmax) wmax=W[j];
         wmin=eps1/100.*wmax;
         for (j=1;j<=nnq;j++) if (W[j]<wmin) W[j]=0.;
         */
         //Method 2
      	 wmax=0.0;for (j=1;j<=nnq;j++) if (W[j]>wmax) wmax=W[j];
         wmin=eps1/100.*wmax;
         for (j=1;j<=nnq;j++){
	   W[j]=W[j]/(W[j]*W[j]+wmin); // Yilmaz version
	   W[j]=1./W[j];
	 }         
         /*Method 3
      	 wmax=0.0;for (j=1;j<=nnq;j++) if (W[j]>wmax) wmax=W[j];
         wmin=eps1/100.*wmax;
         for (j=1;j<=nnq;j++) 
	   if (W[j]<wmin) {
	     W[j]=W[j]/(W[j]*W[j]+wmin);W[j]=1./W[j]; // Tad version
	   }
	 */        
         svbksb2(U,W,V,nnh,nnq,RHS,x);

         for (i=0;i<nq;i++){
	   maux2[i].r=x[i+1];
	   maux2[i].i=x[nq+i+1];           
         } 
    
/////////////////////////////////////////////////
   if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++) maux2[iq]*=wa ;
   for (iq=0;iq<nq;iq++) m2[freq][iq]=maux2[iq];
   }
   for (iq=0;iq<nq;iq++)  m2[0][iq]=czero;  //dc can not be recovered
   for (iq=0;iq<nq;iq++)  m2[1][iq]=czero;  //dc can not be recovered    
   for (freq=maxfreq;freq<nfreq;freq++){
       for (iq=0;iq<nq;iq++)
          m2[freq][iq]=czero;
   }     
   fprintf(stderr,"w=%f\n",w);
   exit=fftback(1,m,m2,nq,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
   fprintf(stderr,"After fftgo -1 nh=%d, nt=%d \n",nh,nt);
   if((myfile3=fopen("Jtottemp","w"))==NULL)
        err("cannot open file=%s\n","Jtottemp");    
   fwrite(Jtot,sizeof(float),iter_end,myfile3);
   fclose(myfile3); 

free_vector(x,1,nnq);
free_vector(RHS,1,nnh);
free_vector(W,1,nnh);
free_matrix(V,1,nnq,1,nnq);
free_matrix(U,1,nnh,1,nnh);   
free1float(Cd);
free1complex(Qp);
free1float(Jtot);
free1complex(daux2);
free1complex(daux);
free1double(diagll);
free1complex(maux2);
free1complex(maux);
free1float(dh);
free2complex(lh);
free2complex(l);
free2complex(m2);
free2complex(d2);
return;
}







