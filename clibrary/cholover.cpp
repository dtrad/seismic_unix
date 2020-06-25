#include "su.h"
#include "clibrary.h"
#include <math.h>
#include "nrutil.h"

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void cholover(float *pos, float **d, float dt,float **m, float *q, float dq,
	     float eps1, float eps2, float eps, float fmax)
{
   int ih, iq, freq, maxfreq, nfreq, exit, nf0, i,j, iter, nn;
   complex **m2, **d2,czero, tempcomp;
   complex **l, **lh, **ll, *Qp, **ll0;
   float w, wa, *dh,df,*Jtot,Jmod,Jdata,powd,*Cd;
   float **U, **V, *W,*RHS,*x, wmax, wmin;
   double *diagll;
   complex *maux, *maux2, *daux, *daux2;
   complex *rtoep;
   extern int nt, nh, nq, method, iter_end, rtmethod, norm,freqflag,costflag;
   FILE *myfile3;                
   czero.r=czero.i=0;
   nn=2*nq;
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
   if ((ll=alloc2complex(nq,nq))==NULL)
     err("cannot allocate memory for ll\n");
   if ((ll0=alloc2complex(nq,nq))==NULL)
     err("cannot allocate memory for ll0\n"); 
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
   if ((rtoep=alloc1complex(nq))==NULL)
     err("cannot allocate memory for rtoep\n");
   if ((Qp=alloc1complex(nq))==NULL)
     err("cannot allocate memory for Qp\n");
   if ((Cd=alloc1float(nh))==NULL)
     err("cannot allocate memory for Cd\n");
   U=matrix(1,nn,1,nn);
   V=matrix(1,nn,1,nn);
   W=vector(1,nn);
   RHS=vector(1,nn);
   x=vector(1,nn);


   for (ih=0;ih<nh;ih++) Cd[ih]=1.;
   for (ih=1;ih<nh;ih++) dh[ih]=pos[ih]-pos[ih-1];
   dh[0]=dh[1];
   for (iter=0;iter<iter_end;iter++) Jtot[iter]=0.;
   fftgo(-1,d,d2,nh,nt,dt,&nf0);
   
   nfreq=(nf0-2)/2+1;
   df=1/(nfreq*dt*2);
   floatprint(df);
   maxfreq=(int) ((fmax/df)-1);
   if (freqflag==1) maxfreq=nfreq;
   fprintf(stderr,"maxfreq=%d\n",maxfreq);
   const double  pi=acos(-1.);
   for (freq=2;freq<maxfreq;freq++){
       w=2*pi*freq*df;
       wa=freqweight(freq,df,fmax-10,fmax);
       
       for (ih=0;ih<nh;ih++) daux[ih]=d2[freq][ih];
       matrix_3(rtoep,l,lh,pos,q,nh,nq,w,dh,dq,rtmethod);
       //fprintf(stderr,"rtoep=%e\n",rtoep[0].r);
       powd=rcdot(nh,daux,daux);
       powd/=(nh*nh);
       for (i=0;i<nq;i++) 
	 for (j=0;j<nq;j++)
             if (j>i) ll[i][j]=rtoep[j-i];
	     else ll[i][j]=conjg(rtoep[i-j]);
       //AtimesBm(ll,lh,l,nq,nh,nq);
       for (i=0;i<nq;i++) for (j=0;j<nq;j++) ll0[i][j]=ll[i][j]; 
       Atimesx(maux,lh,daux,nq,nh);

       iter=0;
       while (iter<iter_end){
         iter++;
         fprintf(stderr,"Cholesky....... freq=%d, iter=%d \n",freq,iter);

	 if (iter==1){
           eps2=modgrad(maux,nq,norm,1.,Qp);
           for (i=0;i<nq;i++) ll[i][i]+=eps1/100.*ll[0][0];         
         }
         if (iter!=1){
	   eps2=modgrad(maux2,nq,norm,powd,Qp);
 	   for (i=0;i<nq;i++) for (j=0;j<nq;j++) ll[i][j]=ll0[i][j];
           //for (i=0;i<nq;i++) ll[i][i]+=eps1/100.*ll[0][0];   
           for (i=0;i<nq;i++) ll[i][i]+=Qp[i]+eps;
	 }

         for (i=0;i<nq;i++) 
	   for (j=0;j<nq;j++){ 
	     U[i+1][j+1]=ll[i][j].r;
	     U[nq+i+1][j+1]=ll[i][j].i;
	     U[i+1][nq+j+1]=-ll[i][j].i;
	     U[nq+i+1][nq+j+1]=ll[i][j].r;
           }
        
 
         for (i=0;i<nq;i++){
	   RHS[i+1]=maux[i].r;
           RHS[nq+i+1]=maux[i].i;
         }
       	 choldc(U,nn,W);
	 cholsl(U,nn,W,RHS,x);
 
         for (i=0;i<nq;i++){
	   maux2[i].r=x[i+1];
	   maux2[i].i=x[nq+i+1];           
         } 
           

          
         if (costflag==1){
	   Atimesx(daux2,l,maux2,nh,nq);
	   Jmod=modnorm(maux2,eps2,1,nq);
	   Jdata=misfit(daux2,daux,powd,dh,nh);
           Jtot[iter-1]+=log10(Jdata+Jmod);
	   //fprintf(stderr,"Iter %d, Jd=%e, Jm=%e, J=%e, eps2=%e\n",
	   //	   iter,Jdata,Jmod,(Jdata+Jmod),eps2);
         }
       } //while
    
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
free_vector(x,1,nn);
free_vector(RHS,1,nn);
free_vector(W,1,nn);
free_matrix(V,1,nn,1,nn);
free_matrix(U,1,nn,1,nn);     
free1float(Cd);
free1complex(Qp);
free1complex(rtoep);
free1float(Jtot);
free1complex(daux2);
free1complex(daux);
free1double(diagll);
free1complex(maux2);
free1complex(maux);
free1float(dh);
free2complex(ll0);
free2complex(ll);
free2complex(lh);
free2complex(l);
free2complex(m2);
free2complex(d2);
return;
}







