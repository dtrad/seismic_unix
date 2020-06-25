#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"
#include "/usr/local/cwp/src/Complex/include/Complex.h"
#include "/home/dtrad/radon/clibrary/clibrary.h"
#include "Dcomplex.h"
#include <math.h>
dcomplex sin(dcomplex);
#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void cholover(float *pos, float **d, float dt,float **m, float *q, float dq,
	     float eps1, float eps2, float eps, float fmax)
{
   int ih, iq, freq, maxfreq, nfreq, exit, nf0, i,j, iter;
   complex **m2, **d2,czero, tempcomp;
   complex **l, **lh, **ll, *Qp, **ll0;
   float w, wa, *dh,df,*Jtot,Jmod,Jdata,powd,*Cd;
   double *diagll;
   complex *maux, *maux2, *daux, *daux2;
   complex *rtoep;
   extern int nt, nh, nq, method, iter_end, rtmethod, norm,freqflag,costflag;
   FILE *myfile3;                
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

   for (ih=0;ih<nh;ih++) Cd[ih]=1.;
   for (ih=1;ih<nh;ih++) dh[ih]=pos[ih]-pos[ih-1];
   //for (ih=1;ih<nh;ih++) dh[ih]=1;
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
       
       matrix_3(rtoep,l,lh,pos,q,nh,nq,w,dh,dq,rtmethod,Cd);
       powd=rcdot(nh,daux,daux);
       for (i=0;i<nq;i++) 
	 for (j=0;j<nq;j++)
             if (j>i) ll[i][j]=rtoep[j-i];
	     else ll[i][j]=conjg(rtoep[i-j]);
       //AtimesBm(ll,lh,l,nq,nh,nq);
       fprintf(stderr,"ll[0][0].r=%e\n",ll[0][0].r);
       fprintf(stderr,"ll[0][0].i=%e\n",ll[0][0].i);
       for (i=0;i<nq;i++) for (j=0;j<nq;j++) ll0[i][j]=ll[i][j]; 
       Atimesx(maux,lh,daux,nq,nh);

       iter=0;
       while (iter<iter_end){
         iter++;
	 if (iter==1){
           if (norm==10) eps2=modgrad(maux,nq,norm,1.,Qp);
           else eps2=1.;
           fprintf(stderr,"freq=%d, iter=%d \n",freq,iter); 
           for (i=0;i<nq;i++) ll[i][i]+=eps1/100.*ll[0][0];         
         }
         if (iter!=1){
	   if (norm==10) eps2=modgrad(maux2,nq,norm,1.,Qp,eps2);
           else eps2=1.;
 	   for (i=0;i<nq;i++) for (j=0;j<nq;j++) ll[i][j]=ll0[i][j];
           if (norm==10) for (i=0;i<nq;i++) ll[i][i]+=Qp[i]+eps;
           if (norm==1) for (i=0;i<nq;i++) ll[i][i]+=1./(abs(maux2[i])+1e-9)+eps;
	 }  
    	 //displayA(Qp,nq);
       	 choldc(ll,nq,diagll);
	 cholsl(ll,nq,diagll,maux,maux2);
          
         if (costflag==1){
	   Atimesx(daux2,l,maux2,nh,nq);
	   Jmod=modnorm(maux2,eps2,1,nq);
	   Jdata=misfit(daux2,daux,powd,dh,nh);
           Jtot[iter-1]+=log10(Jdata+Jmod);
	   fprintf(stderr,"Iter %d, Jd=%e, Jm=%e, J=%e, eps2=%e\n",
	   	   iter,Jdata,Jmod,(Jdata+Jmod),eps2);
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







