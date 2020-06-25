#include "su.h"
#include "clibrary.h"
#include <math.h>
/* This is an interface to the WTCGLS method in routine wtcgls.cpp

   This function solves the system of equations 
     (FH WdT Wd FH + WmT Wm ) m = FH WdT Wd d 
   Notice that LH=FH WdT and L= Wd F are computed with matrix_3.
      
   If we assumed noise and model are uncorrelated,
   i.e., data and model covariances matrices Cd and Cm  are diagonal 
   
   Wd is the factorization of the inverse data Covariance matrix
   inv(Cd)=WdT Wd
   
   Wm is the factorization of the inverse model Covariance matrix
   inv(Cm)=WmT Wm

   Wm is useful to increase resolution 
   Wd is useful to taper data at the boundaries and 
      to reduce effect of killed or bad traces

   Wm depends on the model so iteration are necessary
   Wd does not change. 

*/

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void wtcgls0(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj)
{
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,i,j,k,tempint,iter,laststep,flag;
  complex **m2, **d2, czero;
  complex **L, **LH, *R;
  float w,wa,*dh,*dh2,df,*Jtot, Jdata, Jmod, normu;
  complex *d,*u;
  complex *dtemp, *dc, *Wm;
  complex *qcg, *q1cg, *scg, *x1cg, *zcg, *z1cg, *rcg, *Azcg;
  float  *eta, *rho, *gcv;
  float power, Jtotlast, Jtotprev, bb, tempfloat; 
  float *uaux;

  // Foster and Mosher offset function
  float depth=10; // Not used
  float *gFM;  
   
  int costflag=0;
  int freqflag=0;  
  int rtmethod=2;  // Parabolic 


  czero.r=czero.i=0;
  eps1/=100;
   // Note that c allocates memory per column
   // Hence 2 D array with alloc2type requires reversing dimension
   // i.e., first columns, then rows.
   
   if ((d2=alloc2complex(nh,nt))==NULL)
     err("cannot allocate memory for d2\n"); 
   if ((m2=alloc2complex(nq,nt))==NULL)//==> m2(nt x nq) 
     err("cannot allocate memory for m2\n");
   if ((L=alloc2complex(nq,nh))==NULL)   //==> l(nh x nq)
     err("cannot allocate memory for l\n");
   if ((LH=alloc2complex(nh,nq))==NULL)  //==> l(nq xnh)
     err("cannot allocate memory for lh\n");
   if ((dh=alloc1float(nh))==NULL)
     err("cannot allocate memory for dh\n");
   if ((dh2=alloc1float(nh))==NULL)
     err("cannot allocate memory for dh2\n");
   //   if ((u=alloc1complex(nq))==NULL)
   //     err("cannot allocate memory for u\n");
   //    if ((d=alloc1complex(nh))==NULL)
   //      err("cannot allocate memory for d\n");
   if ((Jtot=alloc1float(20))==NULL)
     err("cannot allocate memory for Jtot\n");
   if ((qcg=alloc1complex(nq))==NULL) 
     err("cannot allocate memory for qcg\n");
   if ((q1cg=alloc1complex(nq))==NULL)
     err("cannot allocate memory for q1cg\n");
   if ((scg=alloc1complex(nq))==NULL)
     err("cannot allocate memory for scg\n");
   if ((x1cg=alloc1complex(nq))==NULL)
     err("cannot allocate memory for x1cg\n");
   if ((zcg=alloc1complex(nq))==NULL)
     err("cannot allocate memory for zcg\n");
   if ((z1cg=alloc1complex(nq))==NULL)
     err("cannot allocate memory for z1cg\n");
   if ((rcg=alloc1complex(nh))==NULL)
     err("cannot allocate memory for rcg\n");
   if ((Azcg=alloc1complex(nh))==NULL)
     err("cannot allocate memory for Azcg\n");
   if ((eta=alloc1float(nq))==NULL)
     err("cannot allocate memory for eta\n");
   if ((rho=alloc1float(nq))==NULL)
     err("cannot allocate memory for rho\n");  
   if ((gcv=alloc1float(nq))==NULL)
     err("cannot allocate memory for gcv\n"); 
   if ((Wm=alloc1complex(nq))==NULL)
     err("cannot allocate memory for Wm\n");
   if ((dc=alloc1complex(nh))==NULL)
     err("cannot allocate memory for dc\n");
   if ((R=alloc1complex(nq))==NULL)
     err("cannot allocate memory for R\n");
   if ((uaux=alloc1float(nq))==NULL)
     err("cannot allocate memory for uaux\n");


   for (ih=1;ih<nh;ih++) dh[ih]=pos[ih]-pos[ih-1];
   dh[0]=dh[1];

   if (rtmethod==3){
     if ((gFM=alloc1float(nh))==NULL)
       err("cannot allocate memory for g\n");
     for (i=0;i<nh;i++) gFM[i]=sqrt(pos[i]*pos[i]+depth*depth)-depth;
   }
   //for (ih=0;ih<nh;ih++) dh2[ih]=sqrt(fabs(dh[ih]));
   

   for (i=0;i<=iter_end;i++) Jtot[i]=0.;
   fftgo0(-1,data,d2,nh,nt,dt,&nf0);
   nfreq=nf0/2;
   df=1/(nf0*dt);
   floatprint(df);
   maxfreq=(int) (fmax/df);
   if (freqflag==1) maxfreq=nfreq;
   fprintf(stderr,"maxfreq=%d, dt=%f, df=%f\n",maxfreq,dt,df);
   const double  pi=acos(-1.);
   float test; 
   for (freq=1;freq<maxfreq;freq++){
     w=2*pi*freq*df;
     wa=freqweight(freq,df,fmax-10,fmax);
     for (iq=0;iq<nq;iq++)  {Wm[iq].r=eps1;Wm[iq].i=0.0;}       
     //for (ih=0;ih<nh;ih++) d[ih]=d2[freq][ih];
     matrix_3(L,LH,pos,q,nh,nq,w,dh,Wd,dq,rtmethod,gFM);
     //if (testadj) test=testadj_rad_f(L,LH);
     fprintf(stderr,"freq=%d\n",freq);
     power=rcdot(nh,d2[freq],d2[freq]);

     for (iter=1;iter<=iter_end;iter++){
       for (iq=0;iq<nq;iq++)  {
	 if ((iter!=1)){
       	   normu=abs(m2[freq][iq])+eps2;
	   Wm[iq].r=1/normu+eps1;
	 }
	 else if (0){	   	   
	   normu=abs(uaux[iq])+eps2;
	   Wm[iq].r=1/normu+eps1;
	 }
	 else{ 
	   Wm[iq].r=eps1;
	   Wm[iq].i=0.;
	 }
       }
       wtcgls(d2[freq],L,LH,m2[freq],Wm,Wd,nh,nq,qcg,q1cg,scg,x1cg,zcg,
	      z1cg,rcg,Azcg,eta,rho,gcv,0,step,itercg);
     }
     for (iq=0;iq<nq;iq++) uaux[iq]=abs(m2[freq][iq]);            
     /////////////////////////////////////////////////
     for (iq=0;iq<nq;iq++) m2[freq][iq]/=nq;      
     if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++)  m2[freq][iq]*=wa;
     //for (iq=0;iq<nq;iq++) m2[freq][iq]=u[iq];          
   }
   for (iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++){
       for (iq=0;iq<nq;iq++) m2[freq][iq]=czero;
   }     
   fprintf(stderr,"w=%f\n",w);
   exit=fftback0(1,model,m2,nq,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");

if (rtmethod==3) free1float(gFM);
free1float(uaux);
free1complex(R);
free1complex(dc);
free1complex(Wm);
free1float(gcv);
free1float(rho);
free1float(eta);
free1complex(Azcg);
free1complex(rcg);
free1complex(z1cg);
free1complex(zcg);
free1complex(x1cg);
free1complex(scg);
free1complex(q1cg);
free1complex(qcg);
free1float(Jtot);
//free1complex(d);
//free1complex(u);
free1float(dh2);
free1float(dh);
free2complex(LH);
free2complex(L);
free2complex(m2);
free2complex(d2);
return;
}
























