#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrary.h"
#include <math.h>

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */
void lsqr0(float *pos, float **data, float dt,float **model, float *q, float dq,float eps1, float eps2, float eps, float fmax)
{
  /*
    % LSQR Solution of least squares problems by 
    % Lanczos bidiagonalization with/without reorthogonalization.
    % Reference: C. C. Paige & M. A. Saunders, "LSQR: an algorithm for
    % sparse linear equations and sparse least squares", ACM Trans.
    % Math. Software 8 (1982), 43-71.  */
 
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,i,j,k,tempint,iter,laststep,flag;
  complex **m2, **d2, czero;
  complex **L, **LH;
  float w,wa,*dh,df,*Cd, normu;
  complex *d,*u;
  extern int nt,nh,nq,method,iter_end,rtmethod,norm,itercg,freqflag,costflag;
  extern int reorth;
  complex *dtemp, *dc, *Qp;
  complex **U,**B,**V,*u2,*v,*w2,*q2,*p,*r,*temp2;
  float *eta,*rho,*gcv;  
  float power, bb, tempfloat;
  czero.r=czero.i=0;
  int maxn;
  maxn= (nq > nh) ? nq : nh;  

   // Note that c allocates memory per column
   // Hence 2 D array with alloc2type requires reversing dimension
   // i.e., first columns, then rows.
   if ((d2=alloc2complex(nh,(nt+100)))==NULL)
     err("cannot allocate memory for d2\n"); 
   if ((m2=alloc2complex(nq,(nt+100)))==NULL)//==> m2(nt x nq) 
     err("cannot allocate memory for m2\n");
   if ((L=alloc2complex(nq,nh))==NULL)   //==> l(nh x nq)
     err("cannot allocate memory for l\n");
   if ((LH=alloc2complex(nh,nq))==NULL)  //==> l(nq xnh)
     err("cannot allocate memory for lh\n");
   if ((dh=alloc1float(nh))==NULL)
     err("cannot allocate memory for dh\n"); 
   if ((u=alloc1complex(nq))==NULL)
     err("cannot allocate memory for u\n");
   if ((d=alloc1complex(nh))==NULL)
     err("cannot allocate memory for d\n");
   if ((Cd=alloc1float(nh))==NULL)
     err("cannot allocate memory for Cd\n");
  if ((U=alloc2complex(itercg,nh))==NULL)   //==> l(ny x itercg)
    err("cannot allocate memory for U\n");
  if ((V=alloc2complex(itercg,nq))==NULL)  //==> l(nq xnh)
    err("cannot allocate memory for V\n");
  if ((B=alloc2complex(2,itercg))==NULL)
     err("cannot allocate memory for B\n");
  if ((Qp=alloc1complex(nq))==NULL) 
    err("cannot allocate memory for Qp\n");
  if ((q2=alloc1complex(nq))==NULL) 
    err("cannot allocate memory for q2\n");
  if ((v=alloc1complex(nq))==NULL)
    err("cannot allocate memory for v\n");
  if ((u2=alloc1complex(nh))==NULL)
    err("cannot allocate memory for u\n");
  if ((w2=alloc1complex(nq))==NULL)
    err("cannot allocate memory for w2\n");
  if ((r=alloc1complex(nq))==NULL)
    err("cannot allocate memory for r\n");
  if ((p=alloc1complex(nh))==NULL)
    err("cannot allocate memory for p\n");
  if ((temp2=alloc1complex(maxn))==NULL)
    err("cannot allocate memory for temp2\n");
  if ((eta=alloc1float(itercg))==NULL)
    err("cannot allocate memory for eta\n");
  if ((rho=alloc1float(itercg))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((gcv=alloc1float(itercg))==NULL)
    err("cannot allocate memory for gcv\n");
  //for (i=0;i<nq;i++) fprintf(stderr,"q[%d]=%e\n",i,q[i]); 
   for (ih=0;ih<nh;ih++) Cd[ih]=1.;
   //for (ih=1;ih<nh;ih++) dh[ih]=pos[ih]-pos[ih-1];
   for (ih=1;ih<nh;ih++) dh[ih]=1.;
   dh[0]=dh[1]; 

   //fprintf(stderr,"lsqr0:nh=%d,nt=%d,nq=%d,itercg=%d\n",nh,nt,nq,itercg);
   fftgo(-1,data,d2,nh,nt,dt,&nf0);

   nfreq=(nf0-2)/2+1;
   df=1/(nfreq*dt*2);
   floatprint(df);
   maxfreq=(int) ((fmax/df)-1);
   if (freqflag==1) maxfreq=nfreq;
   fprintf(stderr,"maxfreq=%d,dt=%f,df=%f\n",maxfreq,dt,df);
   const double  pi=acos(-1.);
   if (reorth) fprintf(stderr,"Lanczos with reorthonalization\n");
   else fprintf(stderr,"Lanczos without reorthonalization\n"); 
   for (freq=1;freq<maxfreq;freq++){
     w=2*pi*freq*df;  
     wa=freqweight(freq,df,fmax-10,fmax);
     for (iq=0;iq<nq;iq++)  {Qp[iq].r=eps1;Qp[iq].i=0.0;}    
     for (ih=0;ih<nh;ih++) d[ih]=d2[freq][ih];

  
     matrix_3(L,LH,pos,q,nh,nq,w,dh,dq,rtmethod,flag);     
    
     power=rcdot(nh,d,d);
     for (iter=1;iter<=iter_end;iter++){
       for (iq=0;iq<nq;iq++)  {
	 normu=abs(u[iq]);
	 if ((normu > 1.e-5)&&(iter!=1)) Qp[iq].r=1/normu;
	 else Qp[iq].r=eps1;
	 Qp[iq].i=0.0;
       }
       lsqrc(L,LH,u,d,eps,eps2,reorth,U,B,V,Qp,u2,v,w2,q2,p,r,temp2,eta,rho,gcv); 
     }
     //displayA(u,nq);       
     /////////////////////////////////////////////////
     for (iq=0;iq<nq;iq++) u[iq]/=nq;      
     if ((wa<1)&&(wa>0))
	for (iq=0;iq<nq;iq++)  u[iq]*=wa;
     for (iq=0;iq<nq;iq++)
       m2[freq][iq]=u[iq];     
   }
   for (iq=0;iq<nq;iq++)     m2[0][iq]=czero;  //dc can not be recovered  
   for (freq=maxfreq;freq<nfreq;freq++){
       for (iq=0;iq<nq;iq++)
          m2[freq][iq]=czero;
   }     
   fprintf(stderr,"w=%f\n",w);
   exit=fftback(1,model,m2,nq,nt,dt,nf0);
   if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
   fprintf(stderr,"After fftback -1 nh=%d, nt=%d \n",nh,nt);

  
  free1float(gcv);
  free1float(rho);
  free1float(eta);
  free1complex(temp2);
  free1complex(p);
  free1complex(r);
  free1complex(w2);
  free1complex(u2);
  free1complex(v);
  free1complex(q2); 
  free1complex(Qp);
  free2complex(B);
  free2complex(V); 
  free2complex(U);
  free1float(Cd);
  free1complex(d);
  free1complex(u);
  free1float(dh);
  free2complex(LH);
  free2complex(L);
  free2complex(m2);
  free2complex(d2);
  fprintf(stderr,"After lsqr0 all arrays were freed\n");         
  return;
}














