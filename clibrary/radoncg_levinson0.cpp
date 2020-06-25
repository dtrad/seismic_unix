#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "cwp.h"
#include "clibrary.h"
#include <math.h>

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

void radoncg_levinson(float *pos, int nh, float **data, float *t,int nt,float dt,float **model, float *q, int nq, float dq,float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil)
{
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,k,iter,flag,nf2;
  register int i;
  int flageps2=0;
  complex **m2, **d2, czero;
  complex **L, **LH, *R, *RC, *fc, *gc, **LL;
  float w,wa,*dh,df,*Jtot;
  complex *d,*dw,*u;
  int rtmethod=2;
  int freqflag,costflag;
  complex *uold, *dtemp, *dc;;
  complex *g1, *g1old, *g2,  *Qp;
  complex *gtemp, *gtemp2, *gtemp4, *gtemp5;
  float alfanum, alfaden, alfa, betanum, betaden, beta, Rorig;
  float power, resid, Jtotlast, Jtotprev, bb;
  complex ctemp;
  float resold;
  
  // option=0 ------> circulant matrix-vector multplication
  // option=1 ------> Conventional matrix vector multiplication
  // option=2 ------> Both methods: prints out the numbers to check.
  int option=0;
  freqflag=0;
  czero.r=czero.i=0;
  nf2=npfa((int) 2*nq);
  intprint(nf2);
  
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
  if ((LL=alloc2complex(nq,nq))==NULL)  //==> l(nq xnh)
    err("cannot allocate memory for LL\n");
  if ((R=alloc1complex(nq))==NULL)  //==> l(nq xnh)
    err("cannot allocate memory for R\n");
  if ((RC=alloc1complex(nf2))==NULL)  //==> l(nq xnh)
    err("cannot allocate memory for RC\n");
  if ((fc=alloc1complex(nf2))==NULL)  //==> l(nq xnh)
    err("cannot allocate memory for fc\n");
  if ((gc=alloc1complex(nf2))==NULL)  //==> l(nq xnh)
    err("cannot allocate memory for gc\n");
  if ((dh=alloc1float(nh))==NULL)
    err("cannot allocate memory for dh\n"); 
  if ((u=alloc1complex(nq))==NULL)
    err("cannot allocate memory for u\n");
  if ((d=alloc1complex(nh))==NULL)
    err("cannot allocate memory for d\n");
  if ((dw=alloc1complex(nh))==NULL)
    err("cannot allocate memory for dw\n");
  if ((Jtot=alloc1float(20))==NULL)
    err("cannot allocate memory for Jtot\n");
  if ((uold=alloc1complex(nq))==NULL) 
    err("cannot allocate memory for uold\n");
  if ((g1=alloc1complex(nq))==NULL)
    err("cannot allocate memory for g1\n");
  if ((g1old=alloc1complex(nq))==NULL)
    err("cannot allocate memory for g1old\n");
  if ((g2=alloc1complex(nq))==NULL)
    err("cannot allocate memory for g2\n"); 
  if ((Qp=alloc1complex(nq))==NULL)
    err("cannot allocate memory for Qp\n");
  if ((gtemp=alloc1complex(nq))==NULL)
    err("cannot allocate memory for gtemp\n");
  if ((gtemp2=alloc1complex(nq))==NULL)
    err("cannot allocate memory for gtemp2\n");
  if ((gtemp4=alloc1complex(nq))==NULL)
    err("cannot allocate memory for gtemp4\n");
  if ((gtemp5=alloc1complex(nq))==NULL)
    err("cannot allocate memory for gtemp5\n");
  if ((dc=alloc1complex(nh))==NULL)
    err("cannot allocate memory for dc\n");
  if ((dtemp=alloc1complex(nh))==NULL)
    err("cannot allocate memory for dtemp\n");

  float *mtemp=ealloc1float(nq);
  fprintf(stderr,"quantil=%f\n",quantil);
  for (ih=1;ih<nh;ih++)
    dh[ih]=pos[ih]-pos[ih-1];
  dh[0]=dh[1]; 
  
  fftgo0(-1,data,d2,nh,nt,dt,&nf0);
  
  nfreq=(nf0-2)/2+1;
  df=1/(nfreq*dt*2);
  floatprint(df);
  maxfreq=(int) ((fmax/df)-1);
  if (freqflag==1) maxfreq=nfreq;
  const double  pi=acos(-1.);
  int freq0=5;


  for (freq=freq0;freq<maxfreq;freq++){
    w=2*pi*freq*df;
    wa=freqweight(freq,df,fmax-10,fmax);
    
    for (ih=0;ih<nh;ih++) d[ih]=d2[freq][ih];
    //for (ih=0;ih<nh;ih++) dh[ih]/=nfreq;          
    matrix_3(R,L,LH,pos,q,nh,nq,w,dh,dq,rtmethod,Wd);
    fprintf(stderr,"freq=%d\n",freq);
    power=rcdot(nh,d,d);
    power/=(nh*nh);
    //power=1; 
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //     Unconstrained gradient step
    
    //-----Apply weights LH
    //       dtemp=LH*d
    Rorig=R[0].r;
    R[0].r*=(1+eps1);
    
    for (i=0; i<nq;i++) RC[i]=conjg(R[i]);       
    for (i=nq; i<nf2;i++) RC[i]=czero; // pad with zeros the autocorrelation
    for (i=1;i<=nq-1;i++) RC[nf2-i]=conjg(RC[i]);      //use the DFT
    
    pfacc(-1,nf2,RC);
    //for (i=0;i<nf2;i++) RC[i]/=nf2;
    //pfacc(1,nf2,RC);
    //displayA(RC,nq);
    for (i=0;i<nh;i++) dw[i]=d[i]*Wd[i];  
    Atimesx(gtemp2,LH,d,nq,nh);
    
    bb=rcdot(nq,gtemp2,gtemp2);
    xequaly(u,gtemp2,nq);  // Initial guess
    
    //-   Minimum for u=(eps,eps)
    //for(i=0;i<nq;i++) if(abs(u[i]) < sqrt(eps)) u[i].r=u[i].i=eps;
    if (iter!=1) R[0].r=Rorig;
    iter=1; Jtotlast=0; Jtotprev=10; flag=0; //This flag is set to 1 when J increases
    //R[0].r=Rorig+eps;
    if (norm==1) norm=10; // Norm=1 is not workin well 
    while ((iter<=iter_end)&&(flag==0)){
      power=rcdot(nq,u,u);
      for (i=0;i<nq;i++) mtemp[i]=abs(u[i]);
      float qmodel=quest(quantil,nq,mtemp);     
      for (i=0;i<nq;i++){
	ctemp=(conjg(u[i])*u[i])/(qmodel*qmodel);
	Qp[i].r=power/(ctemp.r+eps2); 
	Qp[i].i=0;                
      }
      //for (iq=1;iq<=5;iq++) Qp[(nq-6+iq)]=iq*iq*Qp[(nq-6+iq)];
      //if (norm==1) {
	//eps2=1e-3;
        //eps1=1e-2;
      //for (i=0;i<nq;i++) Qp[i]=1./(abs(u[i])+1e-9)+eps;
      //}         
      if (option==0) circ_mult(nq,RC,u,gtemp,nf2,fc,gc);
      
      if (option==2){
	circ_mult(nq,RC,u,gtemp,nf2,fc,gc);
	for (i=0;i<10;i++){
	  fprintf(stderr,"Circ_matrix gtemp[%d].r=%e\n",i,gtemp[i].r);
	  fprintf(stderr,"Circ_matrix gtemp[%d].i=%e\n",i,gtemp[i].i);
	}
      }

      //fprintf(stderr,"R[0]=%e\n",R[0].r);
      if (option==1){
	AtimesBm(LL,LH,L,nq,nh,nq);
	Atimesx(gtemp,LL,u,nq,nq);
      }
     
      if (option==2){
	AtimesBm(LL,LH,L,nq,nh,nq);
	Atimesx(gtemp,LL,u,nq,nq);
	for (i=0;i<10;i++){
	  fprintf(stderr,"LL*U gtemp[%d].r=%e\n",i,gtemp[i].r);
	  fprintf(stderr,"LL*U gtemp[%d].i=%e\n",i,gtemp[i].i);
	} 
      }

      xtimesy(gtemp4,Qp,u,nq);
      xplusy(gtemp,gtemp,gtemp4,nq);     
      xminusy(g1,gtemp2,gtemp,nq);
      //fprintf(stderr,"iter=%d,eps2=%f\n",iter,eps2);
      xequaly(uold,u,nq);
      
      //- Gradient steps
      //- NQ loop
      resid=rcdot(nq,g1,g1);
      k=0; //power=1; 
      //while ((k!=nq)&&(k<itercg)){
      while (((k!=nq)&&(sqrt(resid)>(eps*bb))&&(k<itercg)&&(0.9*resold>resid))||(k<2)){
	k++;
	///// Compute Beta
	if (k==1) xequaly(g2,g1,nq);
	else{ 
	  betanum=resid;
	  betaden=alfanum;              
	  if (betaden<eps) {;
	  //fprintf(stderr,"betanum=%e,betaden=%e,k=%d,eps=%e\n",
	  //  betanum,betaden,k,eps); 
	  break;
	  }
	  beta=betanum/betaden;	      
	  for(i=0;i<nq;i++) g2[i]=g1[i]+beta*g2[i];
	}
	///// Compute alfa
	alfanum=rcdot(nq,g1,g1);
	
	if (option==0) circ_mult(nq,RC,g2,gtemp,nf2,fc,gc); 
	else if (option==1) Atimesx(gtemp,LL,g2,nq,nq); 
	else if (option==2) Atimesx(gtemp,LL,g2,nq,nq); 

	xtimesy(gtemp4,Qp,g2,nq);
	xplusy(gtemp5,gtemp,gtemp4,nq);
	alfaden=rcdot(nq,g2,gtemp5);
	//if (alfaden<0) err("alfaden=%e\n",alfaden);
	if (alfaden<eps) {
	  //fprintf(stderr,"alfanum=%e,alfaden=%e,k=%d,eps=%e\n",
	  //	  alfanum,alfaden,k,eps);
	  break;
	} 
	alfa=step*alfanum/alfaden;
	for(i=0;i<nq;i++){
	  u[i]=u[i]+alfa*g2[i];
	  g1[i]=g1[i]-alfa*gtemp5[i];  
	}
	resold=resid;
	resid=rcdot(nq,g1,g1);
	if ((iter==1)&&(k==15)) break;
      } 
      /*    //!loop for k
	    Atimesx(dc,L,u,nh,nq);
	    Jtot[iter]=costfunc(d,dc,u,nh,nq,eps2,norm);
	    fprintf(stderr,"Jiter=%f,%d,%d,%f\n",Jtot[iter],iter,k,freq);
	    Jtotlast=Jtot[iter];
	    if (iter==1) Jtotprev=Jtot[iter]+10; else Jtotprev=Jtot[iter-1];
	    for(i=iter;i<iter_end;i++){
	    Jtot[i]=Jtot[iter] ;//In case of stop iterating keep the last J
	    
	    if ((Jtotlast>(1.0*Jtotprev))&&(iter!=1)){
	    flag=1;
	    for(i=0;i<nq;i++)
	    u[i]=uold[i];
	    }
	    }*/
      iter++;
      intprint(k);
     
      //for (iq=0;iq<5;iq++) {
      //  float scale=1-exp(-(iq+.3));
      //  u[iq]*=scale;
      //  u[nq-iq-1]*=scale;
      // }
    }  //loop for niter        
    /////////////////////////////////////////////////
    
    if ((wa<1)&&(wa>0))
      for (iq=0;iq<nq;iq++)  u[iq]*=wa;
    for (iq=0;iq<nq;iq++)
      m2[freq][iq]=u[iq];
    
    
  }
  if (1){
    for (freq=0;freq<freq0;freq++) 
      for (iq=0;iq<nq;iq++)  
	m2[freq][iq]=czero;  //dc can not be recovered  
    for (freq=maxfreq;freq<nfreq;freq++){
      for (iq=0;iq<nq;iq++)
	m2[freq][iq]=czero;
    }     
  }
  fprintf(stderr,"w=%f\n",w);
  exit=fftback0(1,model,m2,nq,nt,dt,nf0);
  if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
  
  free1float(mtemp);
  free1complex(dtemp);
  free1complex(dc);
  free1complex(gtemp5);
  free1complex(gtemp4);
  free1complex(gtemp2);
  free1complex(gtemp);
  free1complex(Qp);
  free1complex(g2);
  free1complex(g1old);
  free1complex(g1);
  free1complex(uold);      
  
  free1float(Jtot);
  free1complex(dw);
  free1complex(d);
  free1complex(u);
  free1float(dh);
  // circ_mult
  free1complex(gc);
  free1complex(fc);
  free1complex(RC);
  free1complex(R);
  /////////
  free2complex(LH);
  free2complex(L);
  free2complex(LL);
  free2complex(m2);
  free2complex(d2);
  return;
}























