#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "cwp.h"
#include "clibrary.h"
#include <math.h>

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

void radoncg_levinson_alias(float *pos, int nh, float **data, float *t,int nt,float dt,float **model, float *q, int nq, float dq,float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil, int rtmethod)
{
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,k,iter,flag,nf2;
  register int i;
  complex **m2, **d2, czero;
  complex **L, **LH, *R, *RC, *fc, *gc;
  float w,wa,*dh,df;
  complex *u;
  int freqflag;
  complex *g1;    // residuals r0 
  complex *g2;    // Conjugate gradient
  complex *Qp;    // Regularization
  complex *gtemp, *gtemp2;  // Temporal vectors required to compute Ap
  complex *Ag;  // Matrix times conjugate gradient
  complex *rhs; // Right hand side LH*d
  float alphanum, alphaden, alpha, betanum, betaden, beta;
  float power, resid, bb;
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
  if ((g1=alloc1complex(nq))==NULL)
    err("cannot allocate memory for g1\n");
  if ((g2=alloc1complex(nq))==NULL)
    err("cannot allocate memory for g2\n"); 
  if ((Qp=alloc1complex(nq))==NULL)
    err("cannot allocate memory for Qp\n");
  if ((gtemp=alloc1complex(nq))==NULL)
    err("cannot allocate memory for gtemp\n");
  if ((rhs=alloc1complex(nq))==NULL)
    err("cannot allocate memory for rhs\n");
  if ((gtemp2=alloc1complex(nq))==NULL)
    err("cannot allocate memory for gtemp2\n");
  if ((Ag=alloc1complex(nq))==NULL)
    err("cannot allocate memory for Ag\n");

  complex *uold=ealloc1complex(nq); // Save old model to apply antialias approach
  float *mtemp=ealloc1float(nq);
  
  fprintf(stderr,"quantil=%f\n",quantil);

  for (ih=1;ih<nh;ih++) dh[ih]=pos[ih]-pos[ih-1]; dh[0]=dh[1]; 
  
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
    matrix_3(R,L,LH,pos,q,nh,nq,w,dh,dq,rtmethod,Wd);
    fprintf(stderr,"freq=%d\n",freq);
    power=rcdot(nh,d2[freq],d2[freq]);
    //power/=(nh*nh);
    //power=1;
    // Construct the circulant matrix//////////////////////////////////////
    for (i=0; i<nq;i++) RC[i]=conjg(R[i]);       
    for (i=nq; i<nf2;i++) RC[i]=czero; // pad with zeros the autocorrelation
    for (i=1;i<=nq-1;i++) RC[nf2-i]=conjg(RC[i]);      //use the DFT
    pfacc(-1,nf2,RC);
    ///////////////////////////////////////////////////////////////////////
    Atimesx(rhs,LH,d2[freq],nq,nh);
    bb=rcdot(nq,rhs,rhs);
    iter=1;    
    

    
    while (iter<=iter_end){
      xequaly(g1,rhs,nq);  // resid g1 = rhs  
      xequaly(g2,rhs,nq);  // Conj Grad g2 = rhs
      // Antialias approach from Phillipe Herrmann
      // Use previous model to compute Qp
      if (freq>freq0) iter=iter_end;

      // Compute Qp //////////////////////////////////////////////////////
      //power=rcdot(nq,g1,g1);
      // Setting power to larger values increase the non linearity
      // using smalll values make more least squares results. It is equivalent
      // to the Lagrange multiplier.      
      power=1;  
      for (i=0;i<nq;i++) mtemp[i]=abs(u[i]);
      float qmodel=quest(quantil,nq,mtemp);     
      if (iter==1) for (i=0;i<nq;i++){
	Qp[i].r=eps1;
	Qp[i].i=0.0;
      }
      else for (i=0;i<nq;i++){
	ctemp=sqrt((conjg(u[i])*u[i])/(qmodel*qmodel));
	Qp[i].r=power/(ctemp.r+eps2)+1e-7; 
	Qp[i].i=0;                
      }
      ////////////////////////////////////////////////////////////////////
      for (i=0;i<nq;i++) u[i].r=u[i].i=0.0;    

      resid=rcdot(nq,g1,g1);

      k=0;
      while (((sqrt(resid)>(eps*bb))&&(k<itercg)&&(0.9*resold>resid))||(k<2)){
	k++;
	
	///// Compute alpha
	alphanum=rcdot(nq,g1,g1);
	
	if (option==0) circ_mult(nq,RC,g2,gtemp,nf2,fc,gc); 
	xtimesy(gtemp2,Qp,g2,nq);
	xplusy(Ag,gtemp,gtemp2,nq);

	alphaden=rcdot(nq,g2,Ag);
	if (0)
	for (i=0;i<nq;i++) 
	  fprintf(stderr,"gtemp[%d]=%f,gtemp2[%d]=%f\n",i,abs(gtemp[i]),i,abs(gtemp2[i]));
	if (alphaden<eps) break;
	alpha=step*alphanum/alphaden;
	//fprintf(stderr,"alpha=%f\n",alpha);
	for(i=0;i<nq;i++){
	  u[i]=u[i]+alpha*g2[i];
	  g1[i]=g1[i]-alpha*Ag[i];  
	}

	///// Compute Beta
	betanum=rcdot(nq,g1,g1);
	betaden=alphanum;              
	if (betaden<eps) break;
	beta=betanum/betaden;	      
	//fprintf(stderr,"beta=%f\n",beta);
	for(i=0;i<nq;i++) g2[i]=g1[i]+beta*g2[i];
	  
	resold=resid;
	resid=betanum;
	
	if ((iter==1)&&(k==15)) break;
	
      } 
      
      iter++;
      intprint(k);
      //fprintf(stderr,"resid=%f\n",resid);

    }  //loop for niter        
    /////////////////////////////////////////////////
    
    if ((wa<1)&&(wa>0)) for (iq=0;iq<nq;iq++)  m2[freq][iq]=u[iq]*wa;
    else for (iq=0;iq<nq;iq++) m2[freq][iq]=u[iq];
  }
  for (freq=0;freq<freq0;freq++) for (iq=0;iq<nq;iq++)  m2[freq][iq]=czero;  
  for (freq=maxfreq;freq<nfreq;freq++) for (iq=0;iq<nq;iq++) m2[freq][iq]=czero;

  fprintf(stderr,"w=%f\n",w);
  exit=fftback0(1,model,m2,nq,nt,dt,nf0);
  if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");

  free1complex(uold);
  free1float(mtemp);
  free1complex(Ag);
  free1complex(gtemp2);
  free1complex(rhs);
  free1complex(gtemp);
  free1complex(Qp);
  free1complex(g2);
  free1complex(g1);
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
  free2complex(m2);
  free2complex(d2);
  return;
}

/*
      ///////////////////////////////////////////////////////////////////
      // The following options are only to test multiplication with FFT
      // In general option=0
      //////////////////////////////////////////////////////////////////
      if (option==2){
	circ_mult(nq,RC,u,gtemp,nf2,fc,gc);
	for (i=0;i<10;i++){
	  fprintf(stderr,"Circ_matrix gtemp[%d].r=%e\n",i,gtemp[i].r);
	  fprintf(stderr,"Circ_matrix gtemp[%d].i=%e\n",i,gtemp[i].i);
	}
      }

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
      //////////////////////////////////////////////////////////////////////

*/


















