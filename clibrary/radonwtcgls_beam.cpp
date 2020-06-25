#include "su.h"
#include "clibrary.h"
#include <math.h>

void wtcgls(complex *b,complex **L, complex *x,complex *Wm,
	  float *Wd,int nh,int nq,complex *q,complex *q1,complex *s,
	  complex *x1,complex *z,complex *z1,complex *r,complex *Az,
	  float *eta,float *rho, float *gcv, float tol, float step, int itercg);

void matrix(complex **L, float *h, float *q,int nh,int nq,float w,float dq,int rtmethod);


void radon_param_beam(float fmax1, float fmax2, float *h, int nh,  float *q, float nq, 
		      float qmin1, float qmin2, float *q1, 
		      float *q2, float *ph1, float *ph2, int nq1, int nq2, float depth1, 
		      float depth2, int  rtmethod1, int rtmethod2, float factor1, 
		      float factor2, float *pdq1, float *pdq2);

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
void radonwtcgls_beam(float *h, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, float quantil, int nq1, int nq2, int rtmethod1, int rtmethod2, float depth1, float depth2, float qmin1, float qmin2, float factor1, float factor2)
{
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,i,j,k,tempint,iter,laststep,flag;
  complex **m2;
  complex **d2;
  complex czero;
  complex **L;
  complex **L1;
  complex **L2;

  float w,wa,*dh,*dh2,df,*Jtot, Jdata, Jmod, normu;
  complex *d,*u;
  complex *dtemp, *dc, *Wm;
  complex *qcg, *q1cg, *scg, *x1cg, *zcg, *z1cg, *rcg, *Azcg;
  float  *eta, *rho, *gcv;
  float power, Jtotlast, Jtotprev, bb, tempfloat; 
  float *uaux;
  float dq1;
  float dq2;
  float *q1;
  float *q2;
  int costflag=0;
  int freqflag=0;  
  float *ph1;
  float *ph2;
  
  czero.r=czero.i=0;
  //eps1/=100.;
  
  d2=ealloc2complex(nh,nt);
  
  m2=ealloc2complex(nq,nt);//==> m2(nt x nq) 
  
  L=ealloc2complex(nq,nh);   //==> l(nh x nq)
  
  L1=ealloc2complex(nq1,nh);   //==> l(nq xnh)
  
  L2=ealloc2complex(nq2,nh);   //==> l(nq xnh)

  q1=ealloc1float(nq1);  

  q2=ealloc1float(nq2);

  dh=ealloc1float(nh);
  
  dh2=ealloc1float(nh);
  
  Jtot=ealloc1float(20);
  
  qcg=ealloc1complex(nq); 

  q1cg=ealloc1complex(nq);

  scg=ealloc1complex(nq);
  
  x1cg=ealloc1complex(nq);
  
  zcg=ealloc1complex(nq);

  z1cg=ealloc1complex(nq);
  
  rcg=ealloc1complex(nh);
  
  Azcg=ealloc1complex(nh);
  
  eta=ealloc1float(nq);
  
  rho=ealloc1float(nq);

  gcv=ealloc1float(nq);

  Wm=ealloc1complex(nq);
  
  dc=ealloc1complex(nh);
  
  uaux=ealloc1float(nq);
  
  ph1=ealloc1float(nh);

  ph2=ealloc1float(nh);

  float *Wm2=ealloc1float(nq);

  for (iq=0;iq<nq1;iq++) Wm2[iq]=1;
  for (iq=0;iq<nq2;iq++) Wm2[nq1+iq]=1;
  

  ////////////////////////////////////////////////////

  radon_param_beam(fmax,fmax,h,nh,q, nq, qmin1,qmin2,q1,q2,ph1,ph2,nq1,nq2,depth1,depth2,
		   rtmethod1, rtmethod2,factor1,factor2,&dq1,&dq2);

  ///////////////////////////////////////////////////
  
  float *mtemp=ealloc1float(nq);
  
  for (ih=1;ih<nh;ih++) dh[ih]=h[ih]-h[ih-1];
  dh[0]=dh[1];

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
     for (iq=0;iq<nq;iq++){
       Wm[iq].r=eps1*Wm2[iq];
       Wm[iq].i=0.0;
     }    
     
     if (rtmethod1==3) matrix(L1,ph1,q1,nh,nq1,w,dq1,rtmethod1);
     else  matrix(L1,h,q1,nh,nq1,w,dq1,rtmethod1);
     if (rtmethod2==3) matrix(L2,ph2,q2,nh,nq2,w,dq2,rtmethod2);
     else matrix(L2,h,q2,nh,nq2,w,dq2,rtmethod2);

     for (ih=0;ih<nh;ih++){
       for (iq=0;iq<nq1;iq++)
	 L[ih][iq]=L1[ih][iq];
       for (iq=0;iq<nq2;iq++) 
	 L[ih][iq+nq1]=L2[ih][iq];
     }	 
     
     //if (testadj) test=testadj_rad_f(L,LH);
     fprintf(stderr,"freq=%d\n",freq);
     power=rcdot(nh,d2[freq],d2[freq]);
     
     for (iter=1;iter<=iter_end;iter++){

       if (iter>1) for (iq=0;iq<nq;iq++) mtemp[iq]=abs(m2[freq][iq]); 
    
       for (iq=0;iq<nq;iq++)  {
	 if ((iter!=1)){
	   float qmodel=quest(quantil,nq,mtemp);
	   //qmodel=1;
       	   normu=abs(m2[freq][iq])/qmodel+eps2;
	   Wm[iq].r=1/normu+eps1;
	   Wm[iq].r*=Wm2[iq];
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
       wtcgls(d2[freq],L,m2[freq],Wm,Wd,nh,nq,qcg,q1cg,scg,x1cg,zcg,
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

   free1float(Wm2);
   free1float(q1);
   free1float(q2);
   free1float(ph1);
   free1float(ph2);
   free1float(mtemp);
   free1float(uaux);
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
   free1float(dh2);
   free1float(dh);
   free2complex(L2);
   free2complex(L1);
   free2complex(L);
   free2complex(m2);
   free2complex(d2);
   return;
}

void radonwtcgls_beam2(float *h, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, float quantil, int nq1, int nq2, int rtmethod1, int rtmethod2, float depth1, float depth2, float qmin1, float qmin2, float factor1, float factor2, float fmax1, float fmax2)
{
  int ih,iq,freq,maxfreq,nfreq,exit,nf0,i,j,k,tempint,iter,laststep,flag;
  complex **m2;
  complex **d2;
  complex czero;
  complex **L;
  complex **L1;
  complex **L2;

  float w,wa,*dh,*dh2,df,*Jtot, Jdata, Jmod, normu;
  complex *d,*u;
  complex *dtemp, *dc, *Wm;
  complex *qcg, *q1cg, *scg, *x1cg, *zcg, *z1cg, *rcg, *Azcg;
  float  *eta, *rho, *gcv;
  float power, Jtotlast, Jtotprev, bb, tempfloat; 
  float *uaux;
  float dq1;
  float dq2;
  float *q1;
  float *q2;
  int costflag=0;
  int freqflag=0;  
  float *ph1;
  float *ph2;
  
  czero.r=czero.i=0;
  //eps1/=100.;
  
  d2=ealloc2complex(nh,nt);
  
  m2=ealloc2complex(nq,nt);//==> m2(nt x nq) 
  
  L=ealloc2complex(nq,nh);   //==> l(nh x nq)
  
  L1=ealloc2complex(nq1,nh);   //==> l(nq xnh)
  
  L2=ealloc2complex(nq2,nh);   //==> l(nq xnh)

  q1=ealloc1float(nq1);  

  q2=ealloc1float(nq2);

  dh=ealloc1float(nh);
  
  dh2=ealloc1float(nh);
  
  Jtot=ealloc1float(20);
  
  qcg=ealloc1complex(nq); 

  q1cg=ealloc1complex(nq);

  scg=ealloc1complex(nq);
  
  x1cg=ealloc1complex(nq);
  
  zcg=ealloc1complex(nq);

  z1cg=ealloc1complex(nq);
  
  rcg=ealloc1complex(nh);
  
  Azcg=ealloc1complex(nh);
  
  eta=ealloc1float(nq);
  
  rho=ealloc1float(nq);

  gcv=ealloc1float(nq);

  Wm=ealloc1complex(nq);
  
  dc=ealloc1complex(nh);
  
  uaux=ealloc1float(nq);
  
  ph1=ealloc1float(nh);

  ph2=ealloc1float(nh);

  float *Wm2=ealloc1float(nq);

  complex *m2small=ealloc1complex(nq2);
  
  for (iq=0;iq<nq1;iq++) Wm2[iq]=1;
  for (iq=0;iq<nq2;iq++) Wm2[nq1+iq]=1;
  

  ////////////////////////////////////////////////////

  radon_param_beam(fmax1,fmax2,h,nh,q, nq, qmin1,qmin2,q1,q2,ph1,ph2,nq1,nq2,depth1,
		   depth2,rtmethod1, rtmethod2,factor1,factor2,&dq1,&dq2);

  ///////////////////////////////////////////////////
  
  float *mtemp=ealloc1float(nq);
  
  for (ih=1;ih<nh;ih++) dh[ih]=h[ih]-h[ih-1];
  dh[0]=dh[1];

  for (i=0;i<=iter_end;i++) Jtot[i]=0.;
  fftgo0(-1,data,d2,nh,nt,dt,&nf0);

  if (0){
    ////////////////////////////////////////////
    //test for amplitude for fft
    //save_gather(data,nh,nt,0.004,"data1");
    //system("suxwigb  < data1 title=\"data1\" perc=100 & ");
    float **datatemp2=ealloc2float(nt,nh);
    fftback0(1,datatemp2,d2,nh,nt,dt,nf0);
    save_gather(datatemp2,nh,nt,0.004,"datatemp");
    system("suxwigb  < datatemp title=\"datatemp\" perc=100 & ");
    free2float(datatemp2);
    ////////////////////////////////////////////
  }

  nfreq=nf0/2;
  df=1/(nf0*dt);
  floatprint(df);
  int maxfreq1=(int) (fmax1/df);
  int maxfreq2=(int) (fmax2/df);

  fprintf(stderr,"maxfreq1=%d, maxfreq2=%d, , dt=%f, df=%f\n",maxfreq1,
	  maxfreq2,dt,df);

  const double  pi=acos(-1.);
  float test; 

  
  for (freq=1;freq<maxfreq1;freq++){
     w=2*pi*freq*df;
     wa=freqweight(freq,df,fmax1-10,fmax1);
     for (iq=0;iq<nq;iq++){
       Wm[iq].r=eps1;
       Wm[iq].i=0.0;
     }    
     
     if (rtmethod1==3) matrix(L1,ph1,q1,nh,nq1,w,dq1,rtmethod1);
     else  matrix(L1,h,q1,nh,nq1,w,dq1,rtmethod1);
     if (rtmethod2==3) matrix(L2,ph2,q2,nh,nq2,w,dq2,rtmethod2);
     else matrix(L2,h,q2,nh,nq2,w,dq2,rtmethod2);

     for (ih=0;ih<nh;ih++){
       for (iq=0;iq<nq1;iq++)
	 L[ih][iq]=L1[ih][iq];
       for (iq=0;iq<nq2;iq++) 
	 L[ih][iq+nq1]=L2[ih][iq];
     }	 
     
     //if (testadj) test=testadj_rad_f(L,LH);
     fprintf(stderr,"freq=%d===>2 oper\n",freq);
     power=rcdot(nh,d2[freq],d2[freq]);
     
     for (iter=1;iter<=iter_end;iter++){

       if (iter>1) for (iq=0;iq<nq;iq++) mtemp[iq]=abs(m2[freq][iq]); 
    
       for (iq=0;iq<nq;iq++)  {
	 if ((iter!=1)){
	   float qmodel=quest(quantil,nq,mtemp);
	   //qmodel=1;
       	   normu=abs(m2[freq][iq])/qmodel+eps2;
	   Wm[iq].r=1/normu+eps1;
	   //Wm[iq].r*=Wm2[iq];
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
       wtcgls(d2[freq],L,m2[freq],Wm,Wd,nh,nq,qcg,q1cg,scg,x1cg,zcg,
	      z1cg,rcg,Azcg,eta,rho,gcv,0,step,itercg);
     }
     for (iq=0;iq<nq;iq++) uaux[iq]=abs(m2[freq][iq]);            
     /////////////////////////////////////////////////
     //for (iq=0;iq<nq;iq++) m2[freq][iq]/=nq;      
     if ((wa<1)&&(wa>0)) for (iq=0;iq<nq1;iq++)  m2[freq][iq]*=wa;
     //for (iq=0;iq<nq;iq++) m2[freq][iq]=u[iq];          
  }
  for (freq=maxfreq1;freq<nfreq;freq++){
    for (iq=0;iq<nq1;iq++) m2[freq][iq]=czero;
  }     

  //////////////////////////////////////////////////////////////
  // For Ground roll,
  // Between maxfreq1 and freqmax2   we use only one operator
  for (freq=maxfreq1;freq<maxfreq2;freq++){
    w=2*pi*freq*df;
    wa=freqweight(freq,df,fmax2-10,fmax2);
    for (iq=0;iq<nq2;iq++){
      Wm[iq].r=eps1;
      Wm[iq].i=0.0;
    }    
     
    if (rtmethod2==3) matrix(L2,ph2,q2,nh,nq2,w,dq2,rtmethod2);
    else matrix(L2,h,q2,nh,nq2,w,dq2,rtmethod2);

     //if (testadj) test=testadj_rad_f(L,LH);
     fprintf(stderr,"freq=%d==>1 oper\n",freq);
     power=rcdot(nh,d2[freq],d2[freq]);
     
     for (iter=1;iter<=iter_end;iter++){
       
       if (iter>1) for (iq=0;iq<nq2;iq++) mtemp[iq]=abs(m2small[iq]); 
    
       for (iq=0;iq<nq2;iq++)  {
	 if ((iter!=1)){
	   float qmodel=quest(quantil,nq2,mtemp);
	   //qmodel=1;
       	   normu=abs(m2small[iq])/qmodel+eps2;
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
       wtcgls(d2[freq],L2,m2small,Wm,Wd,nh,nq2,qcg,q1cg,scg,x1cg,zcg,
	      z1cg,rcg,Azcg,eta,rho,gcv,0,step,itercg);
     }
     for (iq=0;iq<nq1;iq++){
       m2[freq][iq]=czero; 
       //fprintf(stderr,"m2[%d][%d].r=%f m2[%d][%d].i=%f\n",freq,iq,
       //       m2[freq][iq].r,freq,iq,m2[freq][iq].i);       
     }

     for (iq=nq1;iq<nq;iq++) m2[freq][iq]=m2small[iq];            
     //for (iq=0;iq<nq;iq++) m2[freq][iq]=czero;            

     //////////////////////////////////////////////////////////////
     
     for (iq=0;iq<nq;iq++) uaux[iq]=abs(m2[freq][iq]); 
     
     /////////////////////////////////////////////////
     //for (iq=0;iq<nq;iq++) m2[freq][iq]/=nq2;      
     if ((wa<1)&&(wa>0)) for (iq=nq1;iq<nq;iq++)  m2[freq][iq]*=wa;
     //for (iq=0;iq<nq;iq++) m2[freq][iq]=u[iq];          
  }

  
  

  for (iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  
  for (freq=maxfreq2;freq<nfreq;freq++){
    for (iq=0;iq<nq;iq++) m2[freq][iq]=czero;
  }     
  fprintf(stderr,"w=%f\n",w);
  exit=fftback0(1,model,m2,nq,nt,dt,nf0);
  if (exit!= EXIT_SUCCESS) err("***No succesfull fftback\n");
  
  if (1){
    save_gather(model,nq,nt,0.004,"model");
    system("sufft  < model | suamp | suxwigb title=\"d\" perc=99 & ");
  }

  free1complex(m2small);
  free1float(Wm2);
  free1float(q1);
  free1float(q2);
  free1float(ph1);
  free1float(ph2);
  free1float(mtemp);
   free1float(uaux);
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
   free1float(dh2);
   free1float(dh);
   free2complex(L2);
   free2complex(L1);
   free2complex(L);
   free2complex(m2);
   free2complex(d2);
   return;
}




void matrix(complex **L, float *h, float *q,int nh,int nq,float w, float dq, int rtmethod)
{

        register int ih;
	register int iq;  
        complex  arg;

        for (iq=0;iq<nq;iq++){
	  for (ih=0;ih<nh;ih++){
              arg.r=0;
              if ((rtmethod==1)||(rtmethod==3)) arg.i=w*h[ih]*q[iq];
              else if (rtmethod==2) arg.i=w*h[ih]*h[ih]*q[iq];
	      
	      L[ih][iq]=exp(-arg);///sqrt(nq*nh);
	  }
	}
	
        return;
}





































