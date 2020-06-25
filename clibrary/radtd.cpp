#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrarytd.h"

void radtd2(float *t, float *q, float *h, float *m,float *d,float eps,
float qmin,float qmax, float fmax, float *vel)

  /*
    subroutine: radtd.cpp
    RADON TRANSFORM IN TIME DOMAIN
    Three possible curves, linear, parabolic and hyperbolic
    are implemented in the file radon.cpp
    Inside this file the functions are called radonlin, radonpar, radonhyp
    At the beginning a pointer to function is declared and 
    associated with the chosen function. 
    This pointer is used as a reference to the operator 
    in WTCGLS, TESTADJ or any other function. 

    Daniel Trad
    E-mail: dtrad@geop.ubc.ca
  */
{
  extern int nt,nh,nq, nx, ny, method, iter_end, reorth, rtmethod;
  extern float dt,dh,dq,eps1,eps2,thres,theta;
  register int i;
  int  j, iter, ih, iq;
  float *Wm;
  float *Wd;
  float dx_max;
  float dx_av;
  float qmaxt;
  float test;
  extern float factor; // Used to compute dq, i.e., dq=factor*dq;
  extern float step;
  extern int itercg;
  extern int testadj;
  extern int taper;
  //////////////////////////////////////////////////////////////////////////
  void (*radon) (float *,float *,float *,float *,float *,int ,int ,int,int);
  void (*radon2) (float *,float *,float *,float *,float *,float *,int ,int ,int,int);  
  if (rtmethod==3) {
    radon=radonhyp;
    radon2=radonhyp;
  }
  else if (rtmethod==2) radon=radonpar;
  else if (rtmethod==1) radon=radonlin;
  
  /////////////////////////////////////////////////////////////////////////
 
  if ((Wm=ealloc1float(nt*nq))==NULL)
    fprintf(stderr,"***Sorry, space for Wm could not be allocated\n"); 
  if ((Wd=ealloc1float(nt*nh))==NULL)
    fprintf(stderr,"***Sorry, space for Wd could not be allocated\n"); 
  ///////////////////////////////////////////////////////////////////
  // Compute inv(CD)=WdT Wd       
  for (i=0;i<ny;i++) Wd[i]=1.;
  if (taper==1)
    for (ih=0;ih<5;ih++) 
      for (i=0;i<nt;i++){
	Wd[ih*nt+i]=1-exp(-(ih+.3));
	Wd[(nh-ih-1)*nt+i]=Wd[ih*nt+i];
      }
  
  //   Velocity axis //////////////////////////////////////////////////
  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);

  if (fmax==0) fmax=0.6/(2*dt); // sinc8 is valid only until 0.6 Nyquist 
  if ((rtmethod==1)||(rtmethod==2)) 
    radon_param(fmax,h,nh,dx_av,qmin,&qmaxt,&qmax,&dq,nq,rtmethod,factor);
  
  else dq = (qmax-qmin)/(nq-1);

  for (i=0;i<nq;i++) q[i] = qmin+i*dq;

  ///////////////////////////////////////////////////////////////////////  

  fprintf(stderr,"Inside radtd nq=%d, nt=%d, nh=%d dq=%f\n",nq,nt,nh,dq);
  
  for (i=0;i<nx;i++) m[i]=0.;
 
  if (method==0){ // Adjoint and tests
    radon(m,t,h,q,d,1,nt,nh,nq); 
    if (testadj) for (i=0;i<iter_end;i++) test=testadjop(radon,t,h,q,nt,nh,nq);
  }  

  /* Marfurt's method with Semblance analysis */
  else if (method==1) {
    //Semblance computation
    FILE *myfilep;
    if((myfilep=fopen("semblance","w"))==NULL)
      err("cannot open file=%s\n","semblance");    
    semblance(Wm,t,h,q,d,nt,nh,nq,dt);
    fwrite(Wm,sizeof(float),nx,myfilep);
    fclose(myfilep);
    
    // Mask computation
    if ((thres>=0.7)||(thres<=0.01)) thres=0.3;
    FILE *myfilep2;
    if((myfilep2=fopen("mask","w"))==NULL)
      err("cannot open file=%s\n","mask");
    /*Mask is given by losigm function*/
    for (i=0;i<nx;i++) Wm[i]=1+1./(1+exp(100*(Wm[i]-thres)+0.5));
    fwrite(Wm,sizeof(float),nx,myfilep2);
    fclose(myfilep2);

    for (iter=1;iter<=iter_end;iter++){
      fprintf(stderr,"iter_end=%d;iter=%d;\n",iter_end,iter);
      wtcgls(radon,nt,nh,nq,t,h,q,m,d,Wm,Wd,0,step,eps1,eps2,itercg);
      if (iter<iter_end) {
        //thres=0.04; 
	for (i=0;i<nx;i++) 
	  Wm[i]=1+1./(1+exp(100*(m[i]-thres)+0.5));
      }
    }
  }
  /* Sparse approach (Thompson, 1984) 
     First it is computed the adjoint model, and  a Covariance matrix
     is computed from it. WTCGLS solves the system 
     (LH*L + Cm) m = LH d
  */
  else if (method==2) {
    radon(m,t,h,q,d,1,nt,nh,nq);
    if (testadj) test=testadjop(radon,t,h,q,nt,nh,nq);
    for (j=1;j<=iter_end;j++){
      for (i=0;i<nx;i++) Wm[i]=1./(eps2+fabs(m[i]))+eps1;
      for (i=0;i<ny;i++) Wd[i]=1.;
      wtcgls(radon,nt,nh,nq,t,h,q,m,d,Wm,Wd,0,step,eps1,eps2,itercg);  
      }
  }
  /* LSQR does not work properly */
  else if (method==3) {
    lsqr(t,q,h,m,d,eps,reorth);
  }
  /* Method=4 is used for test */ 
  else if (method==4)  radonopi_id(m,t,h,q,d); 
  else if (method==5) {

    // Method 5 is similar to Thorson but the rho filter is applied as
    // a preconditioner to the adjoint model before start WTCGLS
    int lfilter=51;
    int np2=lfilter/2;	
    float *ptrace, *ptrace2;
    float *rho;
    if ((ptrace=alloc1float(nt))==NULL)
      fprintf(stderr,"***Sorry, space for ptrace could not be allocated\n");
    if ((ptrace2=alloc1float(nt))==NULL)
      fprintf(stderr,"***Sorry, space for ptrace could not be allocated\n");   
    if ((rho=ealloc1float(lfilter))==NULL)
      fprintf(stderr,"***Sorry, space for rho could not be allocated\n"); 
    rho_filter(lfilter,nt,dt,rho);
    //for (i=0;i<lfilter;i++) fprintf(stderr,"rho[%d]=%f\n",i,rho[i]);
    fprintf(stderr,"RHO filter......................\n");

    radon(m,t,h,q,d,1,nt,nh,nq);
    
    for (int iq=0;iq<nq;iq++){ 
      for (i=0;i<nt;i++) ptrace[i]=m[iq*nt+i];
      conv(nt,-np2,ptrace,lfilter,0,rho,nt,0,ptrace2);
      for (i=0;i<nt;i++) m[iq*nt+i]=ptrace2[i];
    }    
  }
  else if (method==6) {
    radon(m,t,h,q,d,1,nt,nh,nq);
    hilbert(nt*nq,m,Wm);
    for (j=1;j<=iter_end;j++){
      for (i=0;i<nx;i++) Wm[i]=1./(eps2+Wm[i])+eps1;
      wtcgls(radon,nt,nh,nq,t,h,q,m,d,Wm,Wd,0,step,eps1,eps2,itercg);  
    }
  }
  else if (method==7) {
    radon2(m,t,h,q,d,vel,1,nt,nh,nq);
    if (testadj) test=testadjop(radon2,t,h,q,vel,nt,nh,nq);
    for (j=1;j<=iter_end;j++){
      for (i=0;i<nx;i++) Wm[i]=1./(eps2+fabs(m[i]))+eps1;
      for (i=0;i<ny;i++) Wd[i]=1.;
      wtcgls(radon2,nt,nh,nq,t,h,q,m,d,Wm,Wd,vel,eps,step,eps1,eps2,itercg);  
    }
  }
  
  
  free1float(Wm);
  free1float(Wd);  
  return;
  
}









