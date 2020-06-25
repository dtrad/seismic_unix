#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrarytd.h"

void velanal(float *t, float *q, float *h, float *m,float *d,float eps,
float Vmin,float Vmax)

  /*
    subroutine: radtd.cpp
    HYPERBOLIC RADON TRANSFORM IN TIME DOMAIN
    Mauricio D. Sacchi, Daniel Trad
    E-mail: dtrad@geop.ubc.ca
    *********************************************************
    GENERAL DESCRIPTION:
    Given the seimic data  computes a hyperbolic stack using CG.
  */
{
  extern int nt,nh,nq, nx, ny, method, iter_end, reorth;
  extern float dt,dh,dq,eps1,eps2,thres,theta;
  int i, j, iter;
  float t0=0;
  float *ss;
  
  if ((ss=ealloc1float(nt*nq))==NULL)
    fprintf(stderr,"***Sorry, space for ss could not be allocated\n");  
  //   Velocity axis 
  dV = (Vmax-Vmin)/(nq-1);
  for (i=0;i<nq;i++) q[i] = qmin+i*dq;
  
  //   Time axis 
  for (i=0;i<nt;i++) t[i]=t0+i*dt;
  fprintf(stderr,"Inside radtd nq=%d, nt=%d, nh=%d\n",nq,nt,nh);
  // Read binary data
  //  Compute a velocity panel via the adjoint
  
  for (i=0;i<nx;i++) m[i]=0.;
  //  Invert the velocity panel using conjugate gradients 
  if (method==0)  radonopi(m,t,h,q,d);

  
  else if (method==1) {
    //Semblance computation
    FILE *myfilep;
    if((myfilep=fopen("semblance","w"))==NULL)
      err("cannot open file=%s\n","semblance");    
    semblance(ss,t,h,q,d,nt,nh,nq,dt);
    fwrite(ss,sizeof(float),nx,myfilep);
    fclose(myfilep);
    
    // Mask computation
    if ((thres>=0.7)||(thres<=0.01)) thres=0.3;
    FILE *myfilep2;
    if((myfilep2=fopen("mask","w"))==NULL)
      err("cannot open file=%s\n","mask");
    /*Mask is given by losigm function*/
    for (i=0;i<nx;i++) ss[i]=1+1./(1+exp(100*(ss[i]-thres)+0.5));
    fwrite(ss,sizeof(float),nx,myfilep2);
    fclose(myfilep2);
    for (iter=1;iter<=iter_end;iter++){
      fprintf(stderr,"iter_end=%d;iter=%d;\n",iter_end,iter);  
      wtcglstd(t,q,h,m,d,ss,0);
      if (iter<iter_end) {
        thres=0.04; 
	for (i=0;i<nx;i++) 
	  ss[i]=1+1./(1+exp(100*(m[i]-thres)+0.5));
      }
    }
  }
  else if (method==2) {
    radonopi(m,t,h,q,d);
    for (j=1;j<=iter_end;j++){
      for (i=0;i<nx;i++) ss[i]=1./(eps2+fabs(m[i]))+eps1;
      wtcglstd(t,q,h,m,d,ss,0);
      }
  }
  else if (method==3) {
    lsqr(t,q,h,m,d,eps,reorth);
  }
  
  free1float(ss);
  return;
  
}









