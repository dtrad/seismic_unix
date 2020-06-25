#include "su.h"
#include "clibrarytd.h"



void velanrt(float *t, float *q, float *h, float *m,float *d,float eps,
float Vmin,float Vmax)

  /*
    subroutine: velanrt.cpp
    Velocity analysis with HYP. RADON TRANSFORM IN TIME DOMAIN
    Daniel Trad
    E-mail: dtrad@geop.ubc.ca
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
  dq = (Vmax-Vmin)/(nq-1);
  for (i=0;i<nq;i++) q[i] = Vmin+i*dq;
  
  //   Time axis 
  for (i=0;i<nt;i++) t[i]=t0+i*dt;
  fprintf(stderr,"Inside radtd nq=%d, nt=%d, nh=%d\n",nq,nt,nh);
  // Read binary data
  //  Compute a velocity panel via the adjoint
  
  for (i=0;i<nx;i++) m[i]=0.;
  //  Invert the velocity panel using conjugate gradients 
  if (method==0)  velopi(m,t,h,q,d);
  else if (method==2) {
    velopi(m,t,h,q,d);
    for (j=1;j<=iter_end;j++){
      for (i=0;i<nx;i++) ss[i]=1./(eps2+fabs(m[i]))+eps1;
      wtcglsvel(t,q,h,m,d,ss,0);
      }
  }  
  free1float(ss);
  return;
  
}









