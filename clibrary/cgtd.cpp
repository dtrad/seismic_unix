#include "su.h"
#include "clibrarytd.h"
#include <math.h>

void cgtd(float *t, float *q,float *h, float *x, float *y, float eps)
{ 
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq,eps2;
  float ee, bb, beta, gammam, gamma, alfa, rms, alfaden;
  int k,i,j;
  extern int itercg;
  float t0=0, *r, *g, *s, *ss;
  
  if ((r=ealloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for r could not be allocated\n");
  if ((g=ealloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for g could not be allocated\n");
  if ((s=ealloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for s could not be allocated\n");
  if ((ss=ealloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for ss could not be allocated\n");


  for (i=0;i<nx;i++) x[i]=0.;
  for (i=0;i<ny;i++) r[i]=y[i];
 
  rstack(t,q,h,g,r,1);
  bb=dot(ny,y,y);
  for (i=0;i<nx;i++) s[i]=g[i];
  gammam=dot(nx,g,g);
  k=0;
  
  fprintf(stderr,"eps=%e,bb=%e,gammam=%e\n",eps,bb,gammam);
  while ((k!=nx)&&(sqrt(gammam)>(eps*bb))&&(k<itercg)){
    k++;
    rstack(t,q,h,s,ss,0);
    alfaden=dot(ny,ss,ss);
    if (alfaden < 0.) err("alfaden=%e\n",alfaden);
    if (alfaden < eps ){ 
      fprintf(stderr,"gammam=%e,alfaden=%e,k=%d,eps=%e\n",
	      gammam,alfaden,k,eps);
      break;
    }
    alfa=gammam/alfaden;
 
    for(i=0;i<nx;i++) x[i]+=(alfa*s[i]);
    for(i=0;i<ny;i++) r[i]-=(alfa*ss[i]);  
    rstack(t,q,h,g,r,1);

    /*
    fprintf(stderr,"eps=%e,eps2=%e\n",eps,eps2);
    for (i=0;i<nx;i++){ 
      if (fabs(x[i]<eps)) g[i]-=(x[i]/eps);
      else if ((fabs(x[i]>eps))&&(fabs(x[i]<eps2)))  {
	g[i]-=(1./x[i]);
	//fprintf(stderr,"x[%d]=%e\n",i,x[i]);
      }
      else if (fabs(x[i]>eps2)) g[i]-=(x[i]/eps2);
    }
    */
    gamma=dot(nx,g,g);
    
    if (gammam < eps){
      fprintf(stderr,"gamma=%e,gammam=%e,k=%d,eps=%e\n",
	      gamma,gammam,k,eps); 
      break;
    }
    beta=gamma/gammam;
    gammam=gamma;
    ee=dot(ny,r,r);
    rms=sqrt(ee/ny);
    for(i=0;i<nx;i++) s[i]=g[i]+beta*s[i];
    fprintf(stderr,"k=%d, CGLS rms=%e\n",k,rms);           
  }
  free1float(r);
  free1float(g);
  free1float(s);
  free1float(ss);     
  return;
}








