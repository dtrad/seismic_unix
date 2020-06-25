#include "su.h"
#include "clibrarytd.h"
#include <math.h>

void cgtdhr(float *t, float *q,float *h, float *x, float *y,
	  float *r,float *g, float *s, float *ss, float *D, float eps)
{ 
  extern int nt,nh,nq, nx, ny, iter_end;
  extern float dt,dh,dq,eps2;
  float ee, bb, beta, gammam, gamma, alfa, rms, alfaden, sum;
  int k,i,j, iter;
  extern int itercg;
  float *adj;
  if ((adj=ealloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for adj could not be allocated\n");   

  for (i=0;i<nx;i++) x[i]=0.;
  for (i=0;i<ny;i++) r[i]=y[i];
 
  rstack(t,q,h,g,r,1); // g=m_adjoint 
  gammam=dot(nx,g,g);
  for (i=0;i<nx;i++) adj[i]=s[i]=g[i];
  bb=gammam;
 

  for (iter=1;iter<iter_end;iter++){
    fprintf(stderr,"eps=%e,bb=%e,gammam=%e,iter_end=%d\n",
	    eps,bb,gammam,iter_end);
    for (i=0;i<nx;i++){ 
      if (fabs(x[i]<eps)) D[i]=1./eps;
      else if ((fabs(x[i]>eps))&&(fabs(x[i]<eps2)))  
	D[i]=(1./(x[i]*x[i]));
      else if (fabs(x[i]>eps2)) D[i]=1./eps2;
    }
    for (i=0;i<ny;i++) r[i]=adj[
    k=0;
    //while ((k!=nx)&&(sqrt(gammam)>(eps*bb))&&(k<itercg)){
    while (k<itercg){
      k++;
      rstack(t,q,h,s,ss,0);
      for (sum=0,i=0;i<nx;i++) sum+=s[i]*D[i]*s[i];
      alfaden=dot(ny,ss,ss);alfaden+=sum;
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

    
      fprintf(stderr,"eps=%e,eps2=%e\n",eps,eps2);
      for (i=0;i<nx;i++) g[i]-=(x[i]*D[i]);
    
      gamma=dot(nx,g,g);
    
      if (gammam < 1e-10){
	fprintf(stderr,"Betaden < 1e-12 gamma=%e,gammam=%e,k=%d,eps=%e\n",
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
  }
  free1float(adj);       
  return;
}








