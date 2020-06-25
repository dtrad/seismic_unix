#include "su.h"
#include "dan.h"
#include "radonl1.h"
#define USE_TRACE

int karmarker(int nx, int ny, int na, float *c, float *A, float *x, float *b, 
	     int *irow,int *icol,int iter_ext, float delta, float gamma, float PrI, 
	     float DI, float DG,int iter_end, int itercg, int nbig)
{
  
  float PrI_t=1.e-1;
  float DI_t=1.e-1;
  float DG_t=1.e-1;
  int ix, iy, k;
  int np=nx;
  int nd=ny;
  int kx;
  int kz;
  float *D;
  float *t;
  float *r;
  float *z;
  float *y;
  float *rhs;
  float mu;
  float rhod;
  float rhop;
  float x1;
  float z1;
  float *xp;
  float *yp;
  float *dx;
  float *dy;
  float *dz;
  float *v;
  float step_p;
  float step_d;
  int iter_int=0;


  // size of the model
  D=ealloc1float(nx);
  v=ealloc1float(nx);
  t=ealloc1float(nx);
  z=ealloc1float(nx);
  xp=ealloc1float(nx);
  dx=ealloc1float(nx);
  dz=ealloc1float(nx);
  v=ealloc1float(nx);

  // size of data

  r=ealloc1float(ny);
  y=ealloc1float(ny);
  rhs=ealloc1float(ny);
  yp=ealloc1float(ny);
  dy=ealloc1float(ny);

  //  set inital variables
  for (ix=0;ix<nx;ix++) x[ix]=z[ix]=1;
  for (iy=0;iy<ny;iy++) y[iy]=1;
  mu=5.e-2;
  iter_ext=0;
  
  for(k=1;k<=iter_end;k++){
    iter_ext++;
    matmul0(1,na,A,np,xp,nd,y,irow,icol);
    for (ix=0;ix<nx;ix++){
          t[ix]=c[ix]+gamma*gamma*x[ix]-z[ix]-xp[ix];
          v[ix]=mu-x[ix]*z[ix];
          D[ix]=1./(z[ix]/x[ix]+gamma*gamma);
    }

    matmul0(0,na,A,np,x,nd,yp,irow,icol);
    for (iy=0;iy<ny;iy++) r[iy]=b[iy]-yp[iy]-delta*delta*y[iy];

    //     right side term of the system
    for (ix=0;ix<nx;ix++) xp[ix]=D[ix]*(v[ix]/x[ix]-t[ix]);
    matmul0 (0,na,A,np,xp,nd,yp,irow,icol);
    for (iy=0;iy<ny;iy++) rhs[iy]=r[iy]-yp[iy];

    //     solve ADA'dy=rhs
    iter_int=CG_solver(A,D,na,np,nd,irow,icol,dy,rhs,delta,itercg);

    //TRACE("iter=%d\n",iter)

    //    set dx and dz 

    matmul0(1,na,A,np,xp,nd,dy,irow,icol);
    fprintf(stderr,"k=%d,iter_int=%d\n",k,iter_int);

    for (ix=0;ix<nx;ix++){
      dx[ix]=D[ix]*xp[ix]+D[ix]*(v[ix]/x[ix]-t[ix]);
      dz[ix]=v[ix]/x[ix]-z[ix]*dx[ix]/x[ix];
    }

    //    primal and dual step sizes 
       
    x1=x[0]+dx[0];
    z1=z[0]+dz[0];
    kx=0;
    kz=0;
    for (ix=1;ix<nx;ix++){ 
      if (z1>=(z[ix]+dz[ix])){
	z1=z[ix]+dz[ix];
	kz=ix;
      }
      if (x1>=(x[ix]+dx[ix])){
	x1=x[ix]+dx[ix];
	kx=ix;
      } 
    }

    if (x1>=0.) step_p=1.;
    else step_p=-x[kx]/dx[kx];

    if(z1>=0.) step_d=1;
    else step_d=-z[kz]/dz[kz];

    rhop=.99*step_p;
    rhod=.99*step_d;

    for (ix=0;ix<nx;ix++){
      x[ix]+=rhop*dx[ix];
      z[ix]+=rhod*dz[ix];
    }

    for (iy=0;iy<ny;iy++) y[iy]+=rhod*dy[iy];

    mu=(1.-min3(rhop,rhod,0.99))*mu;

    //fprintf(stderr,"rhop=%f,rhod=%f,min3=%f\n",rhop,rhod,min3(rhop,rhod,0.99));

    //  primal-dual conditions are tested here

    PrI=dot(nd,r,r)/(dot(np,x,x)+1.);
    DI=dot(np,t,t)/(dot(nd,y,y)+1.);
    DG=dot(np,z,x)/(dot(np,x,x)+1.);

    fprintf(stderr,"PrI=%f, DI=%f, DG=%f \n",PrI,DI,DG);
    if ( (PrI < PrI_t) && (DI < DI_t) && (DG < DG_t) ){
      fprintf(stderr,"Crazy program\n");
      break;
    }
  }


  free1float(D);
  free1float(v);
  free1float(t);
  free1float(z);
  free1float(xp);
  free1float(dx);
  free1float(dz);
  free1float(v);
  free1float(r);
  free1float(y);
  free1float(rhs);
  free1float(yp);
  free1float(dy);

  return EXIT_SUCCESS;

}
