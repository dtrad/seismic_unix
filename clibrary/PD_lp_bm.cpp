#include "su.h"

void radon_PD_lp_bm(float **d, float **m, float *t, float *q, float *h, int nt, 
		    int nv, int nh, int iter_max)
{
  float *c;
  int *irow;
  int *icol;
  float PrI;
  float DI;
  float DG;
  int iter;
  float delta, gamma;
  int nbig;

  nbig=nt*nh*nv;

  A=ealloc1float(nbig);
  irow=ealloc1int(nbig);
  icol=ealloc1int(nbig);
  c=ealloc1float(nx);
  
  na=NMO(nh,h,nt,t,nq,q,A,icol,irow);

  PD_lp_bm(nt*nq,nt*nh,na,c,A,m[0],d[0],irow,icol,iter,delta,gamma,PrI,DI,DG,iter_max,
	   nbig);


  free1float(A);
  free1float(c);
  free1int(irow);
  free1int(icol);

}


//     Minimize c'x +gamma^2x'x+p'p
//     subject Ax+delta.p=b, x.ge.0

int PD_lp_bm(int nx, int ny, int na, float *c, float *A, float *x, float *b, 
	     int *irow,int *icol,int iter, float delta, float gamma, float PrI, float DI,
	     float DG,int iter_max, int nbig)
{
  
  int PrI_t=1.e-1;
  int DI_t=1.e-1;
  int DG_t=1.e-1;

  float *D;
  float *t;
  float *r;
  float *z;
  float *y;
  float *rhs;
  float mu;
  float gamma;
  float delta;
  float rhod;
  float rhop;
  float x1;
  float z1;
  float *xd;
  float *xp;
  float x1;
  float z1;
  float *xp;
  float *yp;
  float *dx;
  float *dy;
  float *dz;

  // size of the model
  D=ealloc1float(nx);
  v=ealloc1float(nx);
  t=ealloc1float(nx);
  z=ealloc1float(nx);
  xp=ealloc1float(nx);
  dx=ealloc1float(nx);
  dz=ealloc1float(nx);

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
  iter=0;
  
  for(k=1;k<iter_max;k++){
    iter++;
    matmul0(1,na,A,np,xp,nd,y,irow,icol);
    for (ix=0;ix<nx;ix++){
          t[ix]=c[ix]+gamma*gamma*x[ix]-z[ix]-xp[ix];
          v[ix]=mu-x[ix]*z[ix];
          d[ix]=1./(z[ix]/x[ix]+gamma*gamma);
    }

    matmul0(0,na,A,np,x,nd,yp,irow,icol);
    for (iy=0;iy<ny;iy++) r[iy]=b[iy]-yp[iy]-delta*delta*y[iy];

    //     right side term of the system
    for (ix=0;ix<nx;ix++) xp[ix]=d[ix]*(v[ix]/x[ix]-t[ix]);
    matmul0 (0,na,A,np,xp,nd,yp,irow,icol);
    for (iy=0;iy<ny;iy++) rhs[iy]=r[iy]-yp[iy];

    //     solve ADA'dy=rhs
    CG_solver(A,D,na,np,nd,irow,icol,dy,rhs,delta,itercg);

    //    set dx and dz 

    matmul0(1,na,A,np,xp,nd,dy,irow,icol)
    fprintf(stderr,"itercg=%d",itercg);
    for (ix=0;ix<nx;ix++){
      dx[ix]=d[ix]*xp[ix]+d[ix]*(v[ix]/x[ix]-t[ix]);
      dz[ix]=v[ix]/x[ix]-z[ix]*dx[ix]/x[ix];
    }

    //    primal and dual step sizes 
       
    x1=x[1]+dx[1];
    z1=z[1]+dz[1];
    kx=1;
    kz=1;
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
    if (x1<0.) step_p=-x[kx]/dx[kx];

    if(z1>=0.) step_d=1;
    if(z1<0.) step_d=-z(kz)/dz(kz);

    rhop=.99*step_p;
    rhod=.99*step_d;

    for (ix=0;ix<nx;ix++){
      x[ix]=x[ix]+rhop*dx[ix];
      z[ix]=z[ix]+rhod*dz[ix];
    }

    for (iy=0;iy<ny;iy++) y[iy]=y[iy]+rhod*dy[iy];

    mu=(1.-min1(rhop,rhod,0.99))*mu;

    //  primal-dual conditions are tested here

    PrI=dot(nd,r,r)/(dot(np,x,x)+1.);
    DI=dot(np,t,t)/(dot(nd,y,y)+1.);
    DG=dot(np,z,x)/(dot(np,x,x)+1.);

    if ( (PrI<PrI_t) && (DI < DI_t) && (DG < DG_t) ) break;

  }


  free1float(D);
  free1float(v);
  free1float(t);
  free1float(z);
  free1float(xp);
  free1float(dx);
  free1float(dz);

  free1float(r);
  free1float(y);
  free1float(rhs);
  free1float(yp);
  free1float(dy);

  return(SUCCESS);

}


void matmul0(conj,nb,bb,nx,x,ny,y,irow,icol)
{
  //    Sparse matrix multiplication....
  
  if (conj.eq.0){
    for (iy=0;iy<ny;iy++) y[iy]=0;
    for (k=0;k<nb;k++) y[irow[k]]+=bb[k]*x[icol[k]);
  }
  else{
    for (ix=0;ix<nx;ix++) x[ix]=0;
    for (k=0;k<nb;k++) x[icol[k]]+=bb[k]*y[irow[k]];
  }
  return;
}        
 

float min3(float a, float b, float c)
{
  float min;

  if (a<b) min=a; else min=b;
  
  if (c<min) min=c;
  
  return(min);
}




//     ----------------------------------------------------------------


int NMO(int nh, float *h, int nt, float *t, int nv, float *v, float *A, int *icol, 
	int *irow)
{
  //   The NMO operator is load in the matrix A,icol,irow,with ndata elements.
  
  tol=1.e-5;
  float dt,t0,v0,dv,h0,dh,t,tau,v,h;
  ny=nt*nh;
  nx=nt*nv;

  for (itau=0;itau<nt;itau++){
    tau=t[itau];
    for (ih=0;ih<nh;ih++){      
      for (iv=0;iv<nv;iv++){
	v=v0+iv*dv;
	t=(sqrt(tau*tau+h[ih]*h[ih]/(v[iv]*v[iv]))-t[0])/dt;
	it=0.5+t;
	l=l+1;
	j =ih*nt+it;
	k =iv*nt+itau;
	irow[l]=j;
	icol[l]=k;
	A[l]=1.;
      }	
    }
  }   
  return(l);
}




int  CG_solver(float *A,float *D,int na,int np,int nd,float *irow,float *icol,
	       float *x, float *b, float delta)
{
  /******************************************
    Solve the system
    ADA'x=y whenre A is sparse and D is
    diagonal
  *******************************************/      

  float tol=1.e-3;
  float *r;
  float *p;
  float *w;
  float *xp;
  
  int np=nx;
  int nd=ny;
  
  r=ealloc1float(ny);
  p=eallo1float(ny);
  w=ealloc1float(ny);
  xp=ealloc1float(nx);
  
  for (iy=0;iy<ny;iy++){
    x[iy]=0;
    r[iy]=b[iy];
  }  

  rho=dot(ny,r,r);
  
  for (iter=1;iter<iter_max;iter++){
    if(iter==1) for (iy=0;iy<ny;iy++) p[iy]=r[iy];
    else{
      beta=rho/rho_old;
      for (iy=0;iy<ny;iy++) p[iy]=r[iy]+beta*p[iy];
    }
         
    //        xp=A'p  ***    D.xp  *** A.D.A'p

    matmul0 (1,na,A,np,xp,nd,p,irow,icol);
    for (ix=0;ix<nx;ix++) xp[ix]=D[ix]*xp[ix];
    matmul0 (0,na,A,np,xp,nd,w,irow,icol);
               
    //       Regularization

    for (iy=0;iy<ny;iy++) w[iy]=w[iy]+p[iy]*delta**2;
    alpha=rho/dot(nd,p,w);
    for (iy=0;iy<ny;iy++){              
      x[iy]=x[iy]+alpha*p[iy];
      r[iy]=r[iy]-alpha*w[iy];
    }
    rho_old=rho;
    rho=dot(nd,r,r);
    if (sqrt(rho/nd)<=tol) return(iter);
  }
  return(iter);
}



             









