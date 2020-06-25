#include "su.h"
#include "dan.h"
#include "radonl1.h"
#define USE_TRACE

void radonl1_freq_Karmarkar(float **d, float **m, float *t, float *q, float *h,
			    int nt, int nq, int nh, int iter_end, int itercg, 
			    int method)
{
  float *c;
  float *A;
  int *irow;
  int *icol;
  float PrI=0;
  float DI=0;
  float DG=0;
  int iter_ext=0;
  float delta=1.e0;
  float gamma=1.e-4;
  int nbig;
  int na;
  int iq;
  int vel=0;
  int nq2;
  float **m2;
  int it;

  // iter_end: maximum ext iterations
  // itercg: maximum cg iterations
  // iter_ext: external iterations performed
  // iter_int: internal iterations performed

  nq2=2*nq;

  if (method==1) vel=0;   // q is square slowness
  else if (!method || method==2) vel=1; // q is velocities 
  
  nbig=nt*nh*nq2;

  A=ealloc1float(nbig);
  irow=ealloc1int(nbig);
  icol=ealloc1int(nbig);
  c=ealloc1float(nt*nq2);
  m2=ealloc2float(nt,nq2);

  memset( (void *) m2[0], (int) '\0', nq2 * nt *FSIZE);
  for (iq=0;iq<nt*nq2;iq++) c[iq]=1;
  



  fprintf(stderr,"na=%d,method=%d\n",na,method);
  else  PD_lp_bm(nt*nq2,nt*nh,na,c,A,m2[0],d[0],irow,icol,iter_ext,delta,
		 gamma,PrI,DI,DG,iter_end,itercg,nbig);
  
  for (it=0;it<nt;it++) 
    for (iq=0;iq<nq;iq++) 
      m[iq][it]=m2[iq][it]-m2[nq+iq][it];

  
  matmul0(0,A,nq/2,m[0],nh,d[0]);  


  free2float(m2);
  free1float(A);
  free1float(c);
  free1int(irow);
  free1int(icol);

}


//     Minimize c'x +gamma^2x'x+p'p
//     subject Ax+delta.p=b, x.ge.0

int PD_lp_bm(int nx, int ny, float *c, float *A, float *x, float *b, 
	     int iter_ext, float delta, float gamma, 
	     float PrI, float DI, float DG,int iter_end, int itercg, int nbig)
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

    matmul0(0,na,A,np,x,nd,yp);
    for (iy=0;iy<ny;iy++) r[iy]=b[iy]-yp[iy]-delta*delta*y[iy];

    //     right side term of the system
    for (ix=0;ix<nx;ix++) xp[ix]=D[ix]*(v[ix]/x[ix]-t[ix]);
    matmul0 (0,na,A,np,xp,nd,yp);
    for (iy=0;iy<ny;iy++) rhs[iy]=r[iy]-yp[iy];

    //     solve ADA'dy=rhs
    iter_int=CG_solver(A,D,na,np,nd,dy,rhs,delta,itercg);

    //TRACE("iter=%d\n",iter)

    //    set dx and dz 

    matmul0(1,na,A,np,xp,nd,dy);
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
  free1float(r);
  free1float(y);
  free1float(rhs);
  free1float(yp);
  free1float(dy);

  return EXIT_SUCCESS;

}

void matmul0(int adj, float **A, int nx, float *x, int ny, float *y)
{
/*
  From Claerbout's PVI
  subroutine matmult( adj, bb,nx, x, ny, y)
  integer ix, iy, adj, nx, ny
  real bb(ny,nx), x(nx), y(ny)
  call adjnull(adj, 0,x,nx,  y,ny)
  do ix= 1, nx {
  do iy= 1, ny {
  if( adj == 0 ) y(iy) = y(iy) + bb(iy,ix) * x(ix)
  else x(ix) = x(ix) + bb(iy,ix) * y(iy)
  }}
  return; end

*/
  int ix,iy;
  
  if (!adj) for (iy=0;iy<ny;iy++) y[iy]=0.0;
  else for (ix=0;ix<nx;ix++) x[ix]=0.0;
  
  for (ix=0;ix<nx;ix++)
    for (iy=0;iy<ny;iy++)
      if (!adj) y[iy]+=(A[iy][ix]*x[ix]);
      else x[ix]+=(A[iy][ix]*y[iy]);
  
  return;
}



float min3(float a, float b, float c)
{
  float min;

  if (a<b) min=a; else min=b;
  
  if (c<min) min=c;
  
  return(min);
}




int  CG_solver(float **A,float *D,int na,int np,int nd,
	       float *x, float *b, float delta, int iter_end)
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
  int ix, iy;
  int nx=np;
  int ny=nd;
  float rho;
  float rho_old;
  float beta;
  float alpha;
  int iter;

  //fprintf(stderr,"iter_end=%d,nd=%d,np=%d\n",iter_end,nd,np);


  r=ealloc1float(ny);
  p=ealloc1float(ny);
  w=ealloc1float(ny);
  xp=ealloc1float(nx);
  
  for (iy=0;iy<ny;iy++){
    x[iy]=0;
    r[iy]=b[iy];
  }  

  rho=dot(ny,r,r);
  rho_old=1;

  for (iter=1;iter<=iter_end;iter++){
    if(iter==1) for (iy=0;iy<ny;iy++) p[iy]=r[iy];
    else{
      beta=rho/rho_old;
      for (iy=0;iy<ny;iy++) p[iy]=r[iy]+beta*p[iy];
    }
         
    //        xp=A'p  ***    D.xp  *** A.D.A'p

    matmul0 (1,A,np,xp,nd,p);
    for (ix=0;ix<nx;ix++) xp[ix]*=D[ix];
    matmul0 (0,A,np,xp,nd,w);
               
    //       Regularization

    for (iy=0;iy<ny;iy++) w[iy]+=p[iy]*delta*delta;
    alpha=rho/dot(nd,p,w);
    for (iy=0;iy<ny;iy++){              
      x[iy]+=alpha*p[iy];
      r[iy]-=alpha*w[iy];
    }
    rho_old=rho;
    rho=dot(nd,r,r);
    //fprintf(stderr,"rho=%f\n",rho);

    if (sqrt(rho/nd)<tol) return(iter);
  }
  return(iter);
}



             









