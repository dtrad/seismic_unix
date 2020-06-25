#include "su.h"
#include "radonlogbar.h"

#define USE_TRACE

int logbarrier_interface(complex *d,complex **L,complex *m, int nh, int nq, logb_param par)
{
  float *mr=0;
  float *dr=0;
  float **LR=0;
  int nq2;
  int nh2;
  int iter_done;  

  nh2=2*nh;
  nq2=2*nq;

  // Transform the complex system of equation in a real one 
  dr=ealloc1float(nh2);
  mr=ealloc1float(nq2);
  LR=ealloc2float(nq2,nh2);

  complex2real(d,L,dr,LR,nh,nq);

  //Atimesx2(dr,LR,mr,nh2,nq2,1);  

  iter_done=logbarrier_method(dr,LR,mr,nq2,nh2,par);  

  //Atimesx2(dr,LR,mr,nh2,nq2,0);
  real2complex(d,dr,m,mr,nh,nq);  

  free2float(LR);
  free1float(dr);
  free1float(mr);

  return(iter_done);

}


int logbarrier_method(float *d, float **L, float *m, int nq, int nh, logb_param par)
{
  int iq, ih;
  int nq2;
  float *m2;
  float **L2;
  int iter_done;
  //static float m0[400];
  //static int flag=0;

  // iter_end: maximum ext iterations
  // itercg: maximum cg iterations
  // iter_done: external iterations performed

  nq2=2*nq;
  m2=ealloc1float(nq2);
  L2=ealloc2float(nq2,nh);

  memset( (void *) m2, (int) '\0', nq2 * FSIZE);

  
  //  if (flag==0)  for (iq=0;iq<nq2;iq++) m2[iq]=0.01;
  //  else for (iq=0;iq<nq2;iq++) m2[iq]=0.1*m0[iq];
   
  for (ih=0;ih<nh;ih++) 
    for (iq=0;iq<nq;iq++){
      L2[ih][iq]=L[ih][iq];
      L2[ih][nq+iq]=-L[ih][iq];
    }
  // test
  if (0){
    Atimesx2(d,L2,m2,nh,nq2,1);  
    Atimesx2(d,L2,m2,nh,nq2,0);
  }  
  iter_done=PD_lp_bm(d,L2,m2,nq2,nh,par);
  //flag=1;for (iq=0;iq<nq2;iq++) m0[iq]=m2[iq];

  //test
  if (0) matmul0(d,L2,m2,1,nh,nq2);
  //////
  for (iq=0;iq<nq;iq++) m[iq]=m2[iq]-m2[nq+iq];
  
  free1float(m2);
  free2float(L2);

  return(iter_done);

}


//     Minimize c'x +gamma^2x'x+p'p
//     subject Ax+delta.p=b, x.ge.0

int PD_lp_bm(float *b, float **A, float *x, int nx, int ny, logb_param par)   
{
  

  float PI_t=par.PI_t;
  float DI_t=par.DI_t;
  float GI_t=par.GI_t;
  float delta=par.delta;
  float gamma=par.gamma;
  int itercg=par.itercg;
  int iter_end=par.iter_end;
  int iter_done=0;
  float PIn;
  float DIn;
  float GIn;
  
  int ix, iy, k;
  int np=nx;
  int nd=ny;
  int kx;
  int kz;
  float *c;
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
  c=ealloc1float(nx);
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

  for (ix=0;ix<nx;ix++) x[ix]=1;
  for (ix=0;ix<nx;ix++) z[ix]=c[ix]=1;
  for (iy=0;iy<ny;iy++) y[iy]=1;
  mu=5.e-2;
  iter_done=0;
  
  for(k=1;k<=iter_end;k++){
    iter_done++;
    matmul0(y,A,xp,1,nd,np);
    for (ix=0;ix<nx;ix++){
          t[ix]=c[ix]+gamma*gamma*x[ix]-z[ix]-xp[ix];
          v[ix]=mu-x[ix]*z[ix];
          D[ix]=1./(z[ix]/x[ix]+gamma*gamma);
    }

    matmul0(yp,A,x,0,nd,np);
    for (iy=0;iy<ny;iy++) r[iy]=b[iy]-yp[iy]-delta*delta*y[iy];

    //     right side term of the system
    for (ix=0;ix<nx;ix++) xp[ix]=D[ix]*(v[ix]/x[ix]-t[ix]);
    matmul0 (yp,A,xp,0,nd,np);
    for (iy=0;iy<ny;iy++) rhs[iy]=r[iy]-yp[iy];

    //     solve ADA'dy=rhs
    iter_int=CG_solver(A,D,np,nd,dy,rhs,delta,itercg);

    fprintf(stderr,"iter=%d\n",iter_int);

    //    set dx and dz 

    matmul0(dy,A,xp,1,nd,np);
    //fprintf(stderr,"k=%d,iter_int=%d\n",k,iter_int);

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

    PIn=dot(nd,r,r)/(dot(np,x,x)+1.);
    DIn=dot(np,t,t)/(dot(nd,y,y)+1.);
    GIn=dot(np,z,x)/(sqrt(dot(np,x,x)*dot(np,z,z))+1.);

    //fprintf(stderr,"PrI=%f, DI=%f, DG=%f \n",PrI,DI,DG);
    if ( (PIn < PI_t) && (DIn < DI_t) && (GIn < GI_t) ){
      fprintf(stderr,"Convergence \n");
      break;
    }
  }
  free1float(c);
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

  return(iter_done);

}

void matmul0(float *y, float **A, float *x, int adj, int ny, int nx)
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




int  CG_solver(float **A,float *D,int np,int nd,
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

    matmul0 (p,A,xp,1,nd,np);
    for (ix=0;ix<nx;ix++) xp[ix]*=D[ix];
    matmul0 (w,A,xp,0,nd,np);
               
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

    if (sqrt(rho/nd)<tol ) return(iter);
  }
  return(iter);
}



             









