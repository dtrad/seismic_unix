#include "su.h"
#include <math.h>

void Atimesx(complex *b,complex **A,complex *x,int n,int m, int adj);
float rcdot(int n, complex *a, complex *b);
 /* 
     This function solves the system of equations 
     (WmT FH WdT M^{-1} Wd F Wm) m =  FH WdT Wd M^{-1} d
     oper is an operator implemented with the function
     where M is a  preconditioner acting on data space.
     void oper (float *,float *,float *,float *,float *,int ,int ,int,int)
     Example 
     void oper (float x*,float t*,float h*,float p*,float d*,int adj, int nt,
     int nh ,int np);
     When adj=0 ==> forward operator   (adjoint == False)
     When adj=1 ==> adjoint operator   (adjoint == True )
     Wd allows to downweight bad data by using large values. 
     In general is a diagonal size(ndata) 
     M is the preconditioner. Diagonal size(nmodel) Large values
     penalize, small values focus the solution to a desired model value.
     M changes the null space, W does not. 
     Hence prior information about the model must be implemented by M. 
     
     Taken from Yousef Saad, pag 260: 
     Iterative methods for sparse linear systems
     W has bee added to the original algorihm 
     Daniel Trad - UBC March-2000
    
  */

float mpcgne(complex **A, int n, int m, complex *x, complex *b,float *MI, float tol, 
	     float step,int itercg,int restart, float eps1, float eps2, int itermin)
{
  float normb,nit;
  float beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j;
  float J; // Cost function and square residuals
  
  // Temp pointers
  // r  Conjugate Gradient (Residual)
  // g  Gradient
  // z  Preconditioned gradient
  // s  Search direction
  // w  Conjugate search direction
  // M  precondtioner on data space
  // Wd model and data weights respectively.
  
  complex *u,*g,*r,*s,*w;
  float *rho,*eta,rhold;
  complex czero;
  int nx=m;
  int ny=n;  
  float res;
  float sigma=1./eps1;
  czero.r=czero.i=0;

  if ((u=alloc1complex(ny))==NULL) 
    err("cannot allocate memory for u\n");
  if ((g=alloc1complex(nx))==NULL) 
    err("cannot allocate memory for g\n");
  if ((s=alloc1complex(ny))==NULL)
    err("cannot allocate memory for s\n");
  if ((r=alloc1complex(ny))==NULL)
    err("cannot allocate memory for rcg\n");
  if ((w=alloc1complex(ny))==NULL)
    err("cannot allocate memory for Azcg\n");
  if ((rho=alloc1float(itercg+1))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((eta=alloc1float(itercg+1))==NULL)
    err("cannot allocate memory for eta\n"); 
  
  complex *s2=ealloc1complex(ny);
  

  normb=rcdot(ny,b,b);
  //fprintf(stderr,"norm=%f\n",normb);

  for (i=0;i<ny;i++) {
    u[i]=czero;
    s[i]=r[i]=b[i];
  }
  
  k=0;
  rhold=0;
  res=0;
  float resold=1;
  while(((res < resold)&&(k<itercg))||(k<itermin)){

    k++;
    Atimesx(s,A,g,n,m,1);
    for (i=0;i<nx;i++) g[i]*=MI[i];
    Atimesx(w,A,g,n,m,0);
    
    for (i=0;i<ny;i++) s2[i]=s[i]*sigma;
    for (i=0;i<ny;i++) w[i]+=s2[i];
    
    alphanum=rcdot(ny,r,r);
    alphaden=rcdot(ny,w,s);
    alpha=alphanum/alphaden;
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if(1)
      if (alpha < 1e-5 ){ 
	fprintf(stderr,"alphanum=%e,alphaden=%e,alpha=%e\n",
		alphanum,alphaden,alpha);
	break;
      }
    alpha*=step;
    
    for(i=0;i<ny;i++) u[i]+=alpha*s[i];
    for(i=0;i<ny;i++) r[i]-=alpha*w[i];  

    resold=res;

    res=rcdot(ny,r,r);

    rho[k]=res/normb;
    
    beta=res/alphanum;
    //fprintf(stderr,"rho[%d]=%f,beta=%e\n",k,rho[k],beta);    
    for (i=0;i<ny;i++) s[i]=r[i]+beta*s[i];
    eta[k]=rcdot(ny,u,u);
    //fprintf(stderr,"k=%d => res=%f resold=%f itercg=%d \n ", k, res, resold, itercg);
  }

  // Because the number of iterations is acting as regularization 
  // parameter the cost function is || AM^{-1}m-b||^2

  J=rcdot(ny,r,r);  
  Atimesx(u,A,x,n,m,1);
  for (i=0;i<nx;i++) x[i]*=MI[i];
  fprintf(stderr,"iter=%d\n",k);
  //fprintf(stderr,"iter=%d,J=%f\n",k,J);

  free1complex(s2);
  free1float(eta);
  free1float(rho);
  free1complex(w);
  free1complex(r);
  free1complex(s);
  free1complex(g);
  free1complex(u);

  return(J);
}



















