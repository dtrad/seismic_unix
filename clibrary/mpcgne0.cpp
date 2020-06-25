#include "su.h"
#include <math.h>

void Atimesx(float *b,float **A,float *x,int adj, int n,int m);
float dot(int n, float *a, float *b);
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

float mpcgne(float **A, int n, int m, float *x,float *b,float *MI, float tol, float step,int itercg,int restart)
  
{
  float normb,nit;
  double beta,betanum, betaden, alpha, alphanum, alphaden;
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
  
  float *u,*g,*r,*s,*w,*rho,*eta,rhold;
  
  int nx=m;
  int ny=n;  
  
  if ((u=alloc1float(ny))==NULL) 
    err("cannot allocate memory for u\n");
  if ((g=alloc1float(nx))==NULL) 
    err("cannot allocate memory for g\n");
  if ((s=alloc1float(ny))==NULL)
    err("cannot allocate memory for s\n");
  if ((r=alloc1float(ny))==NULL)
    err("cannot allocate memory for rcg\n");
  if ((w=alloc1float(ny))==NULL)
    err("cannot allocate memory for Azcg\n");
  if ((rho=alloc1float(itercg+1))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((eta=alloc1float(itercg+1))==NULL)
    err("cannot allocate memory for eta\n"); 

  normb=dot(ny,b,b);
  
  for (i=0;i<ny;i++) {
    u[i]=0.;
    s[i]=r[i]=b[i];
  }
  
  k=0;
  rhold=tol*2;
  while((rhold>tol)&&(k<itercg)){
    k++;
    Atimesx(s,A,g,1,n,m);
    for (i=0;i<nx;i++) g[i]*=MI[i];
    Atimesx(w,A,g,0,n,m);
    
    alphanum=dot(ny,r,r);
    alphaden=dot(ny,w,s);
    alpha=alphanum/alphaden;
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if(0)
      if (alphaden < 1e+10 ){ 
	fprintf(stderr,"alphanum=%e,alphaden=%e,alpha=%e\n",
		alphanum,alphaden,alpha);
	//break;
      }
    //alpha*=step;
    
    for(i=0;i<ny;i++) u[i]+=step*alpha*s[i];
    for(i=0;i<ny;i++) r[i]-=step*alpha*w[i];  
    
    rho[k]=dot(ny,r,r)/normb;
    
    beta=dot(ny,r,r)/alphanum;
    //fprintf(stderr,"rho[%d]=%f,beta=%e\n",k,rho[k],beta);    
    for (i=0;i<ny;i++) s[i]=r[i]+beta*s[i];
    eta[k]=dot(ny,u,u);

  }

  // Because the number of iterations is acting as regularization 
  // parameter the cost function is || AM^{-1}m-b||^2

  J=dot(ny,r,r);  
  Atimesx(u,A,x,1,n,m);
  for (i=0;i<nx;i++) x[i]*=MI[i];

  fprintf(stderr,"iter=%d,J=%f\n",k,J);
  free1float(eta);
  free1float(rho);
  free1float(w);
  free1float(r);
  free1float(s);
  free1float(g);
  free1float(u);

  return(J);
}



















