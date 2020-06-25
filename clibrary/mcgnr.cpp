#include "su.h"
#include "clibrary.h"
#include "clibrarytd.h"
#include <math.h>

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
    
  */

void mcgnr(void (*oper) (float *,float *,unsigned short **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wm, float *Wd, float *M, unsigned short **index, float tol, float step, float eps1, float eps2, int itercg,int restart)
 
{
  float normb,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j;
  // Temp pointers
  // two letters conjugate variables,
  // one letter variables.
  // rr Conjugate residuals
  // rp Precondtioned Conjugate residuals
  // r  Residuals for the normal equation system (Gradient)
  // s  Search direction
  // ww Conjugte direction
  // M precondtioner on data space
  // Wm , Wd model and data weights respectively.
  float *r,*s,*rp,*rr,*ww,*rho,*eta,rhold;

  int nx=nt*np;
  int ny=nt*nh;  

  if ((r=alloc1float(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((s=alloc1float(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((rp=alloc1float(ny))==NULL)
    err("cannot allocate memory for scg\n");
  if ((rr=alloc1float(ny))==NULL)
    err("cannot allocate memory for rcg\n");
  if ((ww=alloc1float(ny))==NULL)
    err("cannot allocate memory for Azcg\n");
  if ((rho=alloc1float(itercg+1))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((eta=alloc1float(itercg+1))==NULL)
    err("cannot allocate memory for eta\n"); 

  fprintf(stderr,"eps1=%f,eps2=%f,\n",eps1,eps2);

  normb=dot(ny,b,b);

  if (restart) for (i=0;i<nx;i++) x[i]=0.;  
  else{
    for(i=0;i<nx;i++) x[i]/=Wm[i];
    (*oper) (x,rr,index,0,nt,nh,np,nsparse);
    for (i=0;i<ny;i++) rr[i]*=Wd[i];
  }

  for(i=0;i<ny;i++) rr[i]=b[i]-rr[i];

  for(i=0;i<ny;i++) rp[i]=rr[i]/M[i];

  for (i=0;i<ny;i++) rp[i]*=Wd[i]; 
  (*oper) (r,rp,index,1,nt,nh,np,nsparse);
  for(i=0;i<nx;i++) r[i]/=Wm[i];
  
  for(i=0;i<nx;i++) s[i]=r[i];

  k=0;
  while((rhold>tol)&&(k<itercg)){
    k++;

    for(i=0;i<nx;i++) s[i]=s[i]/Wm[i];
    (*oper) (s,ww,index,0,nt,nh,np,nsparse);
    for (i=0;i<ny;i++) ww[i]*=Wd[i];

    alphanum=dot(nx,r,r);
    alphaden=dot(ny,ww,ww);
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < tol ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }

    alpha=alphanum/alphaden;
    alpha*=step;

    for(i=0;i<nx;i++) x[i]+=alpha*s[i];
    for(i=0;i<ny;i++) rr[i]-=alpha*ww[i];  

    rho[k]=dot(ny,rr,rr)/normb;

    for(i=0;i<ny;i++) rp[i]=rr[i]/M[i];

    for (i=0;i<ny;i++) rp[i]*=Wd[i];    
    (*oper) (r,rp,index,1,nt,nh,np,nsparse);
    for(i=0;i<nx;i++) r[i]/=Wm[i];
    
    beta=dot(nx,r,r)/alphanum;
    
    for (i=0;i<nx;i++) s[i]=z[i]+beta*s[i];

    fprintf(stderr,"rho[%d]=%e\n",k,rho[k]);      
    
    eta[k]=dot(nx,x,x);

  }

  fprintf(stderr,"iter=%d\n",k);
  free1float(eta);
  free1float(rho);
  free1float(ww);
  free1float(r);
  free1float(rr);
  free1float(z);
  free1float(s);
  return;
}

