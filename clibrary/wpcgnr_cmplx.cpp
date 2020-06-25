#include "su.h"
#include "Complex.h"
#include "clibrary.h"
#include <math.h>
#include "clibrarytd.h"

int wpcgnr(void (*oper) (complex *,complex *, complex  **, int, int, int), int nh, int np, complex *x,complex *b,float *Wd, float *M, complex  **L, float tol, float step,int itercg,int restart)
{
  float normb,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j;

  // Temp pointers
  // r  Conjugate Gradient (Residual)
  // g  Gradient
  // z  Preconditioned gradient
  // s  Search direction
  // w  Conjugate search direction
  // M  precondtioner on data space
  // Wd model and data weights respectively.

  complex *r,*g,*s,*z,*w;
  float *rho,*eta,rhold;
  complex czero;czero.r=czero.i=0;
  
  int nx=np;
  int ny=nh;  

  if ((g=alloc1complex(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((s=alloc1complex(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((z=alloc1complex(nx))==NULL)
    err("cannot allocate memory for scg\n");
  if ((r=alloc1complex(ny))==NULL)
    err("cannot allocate memory for rcg\n");
  if ((w=alloc1complex(ny))==NULL)
    err("cannot allocate memory for Azcg\n");
  if ((rho=alloc1float(itercg+1))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((eta=alloc1float(itercg+1))==NULL)
    err("cannot allocate memory for eta\n"); 

  //for (i=0;i<ny;i++) fprintf(stderr,"Wd[%d]=%f\n",i,Wd[i]);
  //for (i=0;i<nx;i++) fprintf(stderr,"M[%d]=%f\n",i,M[i]);
  normb=rcdot(ny,b,b);

  if (restart){
    for (i=0;i<nx;i++) x[i]=czero;
    for (i=0;i<ny;i++) r[i]=Wd[i]*b[i];
  }
  else{
    (*oper) (x,r,L,0,nh,np);
    for (i=0;i<ny;i++) r[i]=Wd[i]*(b[i]-r[i]); 
  }


  for (i=0;i<ny;i++) r[i]*=Wd[i]; 
  (*oper) (g,r,L,1,nh,np);
  for(i=0;i<nx;i++) z[i]=g[i]*M[i];
  
  for(i=0;i<nx;i++) s[i]=z[i];
  rhold=tol*2;
  k=0;
  while((rhold>tol)&&(k<itercg)){
    k++;
    (*oper) (s,w,L,0,nh,np);
    for (i=0;i<ny;i++) w[i]*=Wd[i];

    alphanum=rcdot(nx,z,g);
    alphaden=rcdot(ny,w,w);
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-10 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,k=%d\n",
	      alphanum,alphaden,k);
      return(k);
    }

    alpha=alphanum/alphaden;
    alpha*=step;

    for(i=0;i<nx;i++) x[i]+=alpha*s[i];
    for(i=0;i<ny;i++) r[i]-=alpha*w[i];  

    rho[k]=rcdot(ny,r,r)/normb;
 
    for (i=0;i<ny;i++) r[i]*=Wd[i];    
    (*oper) (g,r,L,1,nh,np);

    for(i=0;i<nx;i++) z[i]=g[i]*M[i];
    
    beta=rcdot(nx,z,g)/alphanum;
    
    for (i=0;i<nx;i++) s[i]=z[i]+beta*s[i];

    //fprintf(stderr,"rho[%d]=%e\n",k,rho[k]);      
    if (k>1) if (rho[k]>=(rho[k-1]-1e-3)) return(k);

    eta[k]=rcdot(nx,x,x);

  }
  //if (1){
    //float rhoaux;
    //for (i=1;i<=k;i++)  rhoaux=log10(rho[i]);
  //save_vector(&rho[1],k,"rhofile");
  //}

  fprintf(stderr,"iter=%d\n",k);
  free1float(eta);
  free1float(rho);
  free1complex(w);
  free1complex(r);
  free1complex(z);
  free1complex(s);
  free1complex(g);
  return(k);
}






















