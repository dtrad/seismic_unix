#include "su.h"
#include "clibrarytd.h"
#include <math.h>

 /* 
     Preconditioned Conjugate Gradient method for the Normal Residuals (PCGNR)
     This function solves the system of equations 

     Wd A Wmi Wm m = Wd b 

     (for Wm inevertible),
     using least squares method to solve u= Wm m that gives

     (WimT AT WdT Wd A Wim) u =  WimT AT WdT Wd  d

     m= Wmi m;

     This can be also thoght as the transfromation to standard form of 

     mT WmT Wm m + (b - A m)^T WdT Wd (b - A x) 

     as AS= A Wmi ;   ms = Wm m;  bs = b 
     
     oper is an operator implemented with the function
     void oper (float *,float *,float *,float *,float *,int ,int ,int,int)
     Example 
     void oper (float x*,float t*,float h*,float p*,float d*,int adj, int nt,
     int nh ,int np);

     When adj=0 ==> forward operator   (adjoint == False)
     When adj=1 ==> adjoint operator   (adjoint == True )

     where Wmi is the inverse  preconditioner acting on model space.
     where Wd is the preconditioner acting on data space.

     Wd allows to downweight bad data by using large values. 
     In general is a diagonal size(ndata) 
     Wmi is the diagonal model inverse preconditioner. 
     
     Large values focus the solution to a desired model value.
     Reference
     Iterative methods for sparse linear systems

     Daniel Trad - UBC April-2000
    
  */

float pcgnr(void (*oper) (float *,float *,unsigned short **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *Wmi, unsigned short **index, float tol, float step,int itercg,int restart)
 
{
  double normb,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j;
  float J; // Cost function and square residuals

  // Temp pointers
  // r  Residual for original system
  // g  Gradient (residual of LS system)   
  // z  Preconditioned gradient
  // s  Conjugate vectors (search direction)
  // w  auxiliar vector
  // Wmi Inverse  preconditioner on data space
  // Wd model and data weights respectively.

  float *r,*g,*s,*w,*rho,*eta,rhold;
  double temp;
  int nx=nt*np;
  int ny=nt*nh;  

  if ((g=alloc1float(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((s=alloc1float(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((r=alloc1float(ny))==NULL)
    err("cannot allocate memory for rcg\n");
  if ((w=alloc1float(ny))==NULL)
    err("cannot allocate memory for Azcg\n");
  if ((rho=alloc1float(itercg+1))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((eta=alloc1float(itercg+1))==NULL)
    err("cannot allocate memory for eta\n"); 

  //for (i=0;i<ny;i++) fprintf(stderr,"Wd[%d]=%f\n",i,Wd[i]);
  //for (i=0;i<nx;i++) fprintf(stderr,"M[%d]=%f\n",i,M[i]);



  if (restart){
    for (i=0;i<nx;i++) x[i]=0.;
    for (i=0;i<ny;i++) r[i]=Wd[i]*b[i];
  }
  else{
    (*oper) (x,r,index,0,nt,nh,np,nsparse);
    for (i=0;i<ny;i++) r[i]=Wd[i]*(b[i]-r[i]); 
  }



  if (0) for (i=0;i<ny;i++) r[i]*=Wd[i]; 
  (*oper) (g,r,index,1,nt,nh,np,nsparse);
  float gg0=dot(nx,g,g); fprintf(stderr,"gg0=%f\n",gg0);
  for(i=0;i<nx;i++) g[i]*=Wmi[i];
  float gg=dot(nx,g,g); fprintf(stderr,"gg=%f\n",gg);

  // for(i=0;i<nx;i++) { g[i]*=sqrt(gg0/gg);Wmi[i]*=sqrt(gg0/gg); }
  // float gg=dot(nx,g,g); fprintf(stderr,"gTg=%f\n",gg);

  normb=gg;

  for(i=0;i<nx;i++) s[i]=g[i];
  rhold=tol*2;
  k=0;
  while((rhold>tol)&&(k<itercg)){
    k++;
    for (i=0;i<nx;i++) s[i]*=Wmi[i];
    (*oper) (s,w,index,0,nt,nh,np,nsparse);
    if (0) for (i=0;i<ny;i++) w[i]*=Wd[i];

    alphanum=gg;
    alphaden=dot(ny,w,w);
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-10 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,k=%d\n",
	      alphanum,alphaden,k);
      break;
    }

    alpha=alphanum/alphaden;
    //alpha*=step;

    for(i=0;i<nx;i++) x[i]+=alpha*s[i];
    for(i=0;i<ny;i++) r[i]-=alpha*w[i];  
 
    if (0) for (i=0;i<ny;i++) r[i]*=Wd[i];    
    (*oper) (g,r,index,1,nt,nh,np,nsparse);
    for(i=0;i<nx;i++) g[i]*=Wmi[i];

    gg=dot(nx,g,g);
    rho[k]=gg/normb;
    beta=gg/alphanum;
    
    for (i=0;i<nx;i++) s[i]=g[i]+beta*s[i];

    fprintf(stderr,"rho[%d]=%e\n",k,rho[k]);      
     
    eta[k]=dot(nx,x,x);

  }
  if (1){
    //float rhoaux;
    //for (i=1;i<=k;i++)  rhoaux=log10(rho[i]);
    save_vector(&rho[1],k,"rhofile");
  }


  // Because the number of iterations is acting as regularization 
  // parameter the cost function is || AWm^{-1}m-b||^2

  J=dot(ny,r,r)+dot(nx,x,x);  

  for (i=0;i<nx;i++) x[i]*=Wmi[i];

  fprintf(stderr,"iter=%d,J=%f\n",k,J);
  free1float(eta);
  free1float(rho);
  free1float(w);
  free1float(r);
  free1float(s);
  free1float(g);

  return(J);
}


























