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

float wpcgnr(void (*oper) (float *,float *,unsigned short **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned short **index, float tol, float step,int itercg,int restart)
 
{
  float normb,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
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

  float *r,*g,*s,*z,*w,*rho,*eta,rhold;

  int nx=nt*np;
  int ny=nt*nh;  

  if ((g=alloc1float(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((s=alloc1float(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((z=alloc1float(nx))==NULL)
    err("cannot allocate memory for scg\n");
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

  for (i=0;i<ny;i++) r[i]*=Wd[i]; 
  (*oper) (g,r,index,1,nt,nh,np,nsparse);
  for(i=0;i<nx;i++) z[i]=g[i]*M[i];
  normb=dot(nx,z,z);
  
  for(i=0;i<nx;i++) s[i]=z[i];
  rhold=tol*2;
  k=0;
  while((rhold>tol)&&(k<itercg)){
    k++;
    (*oper) (s,w,index,0,nt,nh,np,nsparse);
    for (i=0;i<ny;i++) w[i]*=Wd[i];

    alphanum=dot(nx,z,g);
    alphaden=dot(ny,w,w);
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-10 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }

    alpha=alphanum/alphaden;
    alpha*=step;

    for(i=0;i<nx;i++) x[i]+=alpha*s[i];
    for(i=0;i<ny;i++) r[i]-=alpha*w[i];  


 
    for (i=0;i<ny;i++) r[i]*=Wd[i];    
    (*oper) (g,r,index,1,nt,nh,np,nsparse);

    for(i=0;i<nx;i++) z[i]=g[i]*M[i];
    
    rho[k]=dot(nx,z,z)/normb;    
    fprintf(stderr,"rho[%d]=%e\n",k,rho[k]);      

    beta=dot(nx,z,g)/alphanum;
    
    for (i=0;i<nx;i++) s[i]=z[i]+beta*s[i];
    eta[k]=dot(nx,x,x);
  }

  if (1) save_vector(&rho[1],k,"rhofile");

  // Because the number of iterations is acting as regularization 
  // parameter the cost function is || AWm^{-1}m-b||^2

  J=dot(nx,z,z);  

  fprintf(stderr,"iter=%d,J=%f\n",k,J);
  free1float(eta);
  free1float(rho);
  free1float(w);
  free1float(r);
  free1float(z);
  free1float(s);
  free1float(g);

  return(J);
}


float wpcgnr(void (*oper)(float *,float *,float *,float *,float *,float **,int ,int ,int,int),int nt, int nh, int np, float *t, float *h, float *p, float *x,float *b,float *Wd, float *M, float **vel, float tol, float step, int itercg,int restart)
 
{
  float normb,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j;
  float J; // Cost function 
  // Temp pointers
  // r  Conjugate Gradient (Residual)
  // g  Gradient
  // z  Preconditioned gradient
  // s  Search direction
  // w  Conjugate search direction
  // M  precondtioner on data space
  // Wd model and data weights respectively.

  float *r,*g,*s,*z,*w,*rho,*eta,rhold;

  int nx=nt*np;
  int ny=nt*nh;  

  if ((g=alloc1float(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((s=alloc1float(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((z=alloc1float(nx))==NULL)
    err("cannot allocate memory for scg\n");
  if ((r=alloc1float(ny))==NULL)
    err("cannot allocate memory for rcg\n");
  if ((w=alloc1float(ny))==NULL)
    err("cannot allocate memory for Azcg\n");
  if ((rho=alloc1float(itercg+1))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((eta=alloc1float(itercg+1))==NULL)
    err("cannot allocate memory for eta\n"); 

  fprintf(stderr,"itercg=%d\n",itercg);
  //for (i=0;i<ny;i++) fprintf(stderr,"Wd[%d]=%f\n",i,Wd[i]);
  //for (i=0;i<nx;i++) fprintf(stderr,"M[%d]=%f\n",i,M[i]);


  if (restart){
    for (i=0;i<nx;i++) x[i]=0.;
    for (i=0;i<ny;i++) r[i]=Wd[i]*b[i];
  }
  else{
    (*oper) (x,t,h,p,r,vel,0,nt,nh,np);
    for (i=0;i<ny;i++) r[i]=Wd[i]*(b[i]-r[i]);
  }


  for (i=0;i<ny;i++) r[i]*=Wd[i]; 
  (*oper) (g,t,h,p,r,vel,1,nt,nh,np);

  normb=dot(nx,g,g);

  for(i=0;i<nx;i++) z[i]=g[i]*M[i];
  
  for(i=0;i<nx;i++) s[i]=z[i];
  rhold=tol*2;  
  k=0;
  while((rhold>tol)&&(k<itercg)){
    k++;
    (*oper) (s,t,h,p,w,vel,0,nt,nh,np);
    for (i=0;i<ny;i++) w[i]*=Wd[i];

    alphanum=dot(nx,z,g);
    alphaden=dot(ny,w,w);
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-7 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }

    alpha=alphanum/alphaden;
    alpha*=step;

    for(i=0;i<nx;i++) x[i]+=alpha*s[i];
    for(i=0;i<ny;i++) r[i]-=alpha*w[i];  


 
    for (i=0;i<ny;i++) r[i]*=Wd[i];    
    (*oper) (g,t,h,p,r,vel,1,nt,nh,np);

    rho[k]=dot(nx,g,g)/normb;


    for(i=0;i<nx;i++) z[i]=g[i]*M[i];
    
    beta=dot(nx,z,g)/alphanum;
    
    for (i=0;i<nx;i++) s[i]=z[i]+beta*s[i];

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

  J=dot(nx,z,z);

  fprintf(stderr,"iter=%d\n",k);
  free1float(eta);
  free1float(rho);
  free1float(w);
  free1float(r);
  free1float(z);
  free1float(s);
  free1float(g);
  return(J);
}


float wpcgnr(void (*oper) (float *,float *,unsigned int **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned int **index, float tol, float step,int itercg,int restart)
 
{
  float normb,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
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

  float *r,*g,*s,*z,*w,*rho,*eta,rhold;

  int nx=nt*np;
  int ny=nt*nh;  

  if ((g=alloc1float(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((s=alloc1float(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((z=alloc1float(nx))==NULL)
    err("cannot allocate memory for scg\n");
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


  for (i=0;i<ny;i++) r[i]*=Wd[i]; 
  (*oper) (g,r,index,1,nt,nh,np,nsparse);
  for(i=0;i<nx;i++) z[i]=g[i]*M[i];
  normb=dot(nx,z,z);
  
  for(i=0;i<nx;i++) s[i]=z[i];
  rhold=tol*2;
  k=0;
  while((rhold>tol)&&(k<itercg)){
    k++;
    (*oper) (s,w,index,0,nt,nh,np,nsparse);
    for (i=0;i<ny;i++) w[i]*=Wd[i];

    alphanum=dot(nx,z,g);
    alphaden=dot(ny,w,w);
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-10 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }

    alpha=alphanum/alphaden;
    alpha*=step;

    for(i=0;i<nx;i++) x[i]+=alpha*s[i];
    for(i=0;i<ny;i++) r[i]-=alpha*w[i];  
 
    for (i=0;i<ny;i++) r[i]*=Wd[i];    
    (*oper) (g,r,index,1,nt,nh,np,nsparse);

    for(i=0;i<nx;i++) z[i]=g[i]*M[i];
    
    rho[k]=dot(nx,z,z)/normb;    
    fprintf(stderr,"rho[%d]=%e\n",k,rho[k]);      

    beta=dot(nx,z,g)/alphanum;
    
    for (i=0;i<nx;i++) s[i]=z[i]+beta*s[i];
     
    eta[k]=dot(nx,x,x);

  }
  if (1){
    //float rhoaux;
    //for (i=1;i<=k;i++)  rhoaux=log10(rho[i]);
    save_vector(&rho[1],k,"rhofile");
  }

  // Because the number of iterations is acting as regularization 
  // parameter the cost function is || AWm^{-1}m-b||^2

  J=dot(nx,z,z);  

  fprintf(stderr,"iter=%d,J=%f\n",k,J);
  free1float(eta);
  free1float(rho);
  free1float(w);
  free1float(r);
  free1float(z);
  free1float(s);
  free1float(g);

  return(J);
}



float wpcgnr(void (*oper) (float *,float *,float **, int, int ,int ,int), int nt, int nh, int np, float *x,float *b,float *Wd, float *M, float **index, float tol, float step,int itercg,int restart)
 
{
  float normb,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
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

  float *r,*g,*s,*z,*w,*rho,*eta,rhold;

  int nx=nt*np;
  int ny=nt*nh;  

  if ((g=alloc1float(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((s=alloc1float(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((z=alloc1float(nx))==NULL)
    err("cannot allocate memory for scg\n");
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
    (*oper) (x,r,index,0,nt,nh,np);
    for (i=0;i<ny;i++) r[i]=Wd[i]*(b[i]-r[i]); 
  }


  for (i=0;i<ny;i++) r[i]*=Wd[i]; 
  (*oper) (g,r,index,1,nt,nh,np);
  for(i=0;i<nx;i++) z[i]=g[i]*M[i];
  normb=dot(nx,z,z);
  
  for(i=0;i<nx;i++) s[i]=z[i];
  rhold=tol*2;
  k=0;
  while((rhold>tol)&&(k<itercg)){
    k++;
    (*oper) (s,w,index,0,nt,nh,np);
    for (i=0;i<ny;i++) w[i]*=Wd[i];

    alphanum=dot(nx,z,g);
    alphaden=dot(ny,w,w);
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-10 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }

    alpha=alphanum/alphaden;
    alpha*=step;

    for(i=0;i<nx;i++) x[i]+=alpha*s[i];
    for(i=0;i<ny;i++) r[i]-=alpha*w[i];  


 
    for (i=0;i<ny;i++) r[i]*=Wd[i];    
    (*oper) (g,r,index,1,nt,nh,np);

    for(i=0;i<nx;i++) z[i]=g[i]*M[i];
    
    rho[k]=dot(nx,z,z)/normb;    
    //fprintf(stderr,"rho[%d]=%e\n",k,rho[k]);      
    fprintf(stderr,"rho[%d]=%e\n",k,rho[k]); 
    beta=dot(nx,z,g)/alphanum;
    
    for (i=0;i<nx;i++) s[i]=z[i]+beta*s[i];


     
    eta[k]=dot(nx,x,x);

  }
  if (1){
    //float rhoaux;
    //for (i=1;i<=k;i++)  rhoaux=log10(rho[i]);
    save_vector(&rho[1],k,"rhofile");
  }

  // Because the number of iterations is acting as regularization 
  // parameter the cost function is || AWm^{-1}m-b||^2

  J=dot(nx,z,z);  

  fprintf(stderr,"iter=%d,J=%f\n",k,J);
  free1float(eta);
  free1float(rho);
  free1float(w);
  free1float(r);
  free1float(z);
  free1float(s);
  free1float(g);

  return(J);
}


float wpcgnr(void (*oper) (float *,float *,float **, int, int, int ,int ,int), int nt, int nh, int np, int ns, float *x,float *b,float *Wd, float *M, float **index, float tol, float step,int itercg,int restart)
 
{
  float normb,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
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

  float *r,*g,*s,*z,*w,*rho,*eta,rhold;

  int nx=nt*np;
  int ny=nt*nh;  

  if ((g=alloc1float(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((s=alloc1float(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((z=alloc1float(nx))==NULL)
    err("cannot allocate memory for scg\n");
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
    (*oper) (x,r,index,0,nt,nh,np,ns);
    for (i=0;i<ny;i++) r[i]=Wd[i]*(b[i]-r[i]); 
  }


  for (i=0;i<ny;i++) r[i]*=Wd[i]; 
  (*oper) (g,r,index,1,nt,nh,np,ns);
  for(i=0;i<nx;i++) z[i]=g[i]*M[i];
  normb=dot(nx,z,z);
  
  for(i=0;i<nx;i++) s[i]=z[i];
  rhold=tol*2;
  k=0;
  while((rhold>tol)&&(k<itercg)){
    k++;
    (*oper) (s,w,index,0,nt,nh,np,ns);
    for (i=0;i<ny;i++) w[i]*=Wd[i];

    alphanum=dot(nx,z,g);
    alphaden=dot(ny,w,w);
    
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-10 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }

    alpha=alphanum/alphaden;
    alpha*=step;

    for(i=0;i<nx;i++) x[i]+=alpha*s[i];
    for(i=0;i<ny;i++) r[i]-=alpha*w[i];  
 
    for (i=0;i<ny;i++) r[i]*=Wd[i];    
    (*oper) (g,r,index,1,nt,nh,np,ns);

    for(i=0;i<nx;i++) z[i]=g[i]*M[i];
    
    rho[k]=dot(nx,z,z)/normb;    
    //fprintf(stderr,"rho[%d]=%e\n",k,rho[k]);      
    fprintf(stderr,"rho[%d]=%e\n",k,rho[k]); 
    beta=dot(nx,z,g)/alphanum;
    
    for (i=0;i<nx;i++) s[i]=z[i]+beta*s[i];
     
    eta[k]=dot(nx,x,x);

  }
  if (1){
    //float rhoaux;
    //for (i=1;i<=k;i++)  rhoaux=log10(rho[i]);
    save_vector(&rho[1],k,"rhofile");
  }

  // Because the number of iterations is acting as regularization 
  // parameter the cost function is || AWm^{-1}m-b||^2

  J=dot(nx,z,z);  

  fprintf(stderr,"iter=%d,J=%f\n",k,J);
  free1float(eta);
  free1float(rho);
  free1float(w);
  free1float(r);
  free1float(z);
  free1float(s);
  free1float(g);

  return(J);
}

















