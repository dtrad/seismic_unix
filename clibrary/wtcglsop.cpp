#include "su.h"
#include "clibrary.h"
#include "clibrarytd.h"
#include <math.h>

 /* 
     This function solves the system of equations 
     (FH WdT Wd FH + WmT Wm ) m = FH WdT Wd d
     oper is an operator implemented with the function
     void oper (float *,float *,float *,float *,float *,int ,int ,int,int)
     Example 
     void oper (float x*,float t*,float h*,float p*,float d*,int adj, int nt,
     int nh ,int np);
     When adj=0 ==> forward operator   (adjoint == False)
     When adj=1 ==> adjoint operator   (adjoint == True )
    
  */

void wtcgls(void (*oper) (float *,float *,float *,float *,float *,int ,int ,int,int), int nt, int nh, int np, float *t, float *h, float *p, float *x,float *b,float *Wm, float *Wd, float tol, float step, float eps1, float eps2, int itercg)
 
{
  float normb,dq,dq2,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j,in,num;
  // Temp pointers
  float *q,*q1,*s,*x1,*z,*z1,*r,*Az,*eta,*rho,*gcv;

  int nx=nt*np;
  int ny=nt*nh;  
  //////////////////////////////////////////////////////////////////////
  // These definitions can be different from other versions of wtcgls

  //////////////////////////////////////////////////////////////////////  


  if ((q=alloc1float(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((q1=alloc1float(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((s=alloc1float(nx))==NULL)
    err("cannot allocate memory for scg\n");
  if ((x1=alloc1float(nx))==NULL)
    err("cannot allocate memory for x1cg\n");
  if ((z=alloc1float(nx))==NULL)
    err("cannot allocate memory for zcg\n");
  if ((z1=alloc1float(nx))==NULL)
    err("cannot allocate memory for z1cg\n");
  if ((r=alloc1float(ny))==NULL)
    err("cannot allocate memory for rcg\n");
  if ((Az=alloc1float(ny))==NULL)
    err("cannot allocate memory for Azcg\n");
  if ((eta=alloc1float(nx))==NULL)
    err("cannot allocate memory for eta\n");
  if ((rho=alloc1float(nx))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((gcv=alloc1float(nx))==NULL)
     err("cannot allocate memory for gcv\n");   

  fprintf(stderr,"eps1=%f,eps2=%f,\n",eps1,eps2);
  for (i=0;i<nx;i++) x[i]=0.;
  normb=sqrt(dot(ny,b,b));

  for (i=0;i<ny;i++) r[i]=b[i]*Wd[i];
  (*oper) (s,t,h,p,r,1,nt,nh,np);
  
  nit=itercg;
  for(i=0;i<nx;i++){
    q1[i]=s[i]/Wm[i];
    q[i]=q1[i]/Wm[i];
  }
  for (i=0;i<nx;i++) z[i]=q[i];
  dq=dot(nx,s,q);
  for (i=0;i<nx;i++) z1[i]=q1[i];
  for(i=0;i<nx;i++) x1[i]=0.;       
  for (j=0;j<nit;j++){

    (*oper) (z,t,h,p,Az,0,nt,nh,np);
    for (i=0;i<ny;i++) Az[i]=Az[i]*Wd[i];

    alphanum=dq;
    alphaden=dot(ny,Az,Az);
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-7 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }
    alpha=alphanum/alphaden;
    alpha*=step;
    //fprintf(stderr,"j=%d,alpha=%e\n",j,alpha);         
    //// Update model u and residuals
    for(i=0;i<nx;i++) x[i]=x[i]+alpha*z[i];
    for(i=0;i<ny;i++) r[i]=r[i]-alpha*Az[i];  
    
    for (i=0;i<ny;i++) r[i]=r[i]*Wd[i];
    (*oper) (s,t,h,p,r,1,nt,nh,np);

    for(i=0;i<nx;i++){
      q1[i]=s[i]/Wm[i];
      q[i]=q1[i]/Wm[i];
    }
    dq2=dot(nx,s,q);
    beta=dq2/dq;
    //fprintf(stderr,"j=%d,beta=%e\n",j,beta);
    dq=dq2;
    for (i=0;i<nx;i++) z[i]=q[i]+beta*z[i];
    rho[j] = sqrt(dot(ny,r,r))/normb;
    fprintf(stderr,"rho[%d]=%e\n",j,rho[j]);      
    for (i=0;i<nx;i++) {
      x1[i]=x1[i]+alpha*z1[i]; 
      z1[i]=q1[i]+beta*z1[i];
    }
    eta[j]=sqrt(dot(nx,x1,x1));
    if ((tol==0) && (j>2)){ // GCV criteria
       in = j-1;
       for (i=1;i<=in;i++){
       num=(ny-(i-1))*(ny-(i-1)); 
       gcv[i]=(rho[i]*rho[i])/num;
       }      
       if (gcv[j-2]<gcv[j-1]){ 
         fprintf(stderr,"GCV Criteria was reached in iteration %d\n",j-1);
         nit = j-1;
         return;
       } 
       
       else if ((tol!=0) && (rho[j] < tol)){ 
        fprintf(stderr,"Convergence have been acheived at iteration # %d\n",j);
        return;
       }
    }          
  }

  fprintf(stderr,"j=%d\n",j);
  free1float(gcv);
  free1float(rho);
  free1float(eta);
  free1float(Az);
  free1float(r);
  free1float(z1);
  free1float(z);
  free1float(x1);
  free1float(s);
  free1float(q1);
  free1float(q);        
  return;
}

void wtcgls(void (*oper) (float *,float *,float *,float *,float *,float *,int ,int ,int,int), int nt, int nh, int np, float *t, float *h, float *p, float *x,float *b,float *Wm, float *Wd, float *vel,float tol, float step, float eps1, float eps2, int itercg)
 
{
  float normb,dq,dq2,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j,in,num;
  // Temp pointers
  float *q,*q1,*s,*x1,*z,*z1,*r,*Az,*eta,*rho,*gcv;

  int nx=nt*np; 
  int ny=nt*nh;  
  //////////////////////////////////////////////////////////////////////
  // These definitions can be different from other versions of wtcgls

  //////////////////////////////////////////////////////////////////////  
  // Normalize data to 1
  // float maxd=0; for (i=0;i<ny;i++) if (fabs(b[i]) > maxd) maxd=fabs(b[i]);
  // for (i=0;i<ny;i++) b[i]/=maxd;
  
  if ((q=alloc1float(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((q1=alloc1float(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((s=alloc1float(nx))==NULL)
    err("cannot allocate memory for scg\n");
  if ((x1=alloc1float(nx))==NULL)
    err("cannot allocate memory for x1cg\n");
  if ((z=alloc1float(nx))==NULL)
    err("cannot allocate memory for zcg\n");
  if ((z1=alloc1float(nx))==NULL)
    err("cannot allocate memory for z1cg\n");
  if ((r=alloc1float(ny))==NULL)
    err("cannot allocate memory for rcg\n");
  if ((Az=alloc1float(ny))==NULL)
    err("cannot allocate memory for Azcg\n");
  if ((eta=alloc1float(nx))==NULL)
    err("cannot allocate memory for eta\n");
  if ((rho=alloc1float(nx))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((gcv=alloc1float(nx))==NULL)
     err("cannot allocate memory for gcv\n");   

  fprintf(stderr,"eps1=%f,eps2=%f,tol=%f\n",eps1,eps2,tol);
  for (i=0;i<nx;i++) x[i]=0.;
  normb=sqrt(dot(ny,b,b));

  for (i=0;i<ny;i++) r[i]=b[i]*Wd[i];
  (*oper) (s,t,h,p,r,vel,1,nt,nh,np);
  
  nit=itercg;
  for(i=0;i<nx;i++){
    q1[i]=s[i]/Wm[i];
    q[i]=q1[i]/Wm[i];
  }
  for (i=0;i<nx;i++) z[i]=q[i];
  dq=dot(nx,s,q);
  for (i=0;i<nx;i++) z1[i]=q1[i];
  for(i=0;i<nx;i++) x1[i]=0.;       
  for (j=0;j<nit;j++){

    (*oper) (z,t,h,p,Az,vel,0,nt,nh,np);
    for (i=0;i<ny;i++) Az[i]=Az[i]*Wd[i];

    alphanum=dq;
    alphaden=dot(ny,Az,Az);
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-7 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }
    alpha=alphanum/alphaden;
    alpha*=step;
    //fprintf(stderr,"j=%d,alpha=%e\n",j,alpha);         
    //// Update model u and residuals
    for(i=0;i<nx;i++) x[i]=x[i]+alpha*z[i];
    for(i=0;i<ny;i++) r[i]=r[i]-alpha*Az[i];  
    
    for (i=0;i<ny;i++) r[i]=r[i]*Wd[i];
    (*oper) (s,t,h,p,r,vel,1,nt,nh,np);

    for(i=0;i<nx;i++){
      q1[i]=s[i]/Wm[i];
      q[i]=q1[i]/Wm[i];
    }
    dq2=dot(nx,s,q);
    beta=dq2/dq;
    //fprintf(stderr,"j=%d,beta=%e\n",j,beta);
    dq=dq2;
    for (i=0;i<nx;i++) z[i]=q[i]+beta*z[i];
    rho[j] = sqrt(dot(ny,r,r))/normb;
    fprintf(stderr,"rho[%d]=%e\n",j,rho[j]);      
    for (i=0;i<nx;i++) {
      x1[i]=x1[i]+alpha*z1[i]; 
      z1[i]=q1[i]+beta*z1[i];
    }
    eta[j]=sqrt(dot(nx,x1,x1));
    if ((tol==0) && (j>2)){ // GCV criteria
       in = j-1;
       for (i=1;i<=in;i++){
       num=(ny-(i-1))*(ny-(i-1)); 
       gcv[i]=(rho[i]*rho[i])/num;
       }      
       if (gcv[j-2]<gcv[j-1]){ 
         fprintf(stderr,"GCV Criteria was reached in iteration %d\n",j-1);
         nit = j-1;
         return;
       } 
       
       else if ((tol!=0) && (rho[j] < tol)){ 
        fprintf(stderr,"Convergence have been acheived at iteration # %d\n",j);
        return;
       }
    }          
  }

  fprintf(stderr,"j=%d\n",j);
  free1float(gcv);
  free1float(rho);
  free1float(eta);
  free1float(Az);
  free1float(r);
  free1float(z1);
  free1float(z);
  free1float(x1);
  free1float(s);
  free1float(q1);
  free1float(q);        
  return;
}

void wtcgls(void (*oper) (float *,float *,float *,float *,float *,float **,int ,int ,int,int), int nt, int nh, int np, float *t, float *h, float *p, float *x,float *b,float *Wm, float *Wd, float **vel,float tol, float step, float eps1, float eps2, int itercg)
 
{
  float normb,dq,dq2,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j,in,num;
  // Temp pointers
  float *q,*q1,*s,*x1,*z,*z1,*r,*Az,*eta,*rho,*gcv;

  int nx=nt*np; 
  int ny=nt*nh;  
  //////////////////////////////////////////////////////////////////////
  // These definitions can be different from other versions of wtcgls

  //////////////////////////////////////////////////////////////////////  
  // Normalize data to 1
  // float maxd=0; for (i=0;i<ny;i++) if (fabs(b[i]) > maxd) maxd=fabs(b[i]);
  // for (i=0;i<ny;i++) b[i]/=maxd;
  
  if ((q=alloc1float(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((q1=alloc1float(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((s=alloc1float(nx))==NULL)
    err("cannot allocate memory for scg\n");
  if ((x1=alloc1float(nx))==NULL)
    err("cannot allocate memory for x1cg\n");
  if ((z=alloc1float(nx))==NULL)
    err("cannot allocate memory for zcg\n");
  if ((z1=alloc1float(nx))==NULL)
    err("cannot allocate memory for z1cg\n");
  if ((r=alloc1float(ny))==NULL)
    err("cannot allocate memory for rcg\n");
  if ((Az=alloc1float(ny))==NULL)
    err("cannot allocate memory for Azcg\n");
  if ((eta=alloc1float(nx))==NULL)
    err("cannot allocate memory for eta\n");
  if ((rho=alloc1float(nx))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((gcv=alloc1float(nx))==NULL)
     err("cannot allocate memory for gcv\n");   

  fprintf(stderr,"eps1=%f,eps2=%f,tol=%f\n",eps1,eps2,tol);
  for (i=0;i<nx;i++) x[i]=0.;
  normb=sqrt(dot(ny,b,b));

  for (i=0;i<ny;i++) r[i]=b[i]*Wd[i];
  (*oper) (s,t,h,p,r,vel,1,nt,nh,np);
  
  nit=itercg;
  for(i=0;i<nx;i++){
    q1[i]=s[i]/Wm[i];
    q[i]=q1[i]/Wm[i];
  }
  for (i=0;i<nx;i++) z[i]=q[i];
  dq=dot(nx,s,q);
  for (i=0;i<nx;i++) z1[i]=q1[i];
  for(i=0;i<nx;i++) x1[i]=0.;       
  for (j=0;j<nit;j++){

    (*oper) (z,t,h,p,Az,vel,0,nt,nh,np);
    for (i=0;i<ny;i++) Az[i]=Az[i]*Wd[i];

    alphanum=dq;
    alphaden=dot(ny,Az,Az);
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-7 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }
    alpha=alphanum/alphaden;
    alpha*=step;
    //fprintf(stderr,"j=%d,alpha=%e\n",j,alpha);         
    //// Update model u and residuals
    for(i=0;i<nx;i++) x[i]=x[i]+alpha*z[i];
    for(i=0;i<ny;i++) r[i]=r[i]-alpha*Az[i];  
    
    for (i=0;i<ny;i++) r[i]=r[i]*Wd[i];
    (*oper) (s,t,h,p,r,vel,1,nt,nh,np);

    for(i=0;i<nx;i++){
      q1[i]=s[i]/Wm[i];
      q[i]=q1[i]/Wm[i];
    }
    dq2=dot(nx,s,q);
    beta=dq2/dq;
    //fprintf(stderr,"j=%d,beta=%e\n",j,beta);
    dq=dq2;
    for (i=0;i<nx;i++) z[i]=q[i]+beta*z[i];
    rho[j] = sqrt(dot(ny,r,r))/normb;
    fprintf(stderr,"rho[%d]=%e\n",j,rho[j]);      
    for (i=0;i<nx;i++) {
      x1[i]=x1[i]+alpha*z1[i]; 
      z1[i]=q1[i]+beta*z1[i];
    }
    eta[j]=sqrt(dot(nx,x1,x1));
    if ((tol==0) && (j>2)){ // GCV criteria
       in = j-1;
       for (i=1;i<=in;i++){
       num=(ny-(i-1))*(ny-(i-1)); 
       gcv[i]=(rho[i]*rho[i])/num;
       }      
       if (gcv[j-2]<gcv[j-1]){ 
         fprintf(stderr,"GCV Criteria was reached in iteration %d\n",j-1);
         nit = j-1;
         return;
       } 
       
       else if ((tol!=0) && (rho[j] < tol)){ 
        fprintf(stderr,"Convergence have been acheived at iteration # %d\n",j);
        return;
       }
    }          
  }

  fprintf(stderr,"j=%d\n",j);
  free1float(gcv);
  free1float(rho);
  free1float(eta);
  free1float(Az);
  free1float(r);
  free1float(z1);
  free1float(z);
  free1float(x1);
  free1float(s);
  free1float(q1);
  free1float(q);        
  return;
}

void wtcgls(void (*oper) (float *,float *,unsigned short **,int ,int ,int, int, int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wm, float *Wd,unsigned short **index,float tol, float step, float eps1, float eps2, int itercg)
 
{
  float normb,dq,dq2,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j,in,num;
  // Temp pointers
  float *q,*q1,*s,*x1,*z,*z1,*r,*Az,*eta,*rho,*gcv;

  int nx=nt*np; 
  int ny=nt*nh;  
  //////////////////////////////////////////////////////////////////////
  // These definitions can be different from other versions of wtcgls

  //////////////////////////////////////////////////////////////////////  
  // Normalize data to 1
  // float maxd=0; for (i=0;i<ny;i++) if (fabs(b[i]) > maxd) maxd=fabs(b[i]);
  // for (i=0;i<ny;i++) b[i]/=maxd;
  
  if ((q=alloc1float(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((q1=alloc1float(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((s=alloc1float(nx))==NULL)
    err("cannot allocate memory for scg\n");
  if ((x1=alloc1float(nx))==NULL)
    err("cannot allocate memory for x1cg\n");
  if ((z=alloc1float(nx))==NULL)
    err("cannot allocate memory for zcg\n");
  if ((z1=alloc1float(nx))==NULL)
    err("cannot allocate memory for z1cg\n");
  if ((r=alloc1float(ny))==NULL)
    err("cannot allocate memory for rcg\n");
  if ((Az=alloc1float(ny))==NULL)
    err("cannot allocate memory for Azcg\n");
  if ((eta=alloc1float(nx))==NULL)
    err("cannot allocate memory for eta\n");
  if ((rho=alloc1float(nx))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((gcv=alloc1float(nx))==NULL)
     err("cannot allocate memory for gcv\n");   

  fprintf(stderr,"eps1=%f,eps2=%f,tol=%f\n",eps1,eps2,tol);
  for (i=0;i<nx;i++) x[i]=0.;


  for (i=0;i<ny;i++) r[i]=b[i]*Wd[i];
  (*oper) (s,r,index,1,nt,nh,np,nsparse);   
  normb=sqrt(dot(nx,s,s));
  nit=itercg;
  for(i=0;i<nx;i++){
    q1[i]=s[i]/Wm[i];
    q[i]=q1[i]/Wm[i];
  }
  for (i=0;i<nx;i++) z[i]=q[i];
  dq=dot(nx,s,q);
  for (i=0;i<nx;i++) z1[i]=q1[i];
  for(i=0;i<nx;i++) x1[i]=0.;       
  for (j=0;j<nit;j++){

    (*oper) (z,Az,index,0,nt,nh,np,nsparse); 

    for (i=0;i<ny;i++) Az[i]=Az[i]*Wd[i];

    alphanum=dq;
    alphaden=dot(ny,Az,Az);
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1e-7 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      //break;
    }
    alpha=alphanum/alphaden;
    alpha*=step;
    //fprintf(stderr,"j=%d,alpha=%e\n",j,alpha);         
    //// Update model u and residuals
    for(i=0;i<nx;i++) x[i]=x[i]+alpha*z[i];
    for(i=0;i<ny;i++) r[i]=r[i]-alpha*Az[i];  
    
    for (i=0;i<ny;i++) r[i]=r[i]*Wd[i];
    (*oper) (s,r,index,1,nt,nh,np,nsparse); 

    for(i=0;i<nx;i++){
      q1[i]=s[i]/Wm[i];
      q[i]=q1[i]/Wm[i];
    }
    dq2=dot(nx,s,q);
    beta=dq2/dq;
    //fprintf(stderr,"j=%d,beta=%e\n",j,beta);
    dq=dq2;
    for (i=0;i<nx;i++) z[i]=q[i]+beta*z[i];
    rho[j] = sqrt(dot(nx,s,s))/normb;
    fprintf(stderr,"rho[%d]=%e\n",j,rho[j]);      
    for (i=0;i<nx;i++) {
      x1[i]=x1[i]+alpha*z1[i]; 
      z1[i]=q1[i]+beta*z1[i];
    }
    eta[j]=sqrt(dot(nx,x1,x1));
    if ((tol==0) && (j>2)){ // GCV criteria
       in = j-1;
       for (i=1;i<=in;i++){
       num=(ny-(i-1))*(ny-(i-1)); 
       gcv[i]=(rho[i]*rho[i])/num;
       }      
       if (gcv[j-2]<gcv[j-1]){ 
         fprintf(stderr,"GCV Criteria was reached in iteration %d\n",j-1);
         nit = j-1;
         return;
       } 
       
       else if ((tol!=0) && (rho[j] < tol)){ 
        fprintf(stderr,"Convergence have been acheived at iteration # %d\n",j);
        return;
       }
    }          
  }

  fprintf(stderr,"j=%d\n",j);
  free1float(gcv);
  free1float(rho);
  free1float(eta);
  free1float(Az);
  free1float(r);
  free1float(z1);
  free1float(z);
  free1float(x1);
  free1float(s);
  free1float(q1);
  free1float(q);        
  return;
}
