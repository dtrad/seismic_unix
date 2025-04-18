#include "su.h"
#include "clibrarytd2.h"
#include <math.h>
void wtcgls(double *t, double *qaxis, double *h,double *x,double *b,double tol)
  
{
  double normb,dq,dq2,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j,in,num;
  extern int itercg,nx,ny;
  double *Qp,*q,*q1,*s,*x1,*z,*z1,*r,*Az,*eta,*rho,*gcv;
  extern double step,eps1,eps2;
  if ((Qp=alloc1double(nx))==NULL) 
    err("cannot allocate memory for Qp\n");
  if ((q=alloc1double(nx))==NULL) 
    err("cannot allocate memory for qcg\n");
  if ((q1=alloc1double(nx))==NULL)
    err("cannot allocate memory for q1cg\n");
  if ((s=alloc1double(nx))==NULL)
    err("cannot allocate memory for scg\n");
  if ((x1=alloc1double(nx))==NULL)
    err("cannot allocate memory for x1cg\n");
  if ((z=alloc1double(nx))==NULL)
    err("cannot allocate memory for zcg\n");
  if ((z1=alloc1double(nx))==NULL)
    err("cannot allocate memory for z1cg\n");
  if ((r=alloc1double(ny))==NULL)
    err("cannot allocate memory for rcg\n");
  if ((Az=alloc1double(ny))==NULL)
    err("cannot allocate memory for Azcg\n");
  if ((eta=alloc1double(nx))==NULL)
    err("cannot allocate memory for eta\n");
  if ((rho=alloc1double(nx))==NULL)
    err("cannot allocate memory for rho\n");  
  if ((gcv=alloc1double(nx))==NULL)
     err("cannot allocate memory for gcv\n");   
  for (i=0;i<nx;i++) 
    Qp[i]=1./(fabs(x[i])+eps2)+eps1;
  //if (sqrt(fabs(x[i]))>1e-3) Qp[i]=1./sqrt(fabs(x[i]));
  //  else Qp[i]=1./1e-3;
  fprintf(stderr,"eps1=%f,eps2=%f,\n",eps1,eps2);
  for (i=0;i<nx;i++) x[i]=0.;
  normb=sqrt(dot(ny,b,b));
  for (i=0;i<ny;i++) r[i]=b[i];
  rstack(t,qaxis,h,s,r,1);
  nit=itercg;
  for(i=0;i<nx;i++){
    q1[i]=s[i]/Qp[i];
    q[i]=q1[i]/Qp[i];
  }
  for (i=0;i<nx;i++) z[i]=q[i];
  dq=dot(nx,s,q);
  for (i=0;i<nx;i++) z1[i]=q1[i];
  for(i=0;i<nx;i++) x1[i]=0.;       
  for (j=0;j<nit;j++){
    rstack(t,qaxis,h,z,Az,0);            
    alphanum=dq;
    alphaden=dot(ny,Az,Az);
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < tol ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      break;
    }
    alpha=alphanum/alphaden;
    alpha*=step;
    //fprintf(stderr,"j=%d,alpha=%e\n",j,alpha);         
    //// Update model u and residuals
    for(i=0;i<nx;i++) x[i]=x[i]+alpha*z[i];
    for(i=0;i<ny;i++) r[i]=r[i]-alpha*Az[i];  
    
    //resold=resid;
    rstack(t,qaxis,h,s,r,1);
    for(i=0;i<nx;i++){
      q1[i]=s[i]/Qp[i];
      q[i]=q1[i]/Qp[i];
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
  free1double(gcv);
  free1double(rho);
  free1double(eta);
  free1double(Az);
  free1double(r);
  free1double(z1);
  free1double(z);
  free1double(x1);
  free1double(s);
  free1double(q1);
  free1double(q); 
  free1double(Qp);       
  return;
}









