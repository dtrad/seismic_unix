#include "su.h"
#include "clibrary.h"
#include <math.h>

void wtcgls(complex *b,complex **L, complex **LH,complex *x,complex *Wm,
	  float *Wd,int nh,int nq,complex *q,complex *q1,complex *s,
	  complex *x1,complex *z,complex *z1,complex *r,complex *Az,
	  float *eta,float *rho, float *gcv, float tol)
{
  /* This function solves the system of equations 
     (FH WdT Wd FH + WmT Wm ) m = FH WdT Wd d 
     Notice that LH=FH WdT and L= Wd F
     so that Wd only needs to be used the first time 
  */
 
  float normb,dq,dq2,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j,in,num;
  extern int itercg;
  complex czero;
  float extern step;
  czero.r=czero.i=0;
  for (i=0;i<nq;i++) x[i]=czero;
  normb=sqrt(rcdot(nh,b,b));
  //xequaly(r,b,nh);
  for (i=0;i<nh;i++) r[i]=b[i]*Wd[i];  
  Atimesx(s,LH,r,nq,nh);
  
  nit=itercg;
  for(i=0;i<nq;i++){
    q1[i]=s[i]/Wm[i];
    q[i]=q1[i]/Wm[i];
  }
  xequaly(z,q,nq);
  dq=rcdot(nq,s,q);
  xequaly(z1,q1,nq);
  for(i=0;i<nq;i++) x1[i]=czero;       
  for (j=0;j<nit;j++){
    Atimesx(Az,L,z,nh,nq);            
    alphanum=dq;
    alphaden=rcdot(nh,Az,Az);
    if (alphaden < 0.) err("alphaden=%e\n",alphaden);
    if (alphaden < 1.e-7 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e\n",
	      alphanum,alphaden);
      break;
    }
    alpha=alphanum/alphaden;
    alpha*=step;
         
    //// Update model u and residuals

    for(i=0;i<nh;i++) r[i]=r[i]-alpha*Az[i];
    for(i=0;i<nq;i++) x[i]=x[i]+alpha*z[i];    

    Atimesx(s,LH,r,nq,nh);
    for(i=0;i<nq;i++){
      q1[i]=s[i]/Wm[i];
      q[i]=q1[i]/Wm[i];
    }
    dq2=rcdot(nq,s,q);
    beta=dq2/dq;
    dq=dq2;
    for (i=0;i<nq;i++) z[i]=q[i]+beta*z[i];
    rho[j] = sqrt(rcdot(nh,r,r))/normb;
    //fprintf(stderr,"rho[%d]=%e\n",j,rho[j]);      
    for (i=0;i<nq;i++) {
      x1[i]=x1[i]+alpha*z1[i]; 
      z1[i]=q1[i]+beta*z1[i];
    }
    eta[j]=sqrt(rcdot(nq,x1,x1));
    if ((tol==0) && (j>2)){ // GCV criteria
       in = j-1;
       for (i=1;i<=in;i++){
       num=(nh-(i-1))*(nh-(i-1)); 
       gcv[i]=(rho[i]*rho[i])/num;
       }      
       if (gcv[j-2]<gcv[j-1]){ 
         fprintf(stderr,"GCV Criteria, iteration %d\n",j-1);
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
return;
}















