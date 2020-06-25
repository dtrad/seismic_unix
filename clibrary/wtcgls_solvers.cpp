#include "su.h"
#include "clibrary.h"
#include <math.h>

float wtcgls(complex *b,complex **L, complex *x,float *Wm,
	    float *Wd,int nh,int nq, float tol, float step, int itercg)
{
  /* This function solves the system of equations 
     (FH WdT Wd FH + WmT Wm ) m = FH WdT Wd d 
     Notice that LH=FH WdT and L= Wd F
     so that Wd only needs to be used the first time 
  */
 
  float normb,dq,dq2,nit,beta, alpha, alphanum, alphaden;
  int j,in,num;
  complex czero;

  register int i;
  complex *q;
  complex *q1;
  complex *s;
  complex *x1;
  complex *z;
  complex *z1;
  complex *r;
  complex *Az;
  complex *xold;
  float *eta;
  float *rho;
  float *gcv;
  float J;
  complex *r2;
  complex *xtemp;

  q=ealloc1complex(nq); 
  q1=ealloc1complex(nq);
  s=ealloc1complex(nq);
  x1=ealloc1complex(nq);
  z=ealloc1complex(nq);
  z1=ealloc1complex(nq);
  r=ealloc1complex(nh);
  Az=ealloc1complex(nh);
  eta=ealloc1float(nq);
  rho=ealloc1float(nq);
  gcv=ealloc1float(nq);
  r2=ealloc1complex(nh);
  xtemp=ealloc1complex(nq);
  xold=ealloc1complex(nq);

  czero.r=czero.i=0;
  for (i=0;i<nq;i++) x[i]=czero;
  normb=sqrt(rcdot(nh,b,b));
  //xequaly(r,b,nh);
  for (i=0;i<nh;i++) r[i]=Wd[i]*b[i];
  for (i=0;i<nh;i++) r2[i]=r[i]*Wd[i];  
  Atimesx(r2,L,s,nh,nq,1);
  
  nit=MIN(itercg,nq);
  for(i=0;i<nq;i++){
    q1[i]=s[i]/Wm[i];
    q[i]=q1[i]/Wm[i];
  }
  xequaly(z,q,nq);
  dq=rcdot(nq,s,q);
  xequaly(z1,q1,nq);
  for(i=0;i<nq;i++) x1[i]=czero;       
  for (j=0;j<nit;j++){
    Atimesx(Az,L,z,nh,nq,0);            
    for (i=0;i<nh;i++) Az[i]*=Wd[i];  
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

    for(i=0;i<nh;i++) r[i]-=alpha*Az[i];
    for(i=0;i<nq;i++) xold[i]=x[i];    
    for(i=0;i<nq;i++) x[i]+=alpha*z[i];    

    for(i=0;i<nh;i++) r2[i]=r[i]*Wd[i];
    Atimesx(r2,L,s,nh,nq,1);

    for(i=0;i<nq;i++){
      q1[i]=s[i]/Wm[i];
      q[i]=q1[i]/Wm[i];
    }
    dq2=rcdot(nq,s,q);
    beta=dq2/dq;
    dq=dq2;
    for (i=0;i<nq;i++) z[i]=q[i]+beta*z[i];
    rho[j] = sqrt(rcdot(nq,s,s))/normb;
    //fprintf(stderr,"rho[%d]=%e\n",j,rho[j]);
    for (i=0;i<nq;i++) xtemp[i]=x[i]/Wm[i];
    J=rcdot(nh,r2,r2)+rcdot(nq,xtemp,xtemp);
    for (i=0;i<nq;i++) {
      x1[i]=x1[i]+alpha*z1[i]; 
      z1[i]=q1[i]+beta*z1[i];
    }
    eta[j]=sqrt(rcdot(nq,x1,x1));
    if ((tol==0) && (j>1)){ // GCV criteria
       in = j;
       for (i=0;i<in;i++){
       num=(nq-i)*(nq-i); 
       gcv[i]=(rho[i]*rho[i])/num;
       }      
       if ((gcv[j-2]<gcv[j-1])&&(gcv[j-1]<gcv[j])) { 
         if (1) fprintf(stderr,"GCV Criteria, iteration %d\n",j-1);
	 for(i=0;i<nq;i++) x[i]=xold[i];    
         nit = j-1;
         break;
       } 
       
       else if ((tol!=0) && (rho[j] < tol)){ 
        if (1) fprintf(stderr,"Convergence have been acheived at iteration # %d\n",j);
        break;
       }
    }          
  }
  if (0) fprintf(stderr,"j=%d\n",j);

  free1complex(xold);
  free1complex(xtemp);
  free1complex(r2);
  free1float(gcv);
  free1float(rho);
  free1float(eta);
  free1complex(Az);
  free1complex(r);
  free1complex(z1);
  free1complex(z);
  free1complex(x1);
  free1complex(s);
  free1complex(q1);
  free1complex(q);

 

  return(J);
}















