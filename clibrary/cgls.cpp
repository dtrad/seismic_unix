#include "su.h"
#include "Complex.h"
#include "clibrary.h"
#include <math.h>
void cgls(complex *d,complex **L, complex **LH,complex *u,complex *Qp,
	  int nh,int nq,complex *q,complex *r,complex *raux,complex *s,
	  complex *p, float eps1,float eps2, float eps)
{
  float bb, resid, resold, beta, betanum, betaden, alfa, alfanum, alfaden;
  int k,i,j;
  extern int itercg;
  complex czero;
  bb=rcdot(nh,d,d);
  Atimesx(u,LH,d,nq,nh);
  Atimesx(s,L,u,nh,nq);
  xminusy(s,d,s,nh);
  //for (i=0;i<nq;i++) u[i]=czero;
  //for (i=0;i<nh;i++) s[i]=czero;
  //for (i=0;i<nq;i++) p[i]=czero;  
  //xequaly(s,d,nh);
  Atimesx(p,LH,s,nq,nh);
  for (i=0;i<nq;i++) p[i]/=Qp[i];  
  xequaly(r,p,nq); 
  resid=rcdot(nq,r,r);
  k=0;
 
  while ((k!=nq)&&(sqrt(resid)>(eps*bb))&&(k<itercg)){
    k++;
    Atimesx(q,L,p,nh,nq);            
    alfanum=rcdot(nq,r,r);
    
    alfaden=rcdot(nh,q,q);
    if (alfaden < 0.) err("alfaden=%e\n",alfaden);
    if (alfaden < eps ){ 
      fprintf(stderr,"alfanum=%e,alfaden=%e,k=%d,eps=%e\n",
	      alfanum,alfaden,k,eps);
      break;
    }
    alfa=alfanum/alfaden;
    //fprintf(stderr,"alfanum=%e,alfaden=%e,k=%d\n",
    //	      alfanum,alfaden,k); 
    //fprintf(stderr,"alfa=%e\n",alfa);
    //// Update model u and residuals
    for(i=0;i<nq;i++) u[i]+=(alfa*p[i]);
    for(i=0;i<nh;i++) s[i]-=(alfa*q[i]);  
    
    resold=resid;
    Atimesx(r,LH,s,nq,nh);
    for (i=0;i<nq;i++) r[i]/=Qp[i];  
    //xtimesy(raux,Qp,u,nq);
    //xminusy(r,r,raux,nq);
    resid=rcdot(nq,r,r);

    betaden=resold;
    betanum=resid;
    if (betaden < eps){
      fprintf(stderr,"betanum=%e,betaden=%e,k=%d,eps=%e\n",
	      betanum,betaden,k,eps); 
      break;
    }
    beta=betanum/betaden;
    //fprintf(stderr,"beta=%e\n",beta);
    for(i=0;i<nq;i++) p[i]=r[i]+beta*p[i];           
  }
  fprintf(stderr,"k=%d\n",k);        
return;
}








