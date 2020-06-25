#include "su.h"
#include "clibrarytd.h"
#include <math.h>
inline double sincin(double x)
{
  double y;
  if (x<1e-1) return(y=1.);
  y=(double) sin(PI*x)/(PI*x);
  return(y);
}
inline double sincin(double x);
void rstack(double *t,double *q,double *h,double *m,double *d,int conj)
{
  int i,k,ih,iq,itau;
  unsigned int j;
  double time,it;
  extern int nt,nh,nq, nx, ny;
  extern double dt,dh,dq;
  
  
  
  //  Compute CMP gathers  when conj = 0
  //  Compute velocity gathers when conj = 1 
  
  if (conj==1) for (i=0;i<(nq*nt);i++) m[i]=0;
  if (conj==0) for (i=0;i<(nh*nt);i++) d[i]=0;
  
  for (ih=0;ih<nh;ih++){
    for (iq=0;iq<nq;iq++){
	  for (itau=0;itau<nt;itau++){
	    k=iq*nt+itau;
	    time=sqrt(pow(t[itau],2)+pow(h[ih],2)*q[iq]);
	    it=time/dt;
	    j=ih*nt+(int) floor(it);
	    if ((it!=nt)&&(j<ny)&&(k<nx)) {
	      if(conj==1) {
                //if (ih>0) dh=h[ih]-h[ih-1];
                //else dh=h[ih];
		m[k]=m[k]+d[j]*sincin(it-floor(it));
	      }
	      if(conj==0) d[j]=d[j]+m[k]*sincin(it-floor(it));
	    }
	  }
    }	
  }
  return;
}

double sinc(double x)
{
  double y;
  if (x<1e-1) return(y=1);
  y=(double) sin(PI*x)/(PI*x);
  //fprintf(stderr,"y=%e,x=%e\n",y,x);
  return(y);
}
















