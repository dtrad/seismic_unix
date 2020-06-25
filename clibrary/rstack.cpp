#include "su.h"
#include "clibrarytd.h"
#include <math.h>
inline float sincin(double x)
{
  float y;
  if (x<1e-1) return(y=1.);
  y=(float) sin(PI*x)/(PI*x);
  return(y);
}
inline float sincin(double x);
void rstack(float *t,float *q,float *h,float *m,float *d,int conj)
{
  int i,k,ih,iq,itau;
  unsigned int j;
  float time,it,hxh, pxhxh;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  int tempit;
  
  
  //  Compute CMP gathers  when conj = 0
  //  Compute velocity gathers when conj = 1 
  
  if (conj==1) for (i=0;i<(nq*nt);i++) m[i]=0;
  if (conj==0) for (i=0;i<(nh*nt);i++) d[i]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    for (iq=0;iq<nq;iq++){
          pxhxh=q[iq]*hxh;
	  for (itau=0;itau<nt;itau++){
	    k=iq*nt+itau;
	    time=sqrt(pow(t[itau],2)+pxhxh);
	    it=time/dt;
            tempit=(int) floor(it);
            if ((it-tempit)>0.5) tempit+=1;  
            //fprintf(stderr,"it=%f,floor(it)=%d\n",it,tempit);
	    j=ih*nt+tempit;
	    if ((it!=nt)&&(j<ny)&&(k<nx)) {
	      if(conj==1) m[k]=m[k]+d[j];//*sincin(it-floor(it));
	      if(conj==0) d[j]=d[j]+m[k];//*sincin(it-floor(it));
	    }
	  }
    }	
  }
  return;
}

float sinc(double x)
{
  float y;
  if (x<1e-1) return(y=1);
  y=(float) sin(PI*x)/(PI*x);
  //fprintf(stderr,"y=%e,x=%e\n",y,x);
  return(y);
}














