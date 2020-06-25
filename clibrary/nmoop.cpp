#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "su.h"


void nmoop(float *m, float *t, float *h, float *q, float *d, float *invnmo) 
{
  int i,k,ih,iq,itau,itint,ihxnt;
  unsigned int j;
  float time,it,hxh,qxhxh,a1,a2;
  extern int nt,nq,nh,nx,ny;
  extern float dt,dh,dq;
    
  for (i=0;i<(nh*nt);i++) d[i]=0;

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      for (itau=0;itau<nt;itau++){
	   //qxhxh=(q[iq]-invnmo[itau])*hxh;
	qxhxh=0;
	k=iq*nt+itau;
	if (fabs(m[k])>1e-5){
	  time=sqrt(t[itau]*t[itau]+qxhxh);
	  it=time/dt;
	  itint=(int) floor(it);
	  j=ihxnt+itint;
          a1=1-(it-itint);
          a2=it-itint;
      	  if ((it!=nt)&&(j<ny-1)&&(k<nx)) {
	    d[j]=d[j]+a1*m[k];
            d[j+1]=d[j+1]+a2*m[k];
 	  }
	}
      }
      
    }
  }
  return;
}


void nmoopold(float *m, float *t, float *h, float *q, float *d, float *invnmo) 
{
  int i,k,ih,iq,itau,itint,ihxnt,ktime;
  unsigned int j;
  float time,it,hxh,qxhxh;
  extern int nt,nq,nh,nx,ny;
  extern float dt,dh,dq;
  float *mint,*mtemp,*itime,t2;


  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((mint=alloc1float(nt))==NULL)
    err("cannot allocate memory for mint \n");
  if ((mtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for mtemp \n");

  for (i=0;i<(nh*nt);i++) d[i]=0;
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      for (i=0;i<nt;i++) mint[i]=0;
      for (i=0;i<nt;i++) mtemp[i]=0;
      for (ktime=0,itau=0;itau<nt;itau++){
	qxhxh=(q[iq]-invnmo[itau])*hxh;
        t2=t[itau]*t[itau];
       	if (t2>=qxhxh){ 
            time=sqrt(t2-qxhxh);	
	    itime[itau]=time/dt;
            ktime+=1;
	}
	mtemp[itau]=m[iq*nt+itau];
	
      }
     
      ints8r(nt,1.,0.,mtemp,0.0,0.0,ktime,itime,mint);
      for (itau=0;itau<ktime;itau++) {
	j=ihxnt+itau;
	d[j]=d[j]+mint[itau];
      }      
    }    
  }

  free1float(mint);
  free1float(mtemp);
  free1float(itime);

  return;
}










