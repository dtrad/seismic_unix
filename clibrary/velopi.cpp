#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "su.h"

void velopi(float *m, float *t, float *h, float *q, float *d) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  float *itime,*dint,*dtemp,time,hxh,pxhxh,a1,a2;
  int iqxnt,ihxnt,itint;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      pxhxh=hxh/(q[iq]*q[iq]);
      iqxnt=iq*nt;
      for (itau=0;itau<nt;itau++){
	time=sqrt(t[itau]*t[itau]+pxhxh);
	/*time=(t[itau]*t[itau])+pxhxh;*/
	itime[itau]=time/dt;
        dtemp[itau]=d[ihxnt+itau];
      }
      ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,itime,dint);
      for (itau=0;itau<nt;itau++){
 	k=iqxnt+itau;        
	m[k]=m[k]+dint[itau];            
      }
    }
  }
  free1float(dtemp);
  free1float(itime);
  free1float(dint);
  return;
}

void velopi(float *m, float *t, float *h, float *q, float *d, float *ww) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  float *itime,*dint,*dtemp,time,hxh,pxhxh,a1,a2;
  int iqxnt,ihxnt,itint;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");  
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      pxhxh=hxh/(q[iq]*q[iq]);
      iqxnt=iq*nt;
      for (itau=0;itau<nt;itau++){
	k=iqxnt+itau;
	//if (fabs(ww[k])>1e-2){
	time=sqrt(t[itau]*t[itau]+pxhxh);
	/*time=(t[itau]*t[itau])+pxhxh;*/
	itime[itau]=time/dt;
	dtemp[itau]=d[ihxnt+itau];
	//}
      }
      ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,itime,dint);
      for (itau=0;itau<nt;itau++){
	k=iqxnt+itau;        
	m[k]=m[k]+dint[itau];            	  
      }
    }
  } 
  free1float(dtemp);
  free1float(itime);
  free1float(dint);
  return;
}



void velopi( float *m, float *t, float *h, float *q, float *d, float *ww, float theta) 
{
  int i,k,ih,iq,itau,itint;;
  unsigned int j;
  float time,it,hxh,qxhxh,qxhxhxs,sintheta,a1,a2;
  int iqxnt,ihxnt;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
  theta=theta/180.*acos(-1.);
  sintheta=sin(theta)*sin(theta);    
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      qxhxh=hxh/(q[iq]*q[iq]);
      qxhxhxs=qxhxh*sintheta;
      iqxnt=iq*nt;
      for (itau=0;itau<nt;itau++){
	k=iqxnt+itau;
	if (fabs(ww[k])>1e-2){
	  time=sqrt(t[itau]*t[itau]+qxhxh-qxhxhxs);
	  /*time=(t[itau]*t[itau])+pxhxh;*/
	  it=time/dt;
	  itint=(int) floor(it);
	  j=ihxnt+itint;
	  a1=1-(it-itint);
	  a2=it-itint;
	  if ((it!=nt)&&(j<ny-1)&&(k<nx))
	    m[k]=m[k]+a1*d[j]+a2*d[j+1];/**sincin(it-floor(it));*/
	  //m[k]=m[k]+d[j];/**sincin(it-floor(it));*/	  
	}      
      }
    }
  }
  return;
}

















