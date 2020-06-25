#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
float sinc(float arg);

void radonop(float *m, float *t, float *h, float *q, float *d) 
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
      qxhxh=q[iq]*hxh;
      for (itau=0;itau<nt;itau++){
	k=iq*nt+itau;
	if (fabs(m[k])>1e-2){
	  time=sqrt(t[itau]*t[itau]+qxhxh);
	  it=time/dt;
	  itint=(int) floor(it);
	  j=ihxnt+itint;
          a1=1-(it-itint);
          a2=it-itint;
      	  if ((it!=nt)&&(j<ny-1)&&(k<nx)) {
	    d[j]=d[j]+a1*m[k];//*sinc(it-itint);
            d[j+1]=d[j+1]+a2*m[k];
 	  }
	}
      }
      
    }
  }
  return;
}
void radonop(float *m, float *t, float *h, float *q, float *d, float theta) 
{
  int i,k,ih,iq,itau,itint,ihxnt;
  unsigned int j;
  float time,it,hxh,qxhxh,qxhxhxs,a1,a2,sintheta;
  extern int nt,nq,nh,nx,ny;
  extern float dt,dh,dq;
  theta=theta/180.*acos(-1.);
  sintheta=sin(theta)*sin(theta); 
    
  for (i=0;i<(nh*nt);i++) d[i]=0;

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      qxhxh=q[iq]*hxh;
      qxhxhxs=qxhxh*sintheta;
      for (itau=0;itau<nt;itau++){
	k=iq*nt+itau;
	if (fabs(m[k])>1e-2){
	  time=sqrt(t[itau]*t[itau]+qxhxh-qxhxhxs);
	  it=time/dt;
	  itint=(int) floor(it);
	  j=ihxnt+itint;
          a1=1-(it-itint);
          a2=it-itint;
      	  if ((it!=nt)&&(j<ny-1)&&(k<nx)) {
	    d[j]=d[j]+a1*m[k];//*sinc(it-itint);
            d[j+1]=d[j+1]+a2*m[k];
 	  }
	}
      }
      
    }
  }
  return;
}



void radonop(float *m, float *t, float *h, float *q, float *d, char flag) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  float time,it,hxh,qxhxh;
  extern int nt,nq,nh,nx,ny;
  extern float dt,dh,dq;
    
  for (i=0;i<(nh*nt);i++) d[i]=0;

  for (ih=0;ih<nh;ih++){
	hxh=h[ih]*h[ih];
    for (iq=0;iq<nq;iq++){
      qxhxh=q[iq]*hxh;
	  for (itau=0;itau<nt;itau++){
	    k=iq*nt+itau;
		if (fabs(m[k])>1e-2){
			time=sqrt(t[itau]*t[itau]+qxhxh);
			it=time/dt;
			j=ih*nt+(int) floor(it);
			if ((it!=nt)&&(j<ny)&&(k<nx)) {
			d[j]=d[j]+m[k];
			}
		}
	  }
	  
    }
  }
  return;
}

void radonop_par(float *m, float *t, float *h, float *q, float *d) 
{
  int i,k,ih,iq,itau,itint,ihxnt;
  unsigned int j;
  float time,it,qxhxh,a1,a2;
  extern int nt,nq,nh,nx,ny;
  extern float dt,dh,dq;
    
  for (i=0;i<(nh*nt);i++) d[i]=0;

  for (ih=0;ih<nh;ih++){
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      qxhxh=q[iq]*h[ih]*h[ih];
      for (itau=0;itau<nt;itau++){
	k=iq*nt+itau;
	if (fabs(m[k])>1e-2){
	  time=t[itau]+qxhxh;
	  it=time/dt;
	  itint=(int) floor(it);
	  j=ihxnt+itint;
          a1=1-(it-itint);
          a2=it-itint;
      	  if ((it!=nt)&&(j<ny-1)&&(k<nx)) {
	    d[j]=d[j]+a1*m[k];//*sinc(it-itint);
            d[j+1]=d[j+1]+a2*m[k];
 	  }
	}
      }
      
    }
  }
  return;
}

void radonop_lin(float *m, float *t, float *h, float *q, float *d) 
{
  int i,k,ih,iq,itau,itint,ihxnt;
  unsigned int j;
  float time,it,qxh,a1,a2;
  extern int nt,nq,nh,nx,ny;
  extern float dt,dh,dq;
    
  for (i=0;i<(nh*nt);i++) d[i]=0;

  for (ih=0;ih<nh;ih++){
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      qxh=q[iq]*h[ih];
      for (itau=0;itau<nt;itau++){
	k=iq*nt+itau;
	if (fabs(m[k])>1e-2){
	  time=t[itau]+qxh;
	  it=time/dt;
	  itint=(int) floor(it);
	  j=ihxnt+itint;
          a1=1-(it-itint);
          a2=it-itint;
      	  if ((it!=nt)&&(j<ny-1)&&(k<nx)) {
	    d[j]=d[j]+a1*m[k];//*sinc(it-itint);
            d[j+1]=d[j+1]+a2*m[k];
 	  }
	}
      }
      
    }
  }
  return;
}


float sinc(float arg)
{
	float y;
	float pi = 3.1415926535;
	if (arg<1e-1) return(y=1.);
	y=sin(pi*arg)/(pi*arg);
	return(y);
}











