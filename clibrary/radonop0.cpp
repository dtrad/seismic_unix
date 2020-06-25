#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void radonop(float *m, float *t, float *h, float *q, float *d) 
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
			time=sqrt(pow(t[itau],2)+qxhxh);
			it=time/dt;
			j=ih*nt+(int) floor(it);
			if ((it!=nt)&&(j<ny)&&(k<nx)) {
			d[j]=d[j]+m[k];/**sinc(it-floor(it));*/
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


