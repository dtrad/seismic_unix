#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void radonopi( float *m, float *t, float *h, float *q, float *d, float *ww) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  float time,it,hxh,pxhxh;
  int iqxnt,ihxnt;
  extern int nt,nh,nq, nx, ny;
  extern float dt,dh,dq;
    
  for (i=0;i<(nq*nt);i++) m[i]=0;

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
	ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
		if (fabs(ww[iq])>1e-2){
			pxhxh=hxh*q[iq];
			iqxnt=iq*nt;
			for (itau=0;itau<nt;itau++){
				k=iqxnt+itau;
			    time=sqrt(pow(t[itau],2)+pxhxh);
                /*time=(t[itau]*t[itau])+pxhxh;*/
				it=time/dt;
				j=ihxnt+(int) floor(it);
				if ((it!=nt)&&(j<ny)&&(k<nx)) 
					m[k]=m[k]+d[j];/**sincin(it-floor(it));*/
		
			}
		}
	  
	}
  }
  return;
}






