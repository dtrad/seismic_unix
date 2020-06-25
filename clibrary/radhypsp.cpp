#include <math.h>
#include "su.h"


int radhypsp(float *t, float *h, float *q, float **vel, int nt, int nh, int nq){
  register unsigned int it;
  unsigned int ih,iq;
  unsigned int j;
  unsigned int nonzero; 
  float ftime;
  float dint;
  float time,hxh,pxhxh;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  float a1;
  float a2;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];


  // First we need to find out how big the sparse matrix will be
  j=0;
  for (it=0;it<nt;it++){
    for (ih=0;ih<nh;ih++){
      hxh=h[ih]*h[ih];
      ihxnt=ih*nt;
      for (iq=0;iq<nq;iq++){    
	pxhxh=hxh/vel[iq][it];
	iqxnt=iq*nt;
        time=sqrt(t[it]*t[it]+pxhxh);
	ftime=time/dt;
	itime=(int) floor(ftime+0.5);
	if (itime<nt) j++;
      }            
    }
  }
  nonzero=j;
  fprintf(stderr,"nonzero=%d\n",nonzero);
  return(nonzero);
}
