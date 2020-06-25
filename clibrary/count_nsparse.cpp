#include <math.h>
#include "su.h"


int count_nsparse(float *t, float *h, float *q, float **vgrid, int nt, int nh, 
		  int nq, float t0)
{
  unsigned int it;
  unsigned int ih,iq;
  unsigned int j;
  unsigned int nonzero; 
  float time,ftime,hxh,moveout,slowness2;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];
  int it0=(int) (t0/dt+0.5);
  // First we need to find out how big the sparse matrix will be
  j=0;

  for (it=it0;it<nt;it++){
    for (iq=0;iq<nq;iq++){
      slowness2=1./(vgrid[iq][it]);
      iqxnt=iq*nt;
      for (ih=0;ih<nh;ih++){
	hxh=h[ih]*h[ih];
	ihxnt=ih*nt;
	moveout=hxh*slowness2;
        time=2*sqrt(t[it]*t[it]/4.+moveout);
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
