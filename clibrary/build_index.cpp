#include <math.h>
#include "su.h"

void build_index(float *t, float *h, float *q, float **vgrid, int nt, int nh, 
		  int nq,unsigned short **index, int nonzero, float t0)
{
  unsigned int it;
  unsigned int ih,iq;
  float time,ftime,hxh,moveout,slowness2;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];
  int it0=(int) (t0/dt+0.5);


  for (j=0;j<2;j++) for (it=0;it<nonzero;it++) index[j][it]=0;

  for (it=it0;it<nt;it++){
      for (iq=0;iq<nq;iq++){    
	slowness2=1./(vgrid[iq][it]);
	iqxnt=iq*nt;
	for (ih=0;ih<nh;ih++){
	  hxh=h[ih]*h[ih];
	  ihxnt=ih*nt;
	  moveout=hxh*slowness2;
	  time=sqrt(t[it]*t[it]+moveout);
	  ftime=time/dt;
	  itime=(int) floor(ftime+0.5);
	  if (itime<nt){
	    index[0][ih*nq*nt+iqxnt+it]=ihxnt+itime;
	    index[1][ih*nq*nt+iqxnt+it]=iqxnt+it;
	  }
	}            
      }
  }
  return;
}
