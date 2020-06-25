//#include <math.h>
#include "su.h"

void nmo(float **m, float *t, float *h, float *q, float **d, float *vel,int adj, int nt, int nh, int nq) 
{
  register int it;
  int ih,iq;
  unsigned int j;
  float ftime;
  float dint;
  float time,hxh,moveout;
  int itime;
  float a1;
  float a2;
  float t2;
  float vel2;

  float dt=t[1]-t[0];

  
  if (adj)  for (iq=0;iq<nq;iq++)  for (it=0;it<nt;it++) m[iq][it]=0;
  else for (ih=0;ih<nh;ih++)  for(it=0;it<nt;it++) d[ih][it]=0;
  
  for (it=0;it<nt;it++){
    vel2=vel[it]*vel[it];
    t2=t[it]*t[it];
    for (ih=0;ih<nh;ih++){
      hxh=h[ih]*h[ih];
      for (iq=0;iq<nq;iq++){
	moveout=hxh/vel2;
	time=sqrt(t2+pxhxh);
	ftime=time/dt;
	if (adj) ints8r(nt,1.0,0,&d[ih],0.0,0.0,1,&ftime,&m[iq][it]);
        else{
	  itime=(int) (ftime);
          a2=ftime-itime;
	  a1=1-a2;	  
	  if (itime<nt)  d[ih][itime]=a1*m[iqxnt+it];
	  if ((itime+1) < nt) d[ih][itime+1]+=a2*m[iqxnt+it];
	}
      }            
    }
  }
  
  return;
}










