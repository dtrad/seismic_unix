//#include <math.h>
#include "su.h"

void nmo(float *d,float *m,float *t,float h,float *vel,int invert,int nt,float dt) 
{
  register int it;
  int ih;
  unsigned int j;
  float ftime;
  float time,hxh,moveout;
  int itime;
  float a1;
  float a2;
  float t2;
  float vel2;
  float dt2=dt*dt;
  float *ttn;
  float *tnt;
  float smute=1.5;
  int lmute=25;
  float *atn;	/* amplitude a(tn) for NMO */
  float *at;	/* amplitude a(t) for inverse NMO */

  ttn=ealloc1float(nt);
  at = ealloc1float(nt);
  atn = ealloc1float(nt);
  tnt=ealloc1float(nt);

  if (invert) for(it=0;it<nt;it++) d[it]=0;
  else for (it=0;it<nt;it++) m[it]=0;

  for (it=0;it<nt;it++){
    vel2=vel[it]*vel[it];
    t2=t[it]*t[it];
    moveout=h*h/vel2;
    ttn[it]=sqrt(t2/dt2+moveout/dt2);
  }
  if (invert){ 
    yxtoxy(nt,1.0,0.0,&ttn[0],nt,1.0,0.0,-nt,nt,&tnt[0]);
    ints8r(nt,1.0,0,m,0.0,0.0,nt,tnt,d);
  }
  else ints8r(nt,1.0,0,d,0.0,0.0,nt,ttn,m);


  free1float(at);
  free1float(atn);
  free1float(tnt);
  free1float(ttn);

  return;
}










