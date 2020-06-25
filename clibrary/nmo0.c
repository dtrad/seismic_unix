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

  if (invert){
    for(it=0;it<nt;it++) d[it]=0;
    if ((ttn=ealloc1float(nt))==NULL)
      fprintf(stderr,"***Sorry, space for ttn could not be allocated\n");
    if ((tnt=ealloc1float(nt))==NULL)
      fprintf(stderr,"***Sorry, space for tnt could not be allocated\n");
    for (it=0;it<nt;it++){
      vel2=vel[it]*vel[it];
      t2=t[it]*t[it];
      moveout=h*h/vel2;
      time=sqrt(t2/dt2+moveout/dt2);
      ttn[it]=time;
    }
    yxtoxy(nt,1.0,0.0,&ttn[0],nt,1.0,0.0,-nt,nt,&tnt[0]);
    ints8r(nt,1.0,0,m,0.0,0.0,nt,tnt,d);

    free1float(tnt);
    free1float(ttn);
  }
  else{ 
    for (it=0;it<nt;it++) m[it]=0;    
    for (it=0;it<nt;it++){
      vel2=vel[it]*vel[it];
      t2=t[it]*t[it];
      moveout=h*h/vel2;
      time=sqrt(t2/dt2+moveout/dt2);
      ints8r(nt,1.0,0,d,0.0,0.0,1,&time,&m[it]);
    }
  }
  return;
}










