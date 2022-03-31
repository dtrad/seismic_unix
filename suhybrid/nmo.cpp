//#include <math.h>
#include "su.h"

void nmo(float *d,float *m,float *t,float h,float *vel,int invert,int nt,float dt,float smute) 
{
  register int it;
  float moveout;
  float t2;
  float vel2;
  float dt2=dt*dt;
  float *ttn;
  float *tnt;
  //  float smute=1.5; /* zero samples with NMO stretch exceeding smute */
  float osmute;	   /* 1/smute */
  int lmute=25;    /* length in samples of linear ramp for mute */
  int itmute=0;	   /* zero samples with indices less than itmute */
  float *atn;	   /* amplitude a(tn) for NMO */
  float *at;	   /* amplitude a(t) for inverse NMO */
  int sscale=0;	   /* if non-zero, apply NMO stretch scaling */
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
  /* compute inverse of stretch factor a(tn) */
  atn[0] = ttn[1]-ttn[0];
  for (it=1; it<nt; ++it)
    atn[it] = ttn[it]-ttn[it-1];
  //fprintf(stderr,"smute=%f\n",smute);
  /* determine index of first sample to survive mute */
  osmute = 1.0/smute;
  for (it=0,itmute=0; it<nt && atn[it]<osmute; ++it)
    itmute++;
  
  if (invert){ 
    yxtoxy(nt,1.0,0.0,&ttn[0],nt,1.0,0.0,-nt,nt,&tnt[0]);
    /* adjust mute time */
    itmute = (int) (1.0+ttn[itmute]);
    itmute = MIN(nt-2,itmute);
    
    /* compute a(t) */
    if (sscale) {
      for (it=itmute+1; it<nt; ++it)
	at[it] = tnt[it]-tnt[it-1];
      at[itmute] = at[itmute+1];
    }

    ints8r(nt,1.0,0,m,0.0,0.0,nt,tnt,d);
    /* apply mute */
    for (it=0; it<itmute; ++it)
      d[it] = 0.0;
			
    /* if specified, undo NMO stretch factor scaling */
    if (sscale)
      for (it=itmute; it<nt; ++it)
	d[it] *= at[it];


  }
  else{
    ints8r(nt,1.0,0,d,0.0,0.0,nt,ttn,m);
    /* apply mute */
    for (it=0; it<itmute; ++it)
      m[it] = 0.0;
    
    /* apply linear ramp */
    for (it=itmute; it<itmute+lmute && it<nt; ++it)
      m[it] *= (float)(it-itmute+1)/(float)lmute;
    
    /* if specified, scale by the NMO stretch factor */
    if (sscale)
      for (it=itmute; it<nt; ++it)
	m[it] *= atn[it];
  }

  free1float(at);
  free1float(atn);
  free1float(tnt);
  free1float(ttn);

  return;
}










