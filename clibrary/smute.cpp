#include "su.h"

void smute_trace(float *d, float *t,float h,float *vel,int nt,float dt,float smute);
 
void smute_gather(float **d, int nt, int nh, float *t, float *h, float *vel, float smute)
  /* Applying a mute to a data gather based on the stretch produced by nmo correction */
{
  int ih;
  float dt=t[1]-t[0];
  for (ih=0;ih<nh;ih++) smute_trace(d[ih],t,h[ih],vel,nt,dt,smute);
  return;
}


void smute_trace(float *d, float *t,float h,float *vel,int nt,float dt,float smute) 
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
  
  /* determine index of first sample to survive mute */
  osmute = 1.0/smute;
  for (it=0,itmute=0; it<nt && atn[it]<osmute; ++it)
    itmute++;
  

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

  /* apply mute */
  for (it=0; it<itmute; ++it)
    d[it] = 0.0;
  
  /* if specified, undo NMO stretch factor scaling */
  if (sscale)
    for (it=itmute; it<nt; ++it)
      d[it] *= at[it];

  
  if (0){
    /* apply mute */
    for (it=0; it<itmute; ++it)
      d[it] = 0.0;
    
    /* apply linear ramp */
    for (it=itmute; it<itmute+lmute && it<nt; ++it)
      d[it] *= (float)(it-itmute+1)/(float)lmute;
    
    /* if specified, scale by the NMO stretch factor */
    if (sscale)
      for (it=itmute; it<nt; ++it)
	d[it] *= atn[it];
  }

  free1float(at);
  free1float(atn);
  free1float(tnt);
  free1float(ttn);

  return;
}

