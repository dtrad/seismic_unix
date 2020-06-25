#include "dan.h"
#include "su.h"

/* SImple mute routine. 
   The parmute parameter defines the boundary between muted and non muted spaces
   A taper is applied for smoothness but it may be too smooth in many cases,
   By default is off.
   the minimum time to mute follows a slope give by it0-(iq-nmute)*slope
   This slope is a function of nt/nq. In general values between 2-5 must work fine */

#define VERBOSE 1
#define SLOPE 2
#define TAPER 0

int taper(float **data, int nt, int nh, int ntaper,int flag);

void muting(float **model, int nq, int nt, float *q, float *t, float parmute, float t0, int side)
{

  // side 2 --> right mute ; side 1 left mute
  int iq, it;
  int nmute=0;
  float dt=t[1]-t[0];
  int it0=(int) (t0/dt); // define a minimum time for mute2
  int slope=SLOPE; // Defines a slope in the muting

  iq=0; while(q[iq]<parmute) iq++; 
  if (side==2) nmute=nq-iq;
  else if (side==1) nmute=iq;

  if (VERBOSE) fprintf(stderr,"MUTING at nmute=%d \n",nmute);
  
  if (TAPER) taper(model,nt,nq,nmute,side);
  else{
    if (side==1)  
      for (iq=0;iq<MIN(nq,nmute);iq++) for (it=0;it<nt;it++) model[iq][it]=0; 
    else   
      for (iq=MIN(nq,nmute);iq<nq;iq++) for (it=0;it<nt;it++) model[iq][it]=0; 
  }

  for (iq=0;iq<nq;iq++) 
    for (it=0;it<(it0-(iq-nmute)*slope);it++) 
      model[iq][it]=0; 
  
  return;
}



 
 


