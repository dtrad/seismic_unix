#include "su.h"
#include "clibrary.h"

int taper(float **data, int nt, int nh, int ntaper,int flag);
void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,
	 float smute);

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int  nh, int nq, int rtmethod);


/*
radon00inv
Input: is a Radon model 
Output: is the data gather

Daniel Trad - September 20- 2000
*/


void radon00inv(float *h, int nh, float **data ,float *t, int nt, float dt, float *vel, float **model, float *q, int nq, float smute,float nmofactor, int pseudohyp, float depth, float fmax,int rtmethod)
{
  int i, it, iq, nqt, ih;
  float *dtemp;
  float dq;
  float *posh;  // Pseudo offset for pseudohyperbolic for h[]
  float *htemp;
  
  if ((dtemp=alloc1float(nt))==NULL)
     err("cannot allocate memory for dtemp\n");
  
  //////////////////////////////////////////////////////////////////////
  // The following lines are for pseudo hyperbolic RT
  // All methods are the same with the only difference that the offsets are
  // transformed to pseudo offsets.

  if (pseudohyp){
    fprintf(stderr,"Using Pseudo hyperbolic RT#######################\n");
    rtmethod=1;
    posh=ealloc1float(nh);  
    for (i=0;i<nh;i++) posh[i]=sqrt(h[i]*h[i]+depth*depth)-depth;
  }

  //////////////////////////////////////////////////////////////////////

  // For pseudo hyperbolic replace h with pos 
  // keeping the original offset h and h2 for the final process
  if (pseudohyp){
    htemp=ealloc1float(nh);
    for (ih=0;ih<nh;ih++) htemp[ih]=h[ih];
    for (ih=0;ih<nh;ih++) h[ih]=posh[ih];
  }
  
  if (1) for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) data[ih][it]=0;

  hrrti(data,h,dt,model,q,fmax,nt,nh,nq,rtmethod);
  
  // Before applying inverse nmo recover the offset and free auxiliar arrays
  if (pseudohyp){
    for (ih=0;ih<nh;ih++) h[ih]=htemp[ih];
    free1float(htemp);
    free1float(posh);
  }     

  /////////////////////////////
  if (1)
    for (ih=0;ih<nh;ih++){
      for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
      nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
    }

  free1float(dtemp);

  return;

}

























