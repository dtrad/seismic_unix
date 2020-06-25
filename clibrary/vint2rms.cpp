
#include "su.h"
 
void vint2rms(int inverse, float vminallow,float dt,float *vint,int nt,
	      float *vrms)
{
  int it, wide;
  float *vis, *sum;

  if ((vis=ealloc1float(nt))==NULL)
    err("***Sorry, space for vis could not be allocated\n");
  if ((sum=ealloc1float(nt))==NULL)
    err("***Sorry, space for vis could not be allocated\n");
  
  if (inverse==0){
    for(it=0;it<nt;it++)
      vis[it]=vint[it]*vint[it];
    for(sum[0]=0,it=1;it<nt;it++)
      sum[it]=sum[it-1]+vis[it]*dt;
    for(vrms[0]=vint[0],it=1;it<nt;it++)
      vrms[it]=sqrt(sum[it]/(it*dt));
  }
  else{
    for(it=0;it<nt;it++)
      sum[it]=it*dt*MAX(vrms[it]*vrms[it],vminallow*vminallow);
    vis[0]=vrms[0]*vrms[0];
    for(it=1;it<nt;it++)
      vis[it]=(sum[it]-sum[it-1])/dt;
    wide=2;
    while(1){
      vmin=vis[0]; for(it=0;it<nt;it++)
    }
  }
  free1float(vis);
  free1float(sum);
  return;
}
