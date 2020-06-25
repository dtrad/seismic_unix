#include "radonfk.h"

#define testpad 1 /* Pad with extra zeros in offset */

void stoltzop2(float **data, float ***model, int nt, int nx, int nv, float *t, 
	       float *x, float *vel, int adj)
{
  int iv;

  if (adj)
    for (iv=0;iv<nv;iv++)
      stoltzop2(data,model[iv],nt,nx,t,x,vel[iv],adj);
  else{
    float **datatemp=ealloc2float(nt,nx);
    zero_array(data,nx,nt);
    for (iv=0;iv<nv;iv++){ 
      zero_array(datatemp,nx,nt);
      stoltzop2(datatemp,model[iv],nt,nx,t,x,vel[iv],adj);
      Aplus_equalB(data,datatemp,nt,nx);
    } 
    free2float(datatemp);
  }

  return;

}







