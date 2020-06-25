#include "dan.h"
void gcvfunction(float *gcv, float rho, float normdata, int nh, int nq, int iter)
{
  float in;
  int i;
  in=iter-1;
  for (i=1;i<=in;i++){
    num=(nh-(i-1))*(nh-(i-1));
    gcv[i]=(rho[i]*rho[i])/num;
  }
  

  return();

}
 
