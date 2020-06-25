#include "su.h"
void conjzero(int conj, int add, int nx, float *x, int ny, float *y) 
{
  int ix;
  int iy;
  if (add == 0)
    if(conj == 0)
      for (iy=0;iy<ny;iy++)
	y[iy]=0.;
    else
      for (ix=0;ix<nx;ix++)
	x[ix]=0;
  return;
}
