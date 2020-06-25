#include "su.h"
void conjzero(int conj, int add, int nx, float *x, int ny, float *y);

void contran(int conj,int sum, int nx, float *xx, int nb, float *bb,float *yy)
{
  // if sum==0 erases the output first
  int ny;
  int ib;
  int ix;
  ny=nx + nb -1;
  conjzero(conj,sum,nb,bb,ny,yy);
  if (conj == 0) 
    for(ib=0;ib<nb;ib++) 
      for(ix=0;ix<nx;ix++) 
	yy[ib+ix]=yy[ib+ix]+bb[ib]*xx[ix];
  else
    for(ib=0;ib<nb;ib++) 
      for(ix=0;ix<nx;ix++) 
	bb[ib]=bb[ib]+yy[ib+ix]*xx[ix];
  return;
}
