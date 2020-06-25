#include "su.h"
#include "segy.h"
#include "Complex.h"
void xequaly(complex *x,complex *y,int n)
{
     memcpy((void *) x,(const void*) y,n*sizeof(complex));
     return;
}
void xequaly(float *x,float *y,int n)
{
     memcpy((void *) x,(const void*) y,n*FSIZE);
     return;
}
void AequalB(float **A,float **B,int n1, int n2)
{
  int i2;
  for (i2=0;i2<n2;i2++)
     memcpy((void *) A[i2],(const void*) B[i2],n1*FSIZE);
  return;
}

void AequalB(complex **A, complex **B,int n1, int n2)
{
  int i1,i2;
  for (i2=0;i2<n2;i2++)
    for (i1=0;i1<n1;i1++)
      A[i2][i1]=B[i2][i1];
  return;
}
