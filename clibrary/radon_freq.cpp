#include "su.h"
#include "Complex.h"

void radon_freq(complex *x, complex *d, complex **A,int adj,int n,int m)
{
  int i,j;
  //fprintf(stderr,">>>>>>>>>>>>>>>\n");
  if (adj) for (i=0;i<m;i++) x[i].r=x[i].i=0;
  else  for (i=0;i<n;i++) d[i].r=d[i].i=0;

  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      if (adj) x[j]+=conjg(A[i][j])*d[i];
      else d[i]+=A[i][j]*x[j];

  return;
}












