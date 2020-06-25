#include "su.h"
#include "Complex.h"

void mask(float *Wm, complex *m, int n, float threshold)
  // Function to compute a mask 
{
  int i;
  // Mask computation
  if ((threshold>=0.7)||(threshold<=0.01)) threshold=0.3;
  /*Mask is given by losigm function*/
  for (i=0;i<nx;i++) Wm[i]=1+1./(1+exp(100*(abs(m[i])-threshold)+0.5));

  return;
}
