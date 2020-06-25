#include "su.h"

void scale(float factor, int n, float *data)
{
  for (int i=0;i<n;i++) data[i]*=factor;
  return;
}

