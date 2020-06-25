#include "su.h"

int findpower2(int nt)
{
  int nn;
  for (nn=1;nn<nt;nn<<=1) ;
  return(nn);
}

void padding_zeros(float **d, int n2, int n1, int n1z)
{
  int i2, i1;

  if (n1 > n1z) err("n1 > n1z\n");
 
  for (i2=0;i2<n2;i2++) 
    for (i1=n1;i1<n1z;i1++) 
      d[i2][i1]=0;
  
  return;

}

