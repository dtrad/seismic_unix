#include <math.h>
void getcritangle(int n,float *t0, float *va, float *vb, float *critangle)
{
  register int i;		/* counter				*/
  
  for (i = 1; i < n; i++) {
    critangle[i]=va[i]/vb[i];
  }


}
