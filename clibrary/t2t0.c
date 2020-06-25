#include <math.h>
int t2t0(float t,float h,float v, float dt)
{
  /* Given t returns the index for t0 */  
  float t0;
  if ((t0=t*t-4*h*h/(v*v)) >= 0)
    return((int) (sqrt(t0)/dt+0.5));
  else return(0);
}

