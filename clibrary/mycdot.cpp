#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"
#include "/usr/local/cwp/src/Complex/include/Complex.h"

complex mycdot(complex *x,complex *y,int n)
{
  //     Compute the inner product
  //     dot=(x,y) for complex x,y
        
   int i;
   complex  cdot, val;
   val.r=0;
   val.i=0;
   for (i=0;i<n;i++)  
       val = val +conjg(x[i])*y[i];
   cdot=val;
   return(val);
}





