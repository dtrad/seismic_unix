#include "su.h"
#include "segy.h"
#include "Complex.h"

complex cdot(complex *x,complex *y,int n)
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





