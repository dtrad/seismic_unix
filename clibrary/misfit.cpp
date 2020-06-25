#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrary.h"
float misfit(complex *a,complex *b,float powd,float *dh, int n)
{
     int i;
     float Jdata;

     for (Jdata=0,i=0;i<n;i++)
         Jdata+=(dh[i]*real((a[i]-b[i])*conjg(a[i]-b[i])));         
     Jdata/=powd; 
     return(Jdata);
}
