#include "su.h"
#include "segy.h"
#include "Complex.h"
void xminusy(complex *z,complex *x,complex *y,int n)
{
     int i;

     for (i=0;i<n;i++)
       z[i]=x[i]-y[i];
     return;
}
void xminusy(float *z,float *x,float *y,int n)
{
     int i;

     for (i=0;i<n;i++)
       z[i]=x[i]-y[i];
     return;
}
