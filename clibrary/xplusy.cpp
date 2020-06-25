#include "su.h"
#include "segy.h"
#include "Complex.h"
void xplusy(complex *z,complex *x,complex *y,int n)
{
     int i;

     for (i=0;i<n;i++)
       z[i]=x[i]+y[i];
     return;
}

void Aplus_equalB(complex **A,complex **B,int n1, int n2)
{
     int i1, i2;
     for (i2=0;i2<n2;i2++)
       for (i1=0;i1<n1;i1++)
	 A[i2][i1]+=B[i2][i1];
       
     return;
}
void Aplus_equalB(float **A,float **B,int n1, int n2)
{
     int i1, i2;
     for (i2=0;i2<n2;i2++)
       for (i1=0;i1<n1;i1++)
	 A[i2][i1]+=B[i2][i1];
       
     return;
}
