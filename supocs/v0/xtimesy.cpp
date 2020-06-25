#include "su.h"
#include "Complex.h"

void xtimesy(complex *z,complex *x,complex *y,int n)
{
     int i;

     for (i=0;i<n;i++)
       z[i]=x[i]*y[i];
     return;
}

void xtimesy(complex *z,complex *x,float *y,int n)
{
     int i;

     for (i=0;i<n;i++)
       z[i]=x[i]*y[i];
     return;
}

void xtimesy(float *z,float *x,float *y,int n)
{
     int i;

     for (i=0;i<n;i++)
       z[i]=x[i]*y[i];
     return;
}
void xtimesy(complex *z,float *x,complex *y,int n)
{
     int i;

     for (i=0;i<n;i++)
       z[i]=x[i]*y[i];
     return;
}
void xtimesy(complex *y,complex *x,float a,int n)
{
     int i;

     for (i=0;i<n;i++)
       y[i]=x[i]*a;
     return;
}
void zero_array(complex **d,int n2, int n1)
{
     int i1,i2;

     for (i2=0;i2<n2;i2++)
       for (i1=0;i1<n1;i1++)
	 d[i2][i1].r=d[i2][i1].i=0;

     return;
}
void zero_array(complex ***d,int n3, int n2, int n1)
{
     int i1,i2,i3;

     for (i3=0;i3<n3;i3++)
       for (i2=0;i2<n2;i2++)
	 for (i1=0;i1<n1;i1++)
	   d[i3][i2][i1].r=d[i3][i2][i1].i=0;
     
     return;
}



void zero_array(float **d,int n2, int n1)
{
     int i2;

     for (i2=0;i2<n2;i2++)
       memset((void*) d[i2],(int) '\0',n1*FSIZE);
 
     return;
}

void zero_array(float ***d, int n3, int n2, int n1)
{
     int i2,i3;
     
     for (i3=0;i3<n3;i3++)
       for (i2=0;i2<n2;i2++)
	 memset((void*) d[i3][i2],(int) '\0',n1*FSIZE);
 
     return;
}

void zero_vector(complex *d,int n)
{
     int i;

     for (i=0;i<n;i++) d[i].r=d[i].i=0.0;

     return;
}

void zero_vector(float *d,int n)
{
     int i;

     for (i=0;i<n;i++) d[i]=0;

     return;
}



