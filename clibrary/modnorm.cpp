#include "su.h"
#include "segy.h"
#include "Complex.h"
#include "clibrary.h"
float modnorm(complex *u,float eps2, float dq, int n)
{
     int i;
     float Jmod;
     extern int norm;
     if (norm==10)
     for (Jmod=0,i=0;i<n;i++)
         Jmod+=(dq*log(1+real(u[i]*conjg(u[i]))/eps2));
     if (norm==1)
     for (Jmod=0,i=0;i<n;i++)
         Jmod+=(dq*(abs(u[i])));

     return(Jmod);
}


float modnorm(complex *u,float eps2, float dq, int n, int norm)
{
     int i;
     float Jmod;

     if (norm==10)
     for (Jmod=0,i=0;i<n;i++)
         Jmod+=(dq*log(1+real(u[i]*conjg(u[i]))/eps2));
     if (norm==1)
     for (Jmod=0,i=0;i<n;i++)
         Jmod+=(dq*(abs(u[i])));

     return(Jmod);
}









