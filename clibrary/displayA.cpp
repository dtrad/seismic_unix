#include "su.h"
#include "segy.h"
#include "Complex.h"
void displayA(complex **A,int nr,int nc)
{
     int i,j;
     for (i=0;i<nr;i++){
         for (j=0;j<nc;j++){
           fprintf(stderr,"A.r[%d][%d]=%e,\t A.i[%d][%d]=%e\n",i,j,A[i][j].r
              ,i,j,A[i][j].i);
           fprintf(stderr,"A.r[%d][%d]=%e,\t A.i[%d][%d]=%e\n",j,i,A[j][i].r
              ,j,i,A[j][i].i);
          }
     } 
     return;
}
void displayA(complex *A,int n)
{
     int i;
     for (i=0;i<n;i++){
           fprintf(stderr,"A.r[%d]=%e,\t A.i[%d]=%e\n",i,A[i].r,i,A[i].i);
     } 
     return;
}
void displayA(float *A,int n)
{
     int i;
     for (i=0;i<n;i++){
           fprintf(stderr,"A[%d]=%e,\n",i,A[i]);
     } 
     return;
}





