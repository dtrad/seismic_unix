#include <stdio.h>
void atimes1(double **A,double x[],double b[],int n,int itrnsp)
{
     int i,j;

     if (itrnsp){
       for (i=1;i<=n;i++){
         b[i]=0;
         for (j=1;j<=n;j++)
	   b[i]+=A[j][i]*x[j];
       }
     }
     else{
     for (i=1;i<=n;i++){
         b[i]=0;
         for (j=1;j<=n;j++)
	   b[i]+=A[i][j]*x[j];
       }
     }
     return;
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */
