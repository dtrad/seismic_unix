/* Reading formatted file data with fscanf(). */
#include "su.h"
#include <stdlib.h>
#include <stdio.h>
char *sdoc[] = {
  " 	   								",
  NULL};
float mpcgne(float **A, int n, int m, float *x,float *b,float *MI, float tol, float step,int itercg,int restart);
void Atimesx(float *b,float **A,float *x,int adj, int n,int m);
int main(int argc, char **argv)
{
     FILE *fp;
     int i,j,N=5,M=6;
     float **A, *b, *x, *MI;
     initargs(argc, argv);
     int itercg=N;
     float step=1;
     float tol=1e-7;
     int restart=0;
     
     //requestdoc(1); 
     if ((A=alloc2float(M,N))==NULL)
       err("Cannot allocate A\n");
     if ((b=alloc1float(N))==NULL)
       err("Cannot allocate b\n");
     if ((x=alloc1float(M))==NULL)
       err("Cannot allocate x\n");
     if ((MI=alloc1float(M))==NULL)
       err("Cannot allocate MI\n");

     if ( (fp = fopen("input.txt", "r")) == NULL)
     {
         fprintf(stderr, "Error opening file.\n");
         exit(1);
     }
     for (i=0;i<M;i++) MI[i]=1;
     MI[5]=1e-5;

     for (i=0;i<N;i++)
       for (j=0;j<M;j++)
	 fscanf(fp, "%f", &A[i][j]);

     for (i=0;i<N;i++){
       for (j=0;j<M;j++)
	 printf("%4.1f ",A[i][j]);
       printf("\n");
     }

     for (i=0;i<N;i++) fscanf(fp,"%f",&b[i]);
     for (i=0;i<N;i++) printf("%4.1f ",b[i]);
     printf("\n");

     if (0){
       Atimesx(b,A,x,1,N,M);
       for (i=0;i<M;i++) printf("x[%d]=%f\n",i,x[i]);
       Atimesx(b,A,x,0,N,M);
       for (i=0;i<N;i++) printf("b[%d]=%f\n",i,b[i]);     
     }
     printf("\n");

     mpcgne(A,N,M,x,b,MI,tol,step,itercg,restart);
     for (i=0;i<M;i++) printf("MI[%d]=%f\n",i,MI[i]); 
     for (i=0;i<M;i++) printf("x[%d]=%f\n",i,x[i]); 
   

     fclose(fp);
     free1float(b);
     free1float(x);
     free1float(MI);
     free2float(A);
     return 0;
}


void Atimesx(float *b,float **A,float *x,int adj, int n,int m)
{
     int i,j;
  
     if (!adj) for (i=0;i<n;i++) b[i]=0;
     else for (i=0;i<m;i++) x[i]=0;

     for (i=0;i<n;i++){
         for (j=0;j<m;j++)
	   if (!adj) b[i]+=A[i][j]*x[j];
	   else x[j]+=A[i][j]*b[i];
     }
     return;
}

float dot(int n, float *a, float *b)
/********************************************************************  
return the  dot product
*********************************************************************/
{
	int j;
	float sum=0.;
	for(j=0;j<n;j++) sum += (a[j]*b[j]);
	return(sum);
}








