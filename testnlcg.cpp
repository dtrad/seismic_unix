/* Reading formatted file data with fscanf(). */
#include "su.h"
#include "dan.h"
#include <stdlib.h>
#include <stdio.h>
char *sdoc[] = {
  " TESTNLCG                                                                       ",
  " testnlcg program for testing different nonlinear solvers in a system of        ",	
  " equations or a cost function. 		         			   ",
  NULL};
float mpcgne(float **A, int n, int m, float *x,float *b,float *MI, float tol, float step,int itercg,int restart);

void Atimesx(float *b,float **A,float *x,int adj, int n,int m);

float nlcgtest_interface(float *b,float **A,float *x, int nd, int nx, float ftol, int itercg,float lambda);

int main(int argc, char **argv)
{
     FILE *fp, *fp2;
     int i,j,N=128,M=2*N;
     float **A, *b, *x, *MI;
     initargs(argc, argv);
     int itercg=N;
     float step=1;
     float tol=1e-7;
     int restart=0;
     float ftol;
     float lambda;
     int method;

     if (!getparfloat("ftol", &ftol))  ftol = 1e-2;
     if (!getparfloat("lambda", &lambda))  lambda = 1e-1;     
     if (!getparint("method", &method))  method = 1;     


     if ( (fp = fopen("inputnlcg.txt", "r")) == NULL)
     {
         fprintf(stderr, "Error opening file.\n");
         exit(1);
     }
     if ( (fp2 = fopen("outputnlcg.txt", "w")) == NULL)
     {
         fprintf(stderr, "Error opening file.\n");
         exit(1);
     }
     float fN,fM;
     fscanf(fp, "%f", &fN);
     fscanf(fp, "%f", &fM);
     N=(int) fN;
     M=(int) fM;
     fprintf(stderr,"N=%d ,M=%d\n",N,M);

     //requestdoc(1); 
     if ((A=alloc2float(M,N))==NULL)
       err("Cannot allocate A\n");
     if ((b=alloc1float(N))==NULL)
       err("Cannot allocate b\n");
     if ((x=alloc1float(M))==NULL)
       err("Cannot allocate x\n");
     if ((MI=alloc1float(M))==NULL)
       err("Cannot allocate MI\n");


     for (i=0;i<M;i++) MI[i]=1;
     //MI[5]=1e-5;

     for (i=0;i<N;i++)
       for (j=0;j<M;j++)
	 fscanf(fp, "%f", &A[i][j]);
     if (0)
       for (i=0;i<N;i++){
	 for (j=0;j<M;j++)
	   printf("%4.1f ",A[i][j]);
	 printf("\n");
       }

     for (i=0;i<N;i++) fscanf(fp,"%f",&b[i]);
     //for (i=0;i<N;i++) printf("%4.1f ",b[i]);
     printf("\n");

     if (0){
       Atimesx(b,A,x,1,N,M);
       for (i=0;i<M;i++) printf("x[%d]=%f\n",i,x[i]);
       Atimesx(b,A,x,0,N,M);
       for (i=0;i<N;i++) printf("b[%d]=%f\n",i,b[i]);     
     }
     printf("\n");
     TRACE;
     if (method==1) nlcgtest_interface(b,A,x,N,M,ftol,itercg,lambda);
     else if (method==2) mpcgne(A,N,M,x,b,MI,tol,step,itercg,restart);
     //for (i=0;i<M;i++) printf("MI[%d]=%f\n",i,MI[i]); 
     //for (i=0;i<M;i++) printf("x[%d]=%f\n",i,x[i]); 
     for (i=0;i<M;i++) fprintf(fp2,"%f\n",x[i]); 

     fclose(fp2);
     fclose(fp);
     free1float(b);
     free1float(x);
     free1float(MI);
     free2float(A);
     return 0;
}








