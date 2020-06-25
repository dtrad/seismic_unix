/* Reading formatted file data with fscanf(). */
#include "su.h"
#include "dan.h"
#include <stdlib.h>
#include <stdio.h>

char *sdoc[] = {
  " TESTNLCG                                                                       ",
  " testnlcg program for testing different linear solvers                          ",	
  " for complex matrices.        		         			   ",
  NULL};


void inverse_matrix_multiply (int nrows1, complex **matrix1, int ncols2,
			      int nrows2, complex **matrix2, complex **out_matrix);
void Atimesx(complex *b, complex **A, complex *x, int n,int m, int adj);

/*** dftls returns the transpose of the LS matrix for use with the Atimesx(..adjoint) */
void dftls(complex **F, complex **FLST, int nh, int nk, float eps);
/*** LSmatrix returns the LS matrix for use with the normal Atimesx */
void LSmatrix(complex **F, complex **FLS, int nh, int nk, float eps);

int main(int argc, char **argv)
{
     initargs(argc, argv);
     //requestdoc(1);
     FILE *fp, *fp2;
     int i,j,N=6,M=5;
     complex **A, *b, *x, *MI;
     complex **ALS, **ALS2;


     int itercg=N;
     float step=1;
     float tol=1e-7;
     int restart=0;
     float ftol;
     float lambda;
     int method;
     float fN,fM;

     if (!getparfloat("ftol", &ftol))  ftol = 1e-2;
     if (!getparfloat("lambda", &lambda))  lambda = 1e-5;     
     if (!getparint("method", &method))  method = 1;     


     fp = efopen("complexinput.txt", "r");
     fp2 = efopen("outputnlcg.txt", "w");


     fscanf(fp, "%f", &fN);
     fscanf(fp, "%f", &fM);
     N=(int) fN;
     M=(int) fM;
     fprintf(stderr,"N=%d ,M=%d\n",N,M);
     
     A=ealloc2complex(M,N);
     b=ealloc1complex(N);
     x=ealloc1complex(M);
     ALS=ealloc2complex(M,N);
     ALS2=ealloc2complex(N,M);
 
     for (j=0;j<M;j++)
       for (i=0;i<N;i++)
	 fscanf(fp, "%f%f", &A[i][j].r, &A[i][j].i);
     

     if (0)
       for (i=0;i<N;i++){
	 for (j=0;j<M;j++){
	   printf("%4.1f ",A[i][j].r);
	   printf("%4.1f ",A[i][j].i);
	 }
	 printf("\n");
       }

     for (i=0;i<N;i++){
       fscanf(fp,"%f",&b[i].r);
       fscanf(fp,"%f",&b[i].i);
     }

     printf("b=[Real Imag]\n");       
     for (i=0;i<N;i++){
       printf("%4.2f ",b[i].r);
       printf("%4.2f ",b[i].i);
     }
     printf("\n");

     if (0){
       Atimesx(b,A,x,N,M,1);
       for (i=0;i<M;i++) printf("x[%d]=%f\n",i,x[i].r);
       Atimesx(b,A,x,N,M,0);
       for (i=0;i<N;i++) printf("b[%d]=%f\n",i,b[i].r);     
     }

     printf("\n");
     if (method==0){
       Atimesx(b,A,x,N,M,1);
     }
     else if (method==1){
       dftls(A,ALS,N,M,lambda);
       Atimesx(b,ALS,x,N,M,1);
     }
     else if (method==2){
       LSmatrix(A,ALS2,N,M,lambda);
       Atimesx(x,ALS2,b,M,N);
     }    

     printf("x=[Real Imag]\n");
     for (i=0;i<M;i++){
       printf("%4.4f ",x[i].r);
       printf("%4.4f ",x[i].i);
     }
     printf("\n");

     if (0){
       for (i=0;i<N;i++){
	 fprintf(fp2,"%f %f",b[i].r,b[i].i); 
	 fprintf(fp2,"\n"); 
       }
     }


     if (1){
       for (i=0;i<M;i++){
	 fprintf(fp2,"%f %f",x[i].r,x[i].i); 
	 fprintf(fp2,"\n"); 
       }
     }


     fclose(fp2);
     fclose(fp);
     free1complex(b);
     free1complex(x);
     free2complex(A);
     free2complex(ALS);
     free2complex(ALS2);
     return 0;
}








