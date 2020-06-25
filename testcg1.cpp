/* Reading formatted file data with fscanf(). */
#include <stdlib.h>
#include <stdio.h>

/* Prototypes */

void **alloc2 (size_t n1, size_t n2, size_t size);
void free2 (void **p);
void *alloc1 (size_t n1, size_t size);
void free1 (void *p);
float *alloc1float(size_t n1);
float **alloc2float(size_t n1, size_t n2);
void free1float(float *p);
void free2float(float **p);
float mpcgne(float **A, int n, int m, float *x,float *b,float *MI, float tol, float step,int itercg,int restart);
void Atimesx(float *b,float **A,float *x,int adj, int n,int m);




int main(int argc, char **argv)
{
     FILE *fp;
     int i,j,N=5,M=6;
     float **A, *b, *x, *MI;

     int itercg=N;
     float step=1;
     float tol=1e-7;
     int restart=0;
     

     if ((A=alloc2float(M,N))==NULL)
       printf("Cannot allocate A\n");
     if ((b=alloc1float(N))==NULL)
       printf("Cannot allocate b\n");
     if ((x=alloc1float(M))==NULL)
       printf("Cannot allocate x\n");
     if ((MI=alloc1float(M))==NULL)
       printf("Cannot allocate MI\n");

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



/* allocate a 2-d array */
void **alloc2 (size_t n1, size_t n2, size_t size)
{
	size_t i2;
	void **p;

	if ((p=(void**)malloc(n2*sizeof(void*)))==NULL) 
		return NULL;
	if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
		free(p);
		return NULL;
	}
	for (i2=0; i2<n2; i2++)
		p[i2] = (char*)p[0]+size*n1*i2;
	return p;
}

/* free a 2-d array */
void free2 (void **p)
{
	free(p[0]);
	free(p);
}



void *alloc1 (size_t n1, size_t size)
{
	void *p;

	if ((p=malloc(n1*size))==NULL)
		return NULL;
	return p;
}
/* free a 1-d array */
void free1 (void *p)
{
	free(p);
}


float *alloc1float(size_t n1)
{
	return (float*)alloc1(n1,sizeof(float));
}

float **alloc2float(size_t n1, size_t n2)
{
	return (float**)alloc2(n1,n2,sizeof(float));
}

void free1float(float *p)
{
	free1(p);
}

/* free a 2-d array of floats */
void free2float(float **p)
{
	free2((void**)p);
}





