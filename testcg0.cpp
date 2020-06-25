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



 /* 
     This function solves the system of equations 
     (WmT FH WdT M^{-1} Wd F Wm) m =  FH WdT Wd M^{-1} d
     oper is an operator implemented with the function
     where M is a  preconditioner acting on data space.
     void oper (float *,float *,float *,float *,float *,int ,int ,int,int)
     Example 
     void oper (float x*,float t*,float h*,float p*,float d*,int adj, int nt,
     int nh ,int np);
     When adj=0 ==> forward operator   (adjoint == False)
     When adj=1 ==> adjoint operator   (adjoint == True )
     Wd allows to downweight bad data by using large values. 
     In general is a diagonal size(ndata) 
     M is the preconditioner. Diagonal size(nmodel) Large values
     penalize, small values focus the solution to a desired model value.
     M changes the null space, W does not. 
     Hence prior information about the model must be implemented by M. 
     
     Taken from Yousef Saad, pag 260: 
     Iterative methods for sparse linear systems
     W has bee added to the original algorihm 
     Daniel Trad - UBC March-2000
    
  */

float mpcgne(float **A, int n, int m, float *x,float *b,float *MI, float tol, float step,int itercg,int restart)
  
{
  float normb,nit;
  double beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j;
  float J; // Cost function and square residuals
  
  // Temp pointers
  // r  Conjugate Gradient (Residual)
  // g  Gradient
  // z  Preconditioned gradient
  // s  Search direction
  // w  Conjugate search direction
  // M  precondtioner on data space
  // Wd model and data weights respectively.
  
  float *u,*g,*r,*s,*w,*rho,*eta,rhold;
  
  int nx=m;
  int ny=n;  
  
  if ((u=alloc1float(ny))==NULL) 
    printf("cannot allocate memory for u\n");
  if ((g=alloc1float(nx))==NULL) 
    printf("cannot allocate memory for g\n");
  if ((s=alloc1float(ny))==NULL)
    printf("cannot allocate memory for s\n");
  if ((r=alloc1float(ny))==NULL)
    printf("cannot allocate memory for rcg\n");
  if ((w=alloc1float(ny))==NULL)
    printf("cannot allocate memory for Azcg\n");
  if ((rho=alloc1float(itercg+1))==NULL)
    printf("cannot allocate memory for rho\n");  
  if ((eta=alloc1float(itercg+1))==NULL)
    printf("cannot allocate memory for eta\n"); 

  normb=dot(ny,b,b);
  
  for (i=0;i<ny;i++) {
    u[i]=0.;
    s[i]=r[i]=b[i];
  }
  
  k=0;
  rhold=tol*2;
  while((rhold>tol)&&(k<itercg)){
    k++;
    Atimesx(s,A,g,1,n,m);
    for (i=0;i<nx;i++) g[i]*=MI[i];
    Atimesx(w,A,g,0,n,m);
    
    alphanum=dot(ny,r,r);
    alphaden=dot(ny,w,s);
    alpha=alphanum/alphaden;
    
    if (alphaden < 0.) printf("alphaden=%e\n",alphaden);
    if (alphaden < 1e+10 ){ 
      fprintf(stderr,"alphanum=%e,alphaden=%e,alpha=%e\n",
	      alphanum,alphaden,alpha);
      //break;
    }


    //alpha*=step;

    for(i=0;i<ny;i++) u[i]+=alpha*s[i];
    for(i=0;i<ny;i++) r[i]-=alpha*w[i];  

    rho[k]=dot(ny,r,r)/normb;

    beta=dot(ny,r,r)/alphanum;
    fprintf(stderr,"rho[%d]=%f,beta=%e\n",k,rho[k],beta);    
    for (i=0;i<ny;i++) s[i]=r[i]+beta*s[i];
    eta[k]=dot(ny,u,u);

  }

  // Because the number of iterations is acting as regularization 
  // parameter the cost function is || AM^{-1}m-b||^2

  J=dot(ny,r,r);  
  Atimesx(u,A,x,1,n,m);
  for (i=0;i<nx;i++) x[i]*=MI[i];

  fprintf(stderr,"iter=%d,J=%f\n",k,J);
  free1float(eta);
  free1float(rho);
  free1float(w);
  free1float(r);
  free1float(s);
  free1float(g);
  free1float(u);

  return(J);
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





