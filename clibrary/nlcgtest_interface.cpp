#include "nrutil.h"
#include "su.h"
#include "dan.h"
#include "radonnlcg.h"
float f1dim(float x);
float df1dim(float x);
void frprmn1(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []), int deriv);
void dfunc1(float *x, float *xc);
float func1(float *x);
void linmin(float p[], float xi[], int n, float *fret,float (*func)(float []));
void dlinmin(float p[], float xi[], int n, float *fret, float (*func)(float []),
	     void (*dfunc)(float [], float []));

// External variables to comunicate the cost function with the cg directly
float **A_E;
float *d_e;
float lambda_e;
int nd_e, nx_e;


float nlcgtest_interface(float *d,float **A,float *x, int nd, int nx, float ftol, int itercg,float lambda)
{
  float fret;
  int itercount, ix, id;
  float *xnr;
  lambda_e=lambda;
  float (*func)(float []);
  void (*dfunc)(float [], float []);
  func=func1;
  dfunc=dfunc1;

  nd_e=nd;
  nx_e=nx;
  intprint(nd_e);
  intprint(nx_e);

  // Transform the variables from zero first index to one first index and extern
  d_e=vector(1,nd_e);
  //x=vector(1,nx_e);
  xnr=x-1;
  A_E=matrix(1,nd_e,1,nx_e);

  for (id=1;id<=nd_e;id++) d_e[id]=d[id-1];
  for (ix=1;ix<=nx_e;ix++) xnr[ix]=0;
  for (id=1;id<=nd_e;id++)
    for (ix=1;ix<=nx_e;ix++)  
      A_E[id][ix]=A[id-1][ix-1];
  
  //Atimesxnr(d_e,A_E,xnr,nd_e,nx_e,1);  
  frprmn1(xnr,nx_e,ftol,&itercount,&fret,func,dfunc,0);
  fprintf(stderr,"itercount=%d,ftol=%f,fret=%f\n",itercount,ftol,fret);
  //Atimesxnr(dr_e,LR_E,mr,nh_e,nq_e,0);

  free_matrix(A_E,1,nd_e,1,nx_e);
  free_vector(d_e,1,nd_e);

  return(fret);

}



float func1(float *x)
{
  extern int nd_e, nx_e;
  float Jd,Jm;
  int i;
  extern float **A_E, *d_e;  
  float *dp;
  float *res;
 
  dp=vector(1,nd_e);
  res=vector(1,nd_e);
  
  Atimesxnr(dp,A_E,x,nd_e,nx_e,0);  
  for (i=1;i<=nd_e;i++) res[i]=(d_e[i]-dp[i]);
  for (Jd=0,i=1;i<=nd_e;i++) Jd+=res[i]*res[i];

  if (1) for (Jm=0,i=1;i<=nx_e;i++) Jm+=fabs(x[i]);
  else for (Jm=0,i=1;i<=nx_e;i++) Jm+=(x[i]*x[i]);

  free_vector(res,1,nd_e);
  free_vector(dp,1,nd_e);

  //fprintf(stderr,"Jd=%f,Jm=%f\n",Jd,Jm);
  return(Jd+lambda_e*Jm);

}



void dfunc1(float *x, float *grad)
{
  extern int nd_e, nx_e;
  int i;
  extern float **A_E, *d_e;  
  float *dp;
  float *res;

  //TRACE; 
  dp=vector(1,nd_e);
  res=vector(1,nd_e);

  Atimesxnr(dp,A_E,x,nd_e,nx_e,0);  
  for (i=1;i<=nd_e;i++) res[i]=(-1)*(d_e[i]-dp[i]); /* gradient=A'(Ax-b)=A'(-res) */
  Atimesxnr(res,A_E,grad,nd_e,nx_e,1); 
  for (i=1;i<=nx_e;i++) grad[i]*=2;
  if (1) for (i=1;i<=nx_e;i++) grad[i]+=(lambda_e * SGN(x[i]));
  else for (i=1;i<=nx_e;i++) grad[i]+=(lambda_e * x[i]);
  free_vector(res,1,nd_e);
  free_vector(dp,1,nd_e);

  return;
}

#define NRANSI
#define ITMAX 500
#define EPS 1.0e-10
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);

void frprmn1(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []), int deriv)
{
	void linmin(float p[], float xi[], int n, float *fret,
		float (*func)(float []));
	void dlinmin(float p[], float xi[], int n, float *fret, float (*func)(float []),
		     void (*dfunc)(float [], float []));

	int j,its;
	float gg,gam,fp,dgg;
	float *g,*h,*xi;

	g=vector(1,n);
	h=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		if (deriv==0) linmin(p,xi,n,fret,func);
		//else dlinmin(p,xi,n,fret,func,dfunc);
		fprintf(stderr,"fret=%f\n",*fret);
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return;
		}
		fp=(*func)(p);
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	nrerror("Too many iterations in frprmn");
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software . */









































