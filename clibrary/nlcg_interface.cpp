#include "nrutil.h"
#include "Complex.h"
#include "su.h"
#include "dan.h"
#include "radonnlcg.h"
float f1dim(float x);
float df1dim(float x);
void frprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []));
void dfrprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []));

void dfunc1(float *x, float *xc);
float func1(float *x);
void linmin(float p[], float xi[], int n, float *fret,float (*func)(float []));
void dlinmin(float p[], float xi[], int n, float *fret, float (*func)(float []),
	     void (*dfunc)(float [], float []));

// External variables to comunicate with the cost function 
float **LR_E;
float *dr_e;
float lambda_e;
int nh_e, nq_e;


float nlcg_interface(complex *d,complex **L,complex *m, int nh, int nq, float ftol, int itercg,float eps2)
{
  float *mr;
  float fret;
  int itercount;
  lambda_e=eps2;
  float (*func)(float []);
  void (*dfunc)(float [], float []);
  func=func1;
  dfunc=dfunc1;

  nh_e=2*nh;
  nq_e=2*nq;

  // Transform the complex system of equation in a real one 
  dr_e=vector(1,nh_e);
  mr=vector(1,nq_e);
  LR_E=matrix(1,nh_e,1,nq_e);
  complex2realnr(d,L,dr_e,LR_E,nh,nq);
  Atimesxnr(dr_e,LR_E,mr,nh_e,nq_e,1);  
  frprmn(mr,nq_e,ftol,&itercount,&fret,func,dfunc);
  fprintf(stderr,"itercount=%d,ftol=%f,fret=%f\n",itercount,ftol,fret);
  Atimesxnr(dr_e,LR_E,mr,nh_e,nq_e,0);
  real2complexnr(d,dr_e,m,mr,nh,nq);  

  free_matrix(LR_E,1,nh_e,1,nq_e);
  free_vector(mr,1,nq_e);
  free_vector(dr_e,1,nh_e);

  return(fret);

}



float func1(float *x)
{
  extern int nh_e, nq_e;
  float Jd,Jm;
  int i;
  extern float **LR_E, *dr_e;  
  float *dp;
  float *res;
 
  dp=vector(1,nh_e);
  res=vector(1,nh_e);
  
  Atimesxnr(dp,LR_E,x,nh_e,nq_e,0);  
  for (i=1;i<=nh_e;i++) res[i]=(dr_e[i]-dp[i]);
  for (Jd=0,i=1;i<=nh_e;i++) Jd+=res[i]*res[i];
  for (Jm=0,i=1;i<=nq_e;i++) Jm+=fabs(x[i]);

  free_vector(res,1,nh_e);
  free_vector(dp,1,nh_e);

  //fprintf(stderr,"Jd=%f,Jm=%f\n",Jd,Jm);
  return(Jd+lambda_e*Jm);

}



void dfunc1(float *x, float *grad)
{
  extern int nh_e, nq_e;
  int i;
  extern float **LR_E, *dr_e;  
  float *dp;
  float *res;

  //TRACE; 
  dp=vector(1,nh_e);
  res=vector(1,nh_e);

  Atimesxnr(dp,LR_E,x,nh_e,nq_e,0);  
  for (i=1;i<=nh_e;i++) res[i]=(-1)*(dr_e[i]-dp[i]); /* gradient=2 A'(Ax-b) = A'(-res) */
  Atimesxnr(res,LR_E,grad,nh_e,nq_e,1); 

  for (i=1;i<=nq_e;i++) grad[i]=2*grad[i]+(lambda_e * SGN(x[i]));

  free_vector(res,1,nh_e);
  free_vector(dp,1,nh_e);

  return;
}

#define NRANSI
#include "nrutil.h"
#define ITMAX 500
#define EPS 1.0e-10
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);

void frprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []))
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
		linmin(p,xi,n,fret,func);
		//dlinmin(p,xi,n,fret,func,dfunc);
		//fprintf(stderr,"fret=%f\n",*fret);
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









































