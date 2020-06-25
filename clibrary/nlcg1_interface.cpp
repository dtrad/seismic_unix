#include "nrutil.h"
#include "Complex.h"
#include "su.h"
#include "dan.h"

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


void nlcg1_interface(complex *d,complex **L,complex *m, float nh, float nq, float ftol, int itercg, )
{
  float *mr;
  

  float (*func)(float []);
  void (*dfunc)(float [], float []);

  nh_e=2*nh;
  nq_e=2*nq;

  func=costfunc;
  dfunc=gradient;

  // Transform the complex system of equation in a real one 
  dr_e=vector(1,2*nh);
  mr=vector(1,2*nq);
  LR_E=matrix(1,2*nq,1,2*nh);

  complex2real(d[freq],L,&dr[1],&LR[1][1],nh,nq);
  Atimesxnr(dp,LR_E,x,nh_e,nq_e,adj);  
  frprmn(mr,2*nq,ftol, itercg,&fret,func,dfunc);
  matmul0 (d2r,LR,m2r,0,2*nh,2*nq);
  real2complex(d2[freq],d2r,m2[freq],m2r,nh,nq);  

  free2float(LR);
  free1float(mr);
  free1float(dr);

}



float func1(float *x)
{
  extern int nh_e, nq_e;
  float Jt,Jd,Jm;
  int i,j, sum;
  extern float **LR_E, *dr_e, lambda;  

  dp=vector(1,nh_e);
  res=vector(1,nh_e);

  Atimesxnr(dp,LR_E,x,nh_e,nq_e,adj);  
  for (i=1;i<=nh_e;nh_e) res[i]=(dr_e[i]-dp[i]);
  for (Jd=0,i=1;i<=nh_e;nh_e) Jd+=res[i];;
  for (Jm=0,i=1;i<=nq_e;nq_e) Jm+=fabs(x[i]);

  free_vector(res,1,nh_e);
  free_vector(dp,1,nh_e);

  return(Jd+lambda_e*Jm);

}



void dfunc1(float *x, float *grad)
{
  extern int nh_e, nq_e;
  float Jt,Jd,Jm;
  int i,j, sum;
  extern float **LR_E, *dr_e, lambda;  

  dp=vector(1,nh_e);
  res=vector(1,nh_e);

  Atimesxnr(dp,LR_E,x,nh_e,nq_e,0);  
  for (i=1;i<=nh_e;nh_e) res[i]=-(dr_e[i]-dp[i]); /* gradient=A'(Ax-b)=A'(-res) */
  Atimesxnr(res,LR_E,grad,nh_e,nq_e,1); 

  for (i=1;i<=nq_e;nq_e) grad-=lambda_e*SGN(x[i]);

  free_vector(res,1,nh_e);
  free_vector(dp,1,nh_e);

  return;
}










































































