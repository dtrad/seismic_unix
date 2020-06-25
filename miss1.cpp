#include "su.h"
#include "segy.h"


void contran(int conj,int sum, int nx, float *xx, int nb, float *bb,float *yy);
void conjzero(int conj, int add, int nx, float *x, int ny, float *y);
void cgstep(int iter,int n,float *x,float *g,float *s,int m,
float *rr, float *gg, float *ss);
 
void miss1(int na,int np,float *data, int niter, float *known)
  /*
Fill in missing data on 1 axis by minimizing power ou of a given filter
Claerbout(1987)- Pag 179
   */
{
  
  float *a; /* filter */
  float *p;
  float *dp; /* dellta model */
  float *sp; 
  float *r;  /*residual */
  float *dr;
  float *sr;
  int nr;
  if ((a = alloc1float(na)) == NULL)
    err("cannot allocate memory for a\n");
  if ((p = alloc1float(np)) == NULL)
    err("cannot allocate memory for p\n");
  if ((dp = alloc1float(np)) == NULL)
    err("cannot allocate memory for dp\n");
  if ((sp = alloc1float(np)) == NULL)
    err("cannot allocate memory for sp\n");
  if ((r = alloc1float(np+na-1)) == NULL)
    err("cannot allocate memory for r\n");  
  if ((dr = alloc1float(np+na-1)) == NULL)
    err("cannot allocate memory for dr\n");  
  if ((sr = alloc1float(np+na-1)) == NULL)
    err("cannot allocate memory for sr\n");  
  //for (int it=0;it<np;it++) p[it]=data[it];
  for (int it=0;it<np;it++) 
    if (known[it]!=0) p[it]=data[it];  // known data with missing values
    else p[it]=0;
  a[0]=1;
  a[1]=-1;
  a[2]=0;
   
  nr= np + na - 1;
  contran(0,0,na,a,np,p,r);  /* r= a*p convolution */
  for (int ir=0;ir<nr;ir++) r[ir]*=-1;           
  for (int iter=0; iter < niter ; iter++){
    contran(1,0,na,a,np,dp,r);  /* dp(a,r) correlation */
    for (register int ip=0; ip < np; ip++)
      if (known[ip] !=0.)  dp[ip] = 0;  /* missing data where known(ip)==0 */
    contran(0,0,na,a,np,dp,dr);  /* dr=a*dp convolution */
    cgstep(iter,np,p,dp,sp,nr,r,dr,sr);  // p=p+dp;  r=r-dr  
  }
  for (int it=0;it<np;it++) data[it]=p[it];
  free1float(sr);
  free1float(dr);
  free1float(r);
  free1float(sp);
  free1float(p);
  free1float(a);

  return;
}
