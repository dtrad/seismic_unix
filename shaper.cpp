#include "su.h"
#include "segy.h"


void contran(int conj,int sum, int nx, float *xx, int nb, float *bb,float *yy);
void conjzero(int conj, int add, int nx, float *x, int ny, float *y);
void cgstep(int iter,int n,float *x,float *g,float *s,int m,
float *rr, float *gg, float *ss);
 
void shaper(int nf,int nx,float *data, int niter)
  /*
Shaping filter:
Minimize SUM rr[i]^2 by finding ff and rr where
rr = yy - xx * ff
Claerbout(1987)- Pag 153
   */
{
  
  float *ff; /* filter */
  float *xx;
  float *yy;
  float *rr; 
  float *df; 
  float *sf;
  float *dr;
  float *sr;
  //if (ny!= nx+nf-1) err("ny is wrong\n");
  
  int ny=nx+nf-1;
  register int it;
  
  if ((ff = alloc1float(nf)) == NULL)
    err("cannot allocate memory for ff\n");
  if ((xx = alloc1float(nx)) == NULL)
    err("cannot allocate memory for xx\n");
  if ((rr = alloc1float(ny)) == NULL)
    err("cannot allocate memory for ff\n");
  if ((yy = alloc1float(ny)) == NULL)
    err("cannot allocate memory for yy\n");


  if ((df = alloc1float(nf)) == NULL)
    err("cannot allocate memory for df\n");
  if ((sf = alloc1float(nf)) == NULL)
    err("cannot allocate memory for sf\n");
  if ((dr = alloc1float(ny)) == NULL)
    err("cannot allocate memory for dr\n");  
  if ((sr = alloc1float(ny)) == NULL)
    err("cannot allocate memory for sr\n");  
  
  for (int it=0;it<nx;it++) xx[it]=data[it];

  for (int it=0;it<nf;it++) ff[it]=0;
  ff[0]=1;
  ff[1]=-0.5;
  ff[2]=0.3;
  
  contran(0,0,nx,xx,nf,ff,yy);

  //for (int it=nx;it<ny;it++) yy[it]=0;
  //for (int it=0;it<nx;it++) yy[it]=data[it];

  for (int it=0;it<nf;it++) ff[it]=0;
  for (int it=0;it<ny;it++) rr[it]=yy[it];    
     
  for (int iter=0; iter < niter ; iter++){
    contran(1,0,nx,xx,nf,df,rr);  /* dp(a,r) correlation */
    contran(0,0,nx,xx,nf,df,dr);  /* dr=a*dp convolution */
    cgstep(iter,nf,ff,df,sf,ny,rr,dr,sr);  // p=p+dp;  r=r-dr  
  }
  for (int it=0;it<nf;it++) fprintf(stderr,"ff[%d]=%f\n",it,ff[it]);
  for (int it=0;it<nx;it++) data[it]=rr[it];
  free1float(sr);
  free1float(dr);
  free1float(sf);
  free1float(df);
  free1float(rr);
  free1float(yy);
  free1float(xx);
  free1float(ff);

  return;
}




