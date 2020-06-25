#include "su.h"
#include "segy.h"
///////////////////////////////////////////////////////////////////////////////
void contran(int conj,int sum, int nx, float *xx, int nb, float *bb,float *yy);
void conjzero(int conj, int add, int nx, float *x, int ny, float *y);
void cgstep(int iter,int n,float *x,float *g,float *s,int m,
float *rr, float *gg, float *ss);
void scale(float factor, int ndata, float *data);
void copy(int n, float *xx, float *yy);
///////////////////////////////////////////////////////////////////////////////
void  missfip(int nt, int na,int np,float *data,float *known,int niter) 
  /*
Fill in missing data on 1 axis by minimizing power out of a  filter
The filter is estimated by using a training data set
Claerbout(1987)- Pag 185
   */
{
  
  float *aa; /* filter */
  float *pp;
  float *tt;
  float *x;
  float *g; /* dellta model */
  float *s;
  float *rr;  /*residual */
  float *gg;
  float *ss;

  if ((aa = alloc1float(na)) == NULL)
    err("cannot allocate memory for aa\n");
  if ((pp = alloc1float(np)) == NULL)
    err("cannot allocate memory for pp\n");
  if ((tt = alloc1float(nt)) == NULL)
    err("cannot allocate memory for tt\n");
  if ((x = alloc1float(np+na)) == NULL)
    err("cannot allocate memory for x\n");
  if ((g = alloc1float(np+na)) == NULL)
    err("cannot allocate memory for g\n");
  if ((s = alloc1float(np+na)) == NULL)
    err("cannot allocate memory for s\n");
  if ((rr = alloc1float(np+na-1+np+na-1)) == NULL)
    err("cannot allocate memory for rr\n");  
  if ((gg = alloc1float(np+na-1+np+nt-1)) == NULL)
    err("cannot allocate memory for gg\n");  
  if ((ss = alloc1float(np+na-1+np+na-1)) == NULL)
    err("cannot allocate memory for ss\n");

  for (int it=0;it<nt;it++) tt[it]=data[it]; // training data set
  for (int it=0;it<np;it++) 
    if (known[it]!=0) pp[it]=data[it];  // known data with missing values
    else pp[it]=0;
  //for (int it=0;it<np;it++) pp[it]=data[it]; //known data with missing values
  for (int it=0;it<na;it++) aa[it]=0; aa[0]=1;  // roughening filter

  // lengths of outputs of filtering
  int npa=np+na-1;
  int nta=nt+na-1;
  // lengths of unknowns and residuals
  int nx=np+na;
  int nr=npa+nta;
  // pointers
  int px=0;
  int qr=0;
  int ax=np;    
  int tr=npa;  

  copy(np,pp,&x[px]);
  copy(na,aa,&x[ax]);
            
  for (int iter=0; iter < niter ; iter++){
    contran(0,0,na,aa,np,pp,&rr[qr]);  /* r= a*p convolution */
    contran(0,0,na,aa,nt,tt,&rr[tr]);  // extend rr with train
    scale(-1.,nr,rr);  
    contran(1,0,na,aa,np,&g[px],&rr[qr]);  // g_pp(aa,rr) correlation 
    contran(1,0,np,pp,na,&g[ax],&rr[qr]);  // g_aa(pp,rr) correlation
    contran(1,1,nt,tt,na,&g[ax],&rr[tr]);  // g_aa(tt,rr) correlation
    // Do not change known elements of the data
    for (register int ip=0; ip < np; ip++)
      if (known[ip] !=0.)  g[ip+px] = 0;  /* missing data where copy(ip)==0 */
    g[ax]=0;  // Do not change first element of the filter
    contran(0,0,na,aa,np,&g[px],&gg[qr]);  // gg=aa*g_pp  convolution
    contran(0,1,np,pp,na,&g[ax],&gg[qr]);  // gg=pp*g_aa convolution
    contran(0,0,nt,tt,na,&g[ax],&gg[tr]);  // gg=tt*g_aa convolution
 
    cgstep(iter,nx,x,g,s,nr,rr,gg,ss);  // p=p+dp;  r=r-dr 
 
    copy(np,&x[px],pp);
    copy(na,&x[ax],aa);
    for (int it=0;it<na;it++) 
      fprintf(stderr,"aa[%d]=%f,iter=%d\n",it,aa[it],iter); 
  } 
  //for (int it=0;it<na;it++)fprintf(stderr,"aa[%d]=%f\n",it,aa[it]);
  for (int it=0;it<np;it++) data[it]=pp[it];  
  for (int it=0;it<na;it++) fprintf(stderr,"aa[%d]=%f\n",it,aa[it]);
  free1float(ss);
  free1float(gg);
  free1float(rr);
  free1float(s);
  free1float(g);
  free1float(x);
  free1float(tt);
  free1float(pp);
  free1float(aa);

  return;
}






