#include "su.h"
#include "segy.h"

void contruc(int conj,int add, int lag, int nx, float *xx, int nb, float *bb,int ny,float *yy);
void scale(float factor, int ndata, float *data);
void conjzero(int conj, int add, int nx, float *x, int ny, float *y);
void cgstep(int iter,int n,float *x,float *g,float *s,int m,
float *rr, float *gg, float *ss);
 
void iner(int nf,int nr,float *data, int niter, int lag, int gap1, int gapn)
  /*
PE filter:
Minimize SUM rr[i]^2 by finding ff and rr where
rr = yy - xx * ff
Claerbout(1987)- Pag 163
   */
{
  
  float *ff; /* filter */
  float *ww;
  float *yy;
  float *rr; 
  float *df; 
  float *sf;
  float *dr;
  float *sr;
  float *wr;
  register int it;  

  if ((lag<gap1)||(lag>gapn)) err("Input fails gap1<=lag<=gapn\n");
    
  if ((ff = alloc1float(nf)) == NULL)
    err("cannot allocate memory for ff\n");
  if ((ww = alloc1float(nr)) == NULL)
    err("cannot allocate memory for xx\n");
  if ((rr = alloc1float(nr)) == NULL)
    err("cannot allocate memory for ff\n");
  if ((yy = alloc1float(nr)) == NULL)
    err("cannot allocate memory for yy\n");
  if ((df = alloc1float(nf)) == NULL)
    err("cannot allocate memory for df\n");
  if ((sf = alloc1float(nf)) == NULL)
    err("cannot allocate memory for sf\n");
  if ((dr = alloc1float(nr)) == NULL)
    err("cannot allocate memory for dr\n");  
  if ((sr = alloc1float(nr)) == NULL)
    err("cannot allocate memory for sr\n");  
  if ((wr = alloc1float(nr)) == NULL)
    err("cannot allocate memory for wr\n");  
  
  for (int it=0;it<nr;it++) yy[it]=data[it];
  
  for (int it=0;it<nf;it++) ff[it]=0;
  ff[lag]=1; // set output lag
  
  contruc(0,0,lag,nr,yy,nf,ff,nr,rr); // set residual
  
  for (int it=0;it<nr;it++) ww[it]=1;
  for (int it=0;it<nr;it++) wr[it]=rr[it]*ww[it];    

  scale(-1,nr,wr);
   
  for (int iter=0; iter < niter ; iter++){
    for(it=0;it<nr;it++) dr[it]=wr[it]*ww[it];
    contruc(1,0,lag,nr,yy,nf,df,nr,dr);  /* df(yy,gr) correlation */
    for(it=gap1;it<gapn;it++) df[it]=0;  // constrained lags
    contruc(0,0,lag,nr,yy,nf,df,nr,dr);  /* dr=yy*df convolution */
    for(it=0;it<nr;it++) dr[it]=dr[it]*ww[it];
    cgstep(iter,nf,ff,df,sf,nr,wr,dr,sr);  // f=f+df;  r=r-dr  
  }
  contruc(0,0,lag,nr,yy,nf,ff,nr,rr);
  for (int it=0;it<nf;it++) fprintf(stderr,"ff[%d]=%f\n",it,ff[it]);
  for (int it=0;it<nr;it++) data[it]=rr[it];

  free1float(wr);
  free1float(sr);
  free1float(dr);
  free1float(sf);
  free1float(df);
  free1float(rr);
  free1float(yy);
  free1float(ww);
  free1float(ff);

  return;
}




