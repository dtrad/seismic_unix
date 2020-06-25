#include "su.h"
float dot(int n, float *a, float *b);

void cgstep(int iter,int n,float *x,float *g,float *s,int m,
float *rr, float *gg, float *ss)
{
  const float TINY=1e-7;
  register int i;
  float sds;
  float gdg;
  float determ;
  float gdr;
  float gds;
  float sdr;
  float alfa;
  float beta;
  if (iter == 0){
    for (i=0;i<n;i++) s[i]=0;
    for (i=0;i<m;i++) ss[i]=0;
    if (dot(m,gg,gg)==0) err("cgstep: grad vanishes identically\n");
    alfa = dot(m,gg,rr) / dot(m,gg,gg);
    beta = 0;
  }
  else{                 // search plane by solving 2 by 2
    gdg=dot(m,gg,gg);   // G.(R-G*alfa -S*beta) = 0
    sds=dot(m,ss,ss);   // S.(R-G*alfa -S*beta) = 0
    gds=dot(m,gg,ss);   
    determ=gdg * sds - gds * gds + TINY;
    gdr = dot(m,gg,rr);
    sdr = dot(m,ss,rr);
    alfa = (sds * gdr - gds * sdr ) / determ;
    beta = (-gds * gdr + gdg * sdr ) / determ;
   }
  for (i=0;i<n;i++) s[i]=alfa * g[i] + beta * s[i];    // model step
  for (i=0;i<m;i++) ss[i]=alfa * gg[i] + beta * ss[i]; // ss=conjugate
  for (i=0;i<n;i++) x[i]+=s[i];                        // update solution
  for (i=0;i<m;i++) rr[i]-=ss[i];                      // update residual
  return;
  
}
