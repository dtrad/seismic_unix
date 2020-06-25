#include "su.h"
#include "radonfrequency.h"
#define SIDE 1    // side 1 --> right mute ; side 1 left mute

/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/
void xplotgather(float **d, int nh, int nt, float dt, char *s, char *s2);

void adaptive_subtract(float **data1, float **data2, int nh, int nt);
void mute(float **model, int nq, int nt, float *q, float *t, float parmute, float t0, int side);

void radonsolver0_subtract(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, int itercg, int iter_end, float step, float eps2, float eps1, float quantil, int norm, float factor, float smute, float nmofactor, int rtmethod, float depth, float fmax, char *solver, float **M, int muteflag)
{
  int i, it, ih;
  float  qmax, qmaxt;
  float  dx_av, dx_max;
  float *dtemp;
  float dq;
  float qmin=q[0];
  float *Wd; 
  float eps=1e-7;
  int testadj=0;
  float **datapred;
  int side=SIDE;   // side 1 --> right mute ; side 1 left mute
  dtemp=ealloc1float(nt);
  Wd=ealloc1float(nh);
  datapred=ealloc2float(nt,nh);
   
  for (ih=0;ih<nh;ih++) Wd[ih]=1;

  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
  radon_param(fmax,h,nh,dx_av,qmin,&qmaxt,&qmax,&dq,nq,rtmethod,factor);
  fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
  fprintf(stderr,"freq max=%f,dq=%e, nh=%d, nt=%d\n", fmax,dq,nh,nt);
  for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 
  
  
  for (ih=0;ih<nh;ih++){
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
  }

  //  plotgather(data,nh,nt,dt,"suxwigb");
  if (0) plotgather_pipe(data,nh,nt,"After_NMO");
  memset( (void *) model[0], (int) '\0', nq * nt *FSIZE);  

  radonsolver(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,
	       norm,step,testadj,rtmethod,depth,solver);
  
  /* Any automatic filtering process can be done here on model[nq][nt] */
  if (1){
    xplotgather(model,nq,nt,dt,"model1","xbox=0 legend=1 perc=98&");
    xplotgather(M,nq,nt,dt,"Mask","xbox=600 legend=1&");
  }
  if (muteflag) AtimesB(model,M,nq,nt);
  if (1){
    xplotgather(model,nq,nt,dt,"model2","xbox=0 legend=1 perc=98 &");
  }  
  /* compute inverse RT or predicted data */
  hrrti(datapred,h,dt,model,q,fmax,nt,nh,nq,rtmethod,depth);
  /* subtract predicted data */
  /* To get predicted rather than difference set this */
  #define subtract 1
  if ((muteflag)&&(subtract)) adaptive_subtract(data,datapred,nh,nt);
  else{
    for (ih=0;ih<nh;ih++)
      for (it=0;it<nt;it++)
	data[ih][it]=datapred[ih][it];
  }  
  

   

  for (ih=0;ih<nh;ih++){
    for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
  }    

  free1float(Wd);
  free1float(dtemp);
  free2float(datapred);

  return;

}


void xplotgather(float **d, int nh, int nt, float dt, char *s, char *s2)
{
  char buf[120];
  save_gather(d,nh,nt,dt,s);
  sprintf(buf,"suximage < %s title=%s curve=curve1 npair=5 hbox=900 wbox=700 %s\n",s,s,s2);
  system(buf);
  return;
}






















