#include <math.h>
#include "su.h"

void radhypsp(float *t, float *h, float *q, float **vel, int nt, int nh, int nq,unsigned int **index)
{
  register unsigned int it;
  unsigned int ih,iq;
  float dint;
  float time,hxh,pxhxh;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  unsigned int nsparse=nt*nq*nh;
  float dt=t[1]-t[0];

  for (it=0;it<nsparse;it++) index[0][it]=index[1][it]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=0;it<nt;it++){
	pxhxh=hxh/vel[iq][it];
        time=sqrt(t[it]*t[it]+pxhxh);
	itime=(int) (time/dt+0.5);
	if (itime<nt){
	  index[0][ih*nq*nt+iqxnt+it]=ihxnt+itime;
          index[1][ih*nq*nt+iqxnt+it]=iqxnt+it;
	}
      }            
    }
  }
  return;
}

void radonhyp(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse)
{
  unsigned long j;
  unsigned int ny=nh*nt;
  unsigned int nx=nq*nt;
  //float nqnh=sqrt(nq*nh);
  
  if (!adj){
    for (j=0;j<ny;j++) d[j]=0;
    for (j=0;j<nsparse;j++) d[index[0][j]]+=m[index[1][j]];
    //for (j=0;j<ny;j++) d[j]/=nqnh;
  }
  else{
    for (j=0;j<nx;j++) m[j]=0;
    for (j=0;j<nsparse;j++) m[index[1][j]]+=d[index[0][j]];
    //for (j=0;j<nx;j++) m[j]/=nqnh;
  }
  d[0]=0;
  m[0]=0;
  return;
}

void radhypsp(float *t, float *h, float *q, float **vel, int nt, int nh, int nq,float **index)
{
  register unsigned int it;
  unsigned int ih,iq;
  unsigned int j;
  float *ttn,*dint,*tnt,hxh,pxhxh;
  unsigned int iqxnt,ihxnt;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  unsigned int ns=nt*nh*nq;

  float dt=t[1]-t[0];
  float dt2=dt*dt;

  if ((ttn=alloc1float(nt))==NULL)
    err("cannot allocate memory for ttn \n");
  if ((tnt=alloc1float(nt))==NULL)
    err("cannot allocate memory for tnt \n");
 
  for (it=0;it<ns;it++) index[0][it]=index[1][it]=0;

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=0;it<nt;it++){
	pxhxh=hxh/vel[iq][it];
        ttn[it]=sqrt(t[it]*t[it]/dt2+pxhxh/dt2);
	if (ttn[it]<nt) index[0][ih*nq*nt+iq*nt+it]=ttn[it];
      }
      yxtoxy(nt,1.0,0.0,ttn,nt,1.0,0.0,-nt,nt,tnt);
      for (it=0;it<nt;it++) if (tnt[it]<nt) index[1][ih*nq*nt+iq*nt+it]=tnt[it];
    }
  }
  free1float(ttn);
  free1float(tnt);
  return;
}

void radonhyp(float *m, float *d, float **index, int adj, int nt, int nh, int nq)
{
  unsigned long j;
  unsigned int it,ih,iq,iqxnt,ihxnt;
  
  unsigned int ny=nh*nt;
  unsigned int nx=nq*nt;
  unsigned int ns=nq*nt*nh;
  float *dint;
  float nqnh=sqrt(nq*nh);
  
  dint=ealloc1float(nt);
  
  if (!adj){
    for (j=0;j<ny;j++) d[j]=0;
    for (ih=0;ih<nh;ih++){
      ihxnt=ih*nt;
      for (iq=0;iq<nq;iq++){
	iqxnt=iq*nt;
	ints8r(nt,1.0,0,&m[iqxnt],0.0,0.0,nt,&index[1][ih*nq*nt+iqxnt],dint);
	for (it=0;it<nt;it++) d[ihxnt+it]+=dint[it];
      }
    }
    //for (j=0;j<ny;j++) d[j]/=nqnh;
  }
  else{
    for (j=0;j<nx;j++) m[j]=0;
    for (iq=0;iq<nq;iq++){
      iqxnt=iq*nt;
      for (ih=0;ih<nh;ih++){
	ihxnt=ih*nt;
	ints8r(nt,1.0,0,&d[ihxnt],0.0,0.0,nt,&index[0][ih*nq*nt+iqxnt],dint);
	for (it=0;it<nt;it++) m[iqxnt+it]+=dint[it];
      }
    }
    //for (j=0;j<nx;j++) m[j]/=nqnh;
  }
  free1float(dint);
  return;
}












