#include <math.h>
#include "su.h"
#include "radontd_sparse.h"

void contran(int conj,int sum, int nx, float *xx, int nb, float *bb,float *yy);
void contruc(int conj,int add, int lag, int nx, float *xx, int nb, float *bb,int ny,float *yy);
void convin(int conj,int sum, int nx, float *xx, int nb, float *bb,float *yy);

void build_index_rad(float *t, float *h, float *q, float **vel, int nt, int nh, int nq,unsigned int **index, int it0)
{
  register unsigned int it;
  unsigned int ih,iq;
  float dint;
  float time,hxh,pxhxh;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  unsigned int nsparse=(nt-it0)*nq*nh;
  float dt=t[1]-t[0];
  int it00;

  for (it=0;it<nsparse;it++) index[0][it]=index[1][it]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=it0;it<nt;it++){
	it00=it-it0;
	pxhxh=hxh/vel[iq][it];
        time=sqrt(t[it]*t[it]+pxhxh);
	itime=(int) (time/dt+0.5);
	if (itime<nt){
	  index[0][ih*nq*(nt-it0)+iq*(nt-it0)+it00]=ihxnt+itime;
          index[1][ih*nq*(nt-it0)+iq*(nt-it0)+it00]=iqxnt+it;
	}
      }            
    }
  }
  return;
}

void radonhyp(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse)
{
  unsigned int j;
  unsigned int ny=nh*nt;
  unsigned int nx=nq*nt;
  
  if (!adj){
     memset((void *) d,(int)'\0',ny*FSIZE);
    for (j=0;j<nsparse;j++) d[index[0][j]]+=m[index[1][j]];
   }
  else{
    memset((void *) m,(int)'\0',nx*FSIZE);
    for (j=0;j<nsparse;j++) m[index[1][j]]+=d[index[0][j]];
  }
  /* 
    A problem appears if some of the values of index are never computed
    because the zero index of d and m are mapped each other for index=0
    I make these two elements equal to zero just to prevent this problem, 
    It does not affect the data or model significantly.  
  */ 
			      
  d[0]=0;
  m[0]=0;

  return;
}

void radonhyp(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse, float *wavelet, int nw)
{
  int it,j,ih;
  int ny=nh*nt;
  int nx=nq*nt;
  int lag=(nw+1)/2;
  float *dtemp;
  //fprintf(stderr,"nw=%d,nt=%d,nh=%d,nq=%d\n",nw,nt,nh,nq);
  d[0]=0;
  m[0]=0;


  dtemp=ealloc1float(nt+nw);
  memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 

  //for (j=0;j<nw;j++) fprintf(stderr,"nw=%d, wavelet[%d]=%f\n",nw,j,wavelet[j]);
  
  if (!adj){
    memset((void *) d,(int)'\0',ny*FSIZE);
    for (j=0;j<nsparse;j++) d[index[0][j]]+=m[index[1][j]];
    d[0]=0;
    for (ih=0;ih<nh;ih++){
      memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 
      contruc_2(0,0,nw,wavelet,nt,&d[ih*nt],dtemp);
      memcpy((void *) &d[ih*nt],(const void *) dtemp,nt*sizeof(float));
    }
  }
  if (adj){
    for (ih=0;ih<nh;ih++){
      memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 
      contruc_2(1,0,nw,wavelet,nt,dtemp,&d[ih*nt]);
      memcpy((void *) &d[ih*nt],(const void *) dtemp,nt*sizeof(float));
    }
    memset((void *) m,(int)'\0',nx*FSIZE);
    for (j=0;j<nsparse;j++) m[index[1][j]]+=d[index[0][j]];
    m[0]=0;
  }  
  
  /* 
    A problem appears if some of the values of index are never computed
    because the zero index of d and m are mapped each other for index=0
    I make these two elements equal to zero just to prevent this problem, 
    It does not affect the data or model significantly.  
  */
  d[0]=0;
  m[0]=0;

  free1float(dtemp);
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
    memset((void *) d,(int)'\0',ny*FSIZE);
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
    memset((void *) m,(int)'\0',nx*FSIZE);   
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

void build_index_rad_inv(float *t, float *h, float *q, float **vel, int nt, int nh, int nq,unsigned int **index, int it0)
{
  register unsigned int it;
  unsigned int ih,iq;
  float dint;
  float time,hxh,pxhxh;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  unsigned int nsparse=(nt-it0)*nq*nh;
  float dt=t[1]-t[0];
  int it00;

  for (it=0;it<nsparse;it++) index[0][it]=index[1][it]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=it0;it<nt;it++){
	if (vel[iq][it]>0){
	  it00=it-it0;
	  pxhxh=hxh*vel[iq][it];
	  time=sqrt(t[it]*t[it]+pxhxh);
	  itime=(int) (time/dt+0.5);
	  if (itime<nt){
	    index[0][ih*nq*(nt-it0)+iq*(nt-it0)+it00]=ihxnt+itime;
	    index[1][ih*nq*(nt-it0)+iq*(nt-it0)+it00]=iqxnt+it;
	  }
	}
      }            
    }
  }
  return;
}



void build_index_rad_inv_FM(float *t, float *h, float *q, float **vel, int nt, int nh, int nq,unsigned int **index, int it0)
{
  register unsigned int it;
  unsigned int ih,iq;
  float dint;
  float time,hxh,pxhxh;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  unsigned int nsparse=(nt-it0)*nq*nh;
  float dt=t[1]-t[0];
  int it00;
  float zref=1000;

  for (it=0;it<nsparse;it++) index[0][it]=index[1][it]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=sqrt(h[ih]*h[ih]+zref*zref)-zref;
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=it0;it<nt;it++){
	if (vel[iq][it]>0){
	  it00=it-it0;
	  pxhxh=hxh*vel[iq][it];
	  time=t[it]+pxhxh;
	  itime=(int) (time/dt+0.5);
	  if (itime<nt){
	    index[0][ih*nq*(nt-it0)+iq*(nt-it0)+it00]=ihxnt+itime;
	    index[1][ih*nq*(nt-it0)+iq*(nt-it0)+it00]=iqxnt+it;
	  }
	}
      }            
    }
  }
  return;
}


void radonhyp_FM(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse)
{
  unsigned int j;
  unsigned int ny=nh*nt;
  unsigned int nx=nq*nt;
  //float nqnh=sqrt(nh);
  
  if (!adj){
     memset((void *) d,(int)'\0',ny*FSIZE);
    for (j=0;j<nsparse;j++) d[index[0][j]]+=m[index[1][j]];
   }
  else{
    memset((void *) m,(int)'\0',nx*FSIZE);
    for (j=0;j<nsparse;j++) m[index[1][j]]+=d[index[0][j]];
  }
  /* 
    A problem appears if some of the values of index are never computed
    because the zero index of d and m are mapped each other for index=0
    I make these two elements equal to zero just to prevent this problem, 
    It does not affect the data or model significantly.  
  */ 
			      
  d[0]=0;
  m[0]=0;
  return;
}















