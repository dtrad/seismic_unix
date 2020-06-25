#include <math.h>
#include "su.h"

// There are 3 implementations of each function build_index and eomig
// 1 - For unsigned short
// 2 - For unsigend int
// 3 - For float using sinc interpolation slower

void build_index(float *t, float *h, float **vgrid, int nt, int nh, 
		  int nq,unsigned short **index, int nonzero, float t0)
{
  unsigned int it;
  unsigned int ih,iq;
  unsigned int j;
  float time,ftime,hxh,moveout,slowness2;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  float a1;
  float a2;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];
  int it0=(int) (t0/dt+0.5);

  for (it=0;it<nonzero;it++) index[0][it]=index[1][it]=0;

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=it0;it<nt;it++){
	moveout=hxh/(vgrid[iq][it]);
	time=sqrt(t[it]*t[it]+4*moveout);
	itime=(int) (time/dt+0.5);
	if (itime<nt){
	  index[0][ih*nq*nt+iq*nt+it]=ihxnt+itime;
	  index[1][ih*nq*nt+iqxnt+it]=iqxnt+it;
	}
      }
    }
  }
  return;
}


void build_index_inv(float *t, float *h, float **vgrid, int nt, int nh, 
		  int nq,unsigned short **index, int nonzero, float t0)
{
  unsigned int it;
  unsigned int ih,iq;
  unsigned int j;
  float time,ftime,hxh,moveout,slowness2;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  float a1;
  float a2;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];
  int it0=(int) (t0/dt+0.5);

  for (it=0;it<nonzero;it++) index[0][it]=index[1][it]=0;

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=it0;it<nt;it++){
	moveout=hxh*(vgrid[iq][it]);
	time=sqrt(t[it]*t[it]+4*moveout);
	itime=(int) (time/dt+0.5);
	if (itime<nt){
	  index[0][ih*nq*nt+iq*nt+it]=ihxnt+itime;
	  index[1][ih*nq*nt+iqxnt+it]=iqxnt+it;
	}
	else{
	  index[0][ih*nq*nt+iq*nt+it]=0;
	  index[1][ih*nq*nt+iq*nt+it]=0;
	}
      }
    }
  }
  return;
}

void build_index(float *t, float *h, float **vgrid, int nt, int nh, 
		  int nq,unsigned int **index, int nonzero, float t0)
{
  unsigned int it;
  unsigned int ih,iq;
  unsigned int j;
  float time,ftime,hxh,moveout,slowness2;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  float a1;
  float a2;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];
  int it0=(int) (t0/dt+0.5);
  int it00;
  for (it=0;it<nonzero;it++) index[0][it]=index[1][it]=0;

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=it0;it<nt;it++){
	if (vgrid[iq][it]>1e6){
	  //it00=it-it0;
	  moveout=hxh/(vgrid[iq][it]);
	  time=sqrt(t[it]*t[it]+4*moveout);
	  itime=(int) (time/dt+0.5);
	  if ((itime<nt)&&(it>0)&&(it<nt)&&(it>it0)){
	    index[0][ih*nq*nt+iqxnt+it]=ihxnt+itime;
	    index[1][ih*nq*nt+iqxnt+it]=iqxnt+it;
	  }
	  else{
	    index[0][ih*nq*nt+iq*nt+it]=0;
	    index[1][ih*nq*nt+iq*nt+it]=0;
	  }
	}
      }
    }
  }
  return;
}

void build_index_inv(float *t, float *h, float **vgrid, int nt, int nh, 
		  int nq,unsigned int **index, int nonzero, float t0)
{
  unsigned int it;
  unsigned int ih,iq;
  unsigned int j;
  float time,ftime,hxh,moveout,slowness2;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  float a1;
  float a2;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];
  int it0=(int) (t0/dt+0.5);
  int it00;
  for (it=0;it<nonzero;it++) index[0][it]=index[1][it]=0;

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=it0;it<nt;it++){
	//if (vgrid[iq][it]>0){
	  //it00=it-it0;
	moveout=hxh*(vgrid[iq][it]);
	time=sqrt(t[it]*t[it]+4*moveout);
	itime=(int) (time/dt+0.5);
	if ((itime<nt)&&(itime>0)&&(it>0)&&(it<nt)){
	  index[0][ih*nq*nt+iqxnt+it]=ihxnt+itime;
	  index[1][ih*nq*nt+iqxnt+it]=iqxnt+it;
	}
	else{
	  index[0][ih*nq*nt+iq*nt+it]=0;
	  index[1][ih*nq*nt+iq*nt+it]=0;
	}
	//}
      }
    }
  }
  return;
}


void eomig(float *m, float *d, unsigned short **index, int adj, int nt, int nh, int nq, int nsparse)
{
  register unsigned int j;
  unsigned int ny=nh*nt;
  unsigned int nx=nq*nt;
  //float nxny=sqrt(nx*ny);
  
  if (!adj){
    for (j=0;j<ny;j++) d[j]=0;
    for (j=0;j<nsparse;j++) d[index[0][j]]+=m[index[1][j]];
    //for (j=0;j<ny;j++) d[j]/=nxny;
  }
  else{
    for (j=0;j<nx;j++) m[j]=0;
    for (j=0;j<nsparse;j++) m[index[1][j]]+=d[index[0][j]];
    //for (j=0;j<nx;j++) m[j]/=nxny;
  }
  d[0]=0;
  m[0]=0;
  return;
}

void eomig(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse)
{
  register unsigned int j;
  unsigned int ny=nh*nt;
  unsigned int nx=nq*nt;
  //float nxny=sqrt(nx*ny);
  
  if (!adj){
    for (j=0;j<ny;j++) d[j]=0;
    for (j=0;j<nsparse;j++) d[index[0][j]]+=m[index[1][j]];
    //for (j=0;j<ny;j++) d[j]/=nxny;
  }
  else{
    for (j=0;j<nx;j++) m[j]=0;
    for (j=0;j<nsparse;j++) m[index[1][j]]+=d[index[0][j]];
    //for (j=0;j<nx;j++) m[j]/=nxny;
  }
  d[0]=0;
  m[0]=0;
  return;
}



void build_index(float *t, float *h, float **vel, int nt, int nh, int nq,float **index, int ns, float t0 )
{
  register unsigned int it;
  unsigned int ih,iq;
  unsigned int j;
  float *ttn,*dint,*tnt,hxh,pxhxh;
  unsigned int iqxnt,ihxnt;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];
  unsigned int it0=(int) (t0/dt+0.5);
  float dt2=dt*dt;
  //fprintf(stderr,"Build nt=%d,nh=%d,nq=%d,ns=%d,t0=%f,dt=%f\n",nt,nh,nq,ns,t0,dt);
  if ((ttn=alloc1float(nt))==NULL)
    err("cannot allocate memory for ttn \n");
  if ((tnt=alloc1float(nt))==NULL)
    err("cannot allocate memory for tnt \n");
 
  for (it=0;it<ns;it++) index[0][it]=index[1][it]=0;

  for (ih=0;ih<nh;ih++){
    hxh=4.0*h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=it0;it<nt;it++){
	pxhxh=hxh/vel[iq][it];
        ttn[it]=sqrt(t[it]*t[it]/dt2+pxhxh/dt2);
        if (ttn[it]<nt) 
	index[0][ih*nq*nt+iq*nt+it]=ttn[it];
      }
      yxtoxy(nt,1.0,0.0,ttn,nt,1.0,0.0,0,0,tnt);
      for (it=0;it<nt;it++) if ((tnt[it]<nt)&&(tnt[it]>=0)) 
	index[1][ih*nq*nt+iq*nt+it]=tnt[it];
    }
  }
  free1float(ttn);
  free1float(tnt);
  return;
}

void build_index_inv(float *t, float *h, float **vel, int nt, int nh, int nq,float **index, int ns, float t0 )
{
  register unsigned int it;
  unsigned int ih,iq;
  unsigned int j;
  float *ttn,*dint,*tnt,hxh,pxhxh;
  unsigned int iqxnt,ihxnt;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];
  unsigned int it0=(int) (t0/dt+0.5);
  float dt2=dt*dt;
  //fprintf(stderr,"Build nt=%d,nh=%d,nq=%d,ns=%d,t0=%f,dt=%f\n",nt,nh,nq,ns,t0,dt);
  if ((ttn=alloc1float(nt))==NULL)
    err("cannot allocate memory for ttn \n");
  if ((tnt=alloc1float(nt))==NULL)
    err("cannot allocate memory for tnt \n");
 
  for (it=0;it<ns;it++) index[0][it]=index[1][it]=0;

  for (ih=0;ih<nh;ih++){
    hxh=4.0*h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=it0;it<nt;it++){
	pxhxh=hxh*vel[iq][it];
        ttn[it]=sqrt(t[it]*t[it]/dt2+pxhxh/dt2);
        if (ttn[it]<nt) 
	index[0][ih*nq*nt+iq*nt+it]=ttn[it];
      }
      yxtoxy(nt,1.0,0.0,ttn,nt,1.0,0.0,0,0,tnt);
      for (it=0;it<nt;it++) if ((tnt[it]<nt)&&(tnt[it]>=0)) 
	index[1][ih*nq*nt+iq*nt+it]=tnt[it];
    }
  }
  free1float(ttn);
  free1float(tnt);
  return;
}


void eomig(float *m, float *d, float **index, int adj, int nt, int nh, int nq, int nsparse)
{
  unsigned long j;
  unsigned int it,ih,iq,iqxnt,ihxnt;
  unsigned int ny=nh*nt;
  unsigned int nx=nq*nt;
  float *dint;
  fprintf(stderr,"Eomig with sinc\n");
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

















