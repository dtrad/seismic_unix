#include <math.h>
#include <stdio.h>
#include <math.h>
#include "su.h"


void radonhyp(float *m, float *t, float *h, float *q, float *d, int adj, int nt, int nh, int nq) 
{
  register int it;
  int ih,iq;
  unsigned int j;
  float *ttn,*dint,*tnt,hxh,pxhxh;
  int iqxnt,ihxnt,itint;
  int nx=nt*nq;
  int ny=nt*nh;
  float dt=t[1]-t[0];
  float dt2=dt*dt;
  
  if ((ttn=alloc1float(nt))==NULL)
    err("cannot allocate memory for ttn \n");
  if ((tnt=alloc1float(nt))==NULL)
    err("cannot allocate memory for tnt \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  fprintf(stderr,"NEW APPROACH+++++++++++++++++++++++++++++++++++\n");

  if (adj) memset((void *) m,(int)'\0',nx*FSIZE);  
  //for (it=0;it<nx;it++) m[it]=0;
  else  memset((void *) d,(int)'\0',ny*FSIZE);  
  //for (it=0;it<ny;it++) d[it]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      pxhxh=hxh*q[iq];
      iqxnt=iq*nt;
      if (adj){
	for (it=0;it<nt;it++) ttn[it]=sqrt(t[it]*t[it]/dt2+pxhxh/dt2);
	ints8r(nt,1.0,0.,&d[ihxnt],0.0,0.0,nt,ttn,dint);	
      }
      else{
	for (it=0;it<nt;it++) ttn[it]=sqrt(t[it]*t[it]/dt2+pxhxh/dt2);
	yxtoxy(nt,1.0,0.0,&ttn[0],nt,1.0,0.0,-nt,nt,&tnt[0]);
	ints8r(nt,1.0,0,&m[iqxnt],0.0,0.0,nt,tnt,dint);
      }
      if (adj) for (it=0;it<nt;it++) m[iqxnt+it]+=dint[it];
      else  for (it=0;it<nt;it++) d[ihxnt+it]+=dint[it];            
    }
  }
  
  free1float(dint);
  free1float(ttn);
  free1float(tnt);

  return;
}

void radonhyp_crude(float *m, float *t, float *h, float *q, float *d, int adj, int nt, int nh, int nq) 
{
  register int it;
  int ih,iq;
  unsigned int j;
  float hxh,pxhxh;
  int iqxnt,ihxnt,itint;
  int nx=nt*nq;
  int ny=nt*nh;
  float dt=t[1]-t[0];
  float dt2=dt*dt;
  int ttn;

  fprintf(stderr,"CRUDE APPROACH+++++++++++++++++++++++++++++++++++\n");

  if (adj)  for (it=0;it<nx;it++) m[it]=0;
  else for (it=0;it<ny;it++) d[it]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      pxhxh=hxh*q[iq];
      iqxnt=iq*nt;
      for (it=0;it<nt;it++){
	ttn=(int) (sqrt(t[it]*t[it]/dt2+pxhxh/dt2)+0.5);
	if (ttn<nt){
	  if (adj) m[iqxnt+it]+=d[ihxnt+ttn];
	  else d[ihxnt+ttn]+=m[iqxnt+it];
	}
      }
    }
  }

  return;
}


void radonpar(float *m, float *t, float *h, float *q, float *d, int adj, int nt, int nh, int nq) 
{
  register int it;
  int ih,iq;
  unsigned int j;
  float *ftime,*dint,*dtemp,time,hxh,pxhxh;
  int iqxnt,ihxnt,itint;
  int itime;
  int nx=nt*nq;
  int ny=nt*nh;

  float dt=t[1]-t[0];

  if ((ftime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");

  if (adj)  for (it=0;it<nx;it++) m[it]=0;
  else for (it=0;it<ny;it++) d[it]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      pxhxh=hxh*q[iq];
      iqxnt=iq*nt;
      for (it=0;it<nt;it++){
        if (adj){
	  time=t[it]+pxhxh;
	  ftime[it]=time/dt;
	  dtemp[it]=d[ihxnt+it];
        }
        else{
	  dtemp[it]=m[iqxnt+it];
	  time=t[it]-pxhxh;
          if (time>0) ftime[it]=time/dt;	  
   	  else ftime[it]=0;
	}
      }
      ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,ftime,dint);
      if (adj) for (it=0;it<nt;it++) m[iqxnt+it]+=dint[it];
      else  for (it=0;it<nt;it++) if (ftime[it]>0) d[ihxnt+it]+=dint[it];            
    }
  }
  
  free1float(dtemp);
  free1float(ftime);
  free1float(dint);
  return;
}


void radonlin(float *m, float *t, float *h, float *q, float *d, int adj, int nt, int nh, int nq) 
{
  register int it;
  int ih,iq;
  unsigned int j;
  float *ftime,*dint,*dtemp,time,hxh,pxhxh;
  int iqxnt,ihxnt,itint;
  int itime;
  int nx=nt*nq;
  int ny=nt*nh;

  float dt=t[1]-t[0];

  if ((ftime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");

  if (adj)  for (it=0;it<nx;it++) m[it]=0;
  else for (it=0;it<ny;it++) d[it]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      pxhxh=hxh*q[iq];
      iqxnt=iq*nt;
      for (it=0;it<nt;it++){
        if (adj){
	  time=t[it]+pxhxh;
	  ftime[it]=time/dt;
	  dtemp[it]=d[ihxnt+it];
        }
        else{
	  dtemp[it]=m[iqxnt+it];
	  time=t[it]-pxhxh;
          if (time>0) ftime[it]=time/dt;	  
   	  else ftime[it]=0;
	}
      }
      ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,ftime,dint);
      if (adj) for (it=0;it<nt;it++) m[iqxnt+it]+=dint[it];
      else  for (it=0;it<nt;it++) if (ftime[it]>0) d[ihxnt+it]+=dint[it];            
    }
  }
  
  free1float(dtemp);
  free1float(ftime);
  free1float(dint);
  return;
}


void radonhyp(float *m, float *t, float *h, float *q, float *d, float *vel,int adj, int nt, int nh, int nq) 
{
  register int it;
  int ih,iq;
  unsigned int j;
  float *ttn,*dint,*tnt,hxh,pxhxh;
  int iqxnt,ihxnt,itint;
  int nx=nt*nq;
  int ny=nt*nh;
  float dt=t[1]-t[0];
  float dt2=dt*dt;
  
  if ((ttn=alloc1float(nt))==NULL)
    err("cannot allocate memory for ttn \n");
  if ((tnt=alloc1float(nt))==NULL)
    err("cannot allocate memory for tnt \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  fprintf(stderr,"NEW APPROACH+++++++++++++++++++++++++++++++++++\n");

  if (adj)  for (it=0;it<nx;it++) m[it]=0;
  else for (it=0;it<ny;it++) d[it]=0;
  
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      iqxnt=iq*nt;
      for (it=0;it<nt;it++){
	pxhxh=hxh/((q[iq]+vel[it])*(q[iq]+vel[it]));
	ttn[it]=sqrt(t[it]*t[it]/dt2+pxhxh/dt2);
      }	
      if (adj) ints8r(nt,1.0,0,&d[ihxnt],0.0,0.0,nt,ttn,dint);
      else{
	yxtoxy(nt,1.0,0.0,ttn,nt,1.0,0.0,-nt,nt,tnt);
	ints8r(nt,1.0,0,&m[iqxnt],0.0,0.0,nt,tnt,dint);
      }
      if (adj) for (it=0;it<nt;it++) m[iqxnt+it]+=dint[it];
      else  for (it=0;it<nt;it++) d[ihxnt+it]+=dint[it];
    }            
  }
  
  free1float(dint);
  free1float(ttn);
  free1float(tnt);  
  return;
}

void radonhyp_old(float *m, float *t, float *h, float *q, float *d, float **vel,int adj, int nt, int nh, int nq) 
{
  register int it;
  int ih,iq;
  unsigned int j;
  float ftime;
  float dint;
  float time,hxh,pxhxh;
  int iqxnt,ihxnt,itint;
  int itime;
  float a1;
  float a2;
  int nx=nt*nq;
  int ny=nt*nh;

  float dt=t[1]-t[0];


  if (adj)  for (it=0;it<nx;it++) m[it]=0;
  else for (it=0;it<ny;it++) d[it]=0;
  
  for (it=0;it<nt;it++){
    for (ih=0;ih<nh;ih++){
      hxh=h[ih]*h[ih];
      ihxnt=ih*nt;
      for (iq=0;iq<nq;iq++){    
	pxhxh=hxh/vel[iq][it];
	iqxnt=iq*nt;
        if (adj){
	  time=sqrt(t[it]*t[it]+pxhxh);
	  ftime=time/dt;
	  ints8r(nt,1.0,0,&d[ihxnt],0.0,0.0,1,&ftime,&dint);
	  m[iqxnt+it]+=dint;
        }
        else{
	  time=sqrt(t[it]*t[it]+pxhxh);
	  ftime=time/dt;
	  itime=(int) floor(ftime);
          a2=ftime-itime;
	  a1=1-a2;	  
	  if (itime<nt)  d[ihxnt+itime]+=a1*m[iqxnt+it];
	  if ((itime+1) < nt) d[ihxnt+itime+1]+=a2*m[iqxnt+it];
	}
      }            
    }
  }
  
  return;
}



void radonhyp(float *m,float *t, float *h, float *q, float *d, float **vel,int adj,int nt, int nh, int nq)
{
  register unsigned int it;
  unsigned int ih,iq;
  float *ttn,*dint,*tnt,hxh,pxhxh;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];
  float dt2=dt*dt;

  if ((ttn=alloc1float(nt))==NULL)
    err("cannot allocate memory for ttn \n");
  if ((tnt=alloc1float(nt))==NULL)
    err("cannot allocate memory for tnt \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  fprintf(stderr,"NEW APPROACH+++++++++++++++++++++++++++++++++++\n");

  if (adj)  for (it=0;it<nx;it++) m[it]=0;
  else for (it=0;it<ny;it++) d[it]=0;   

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=0;it<nt;it++){
	pxhxh=hxh/vel[iq][it];
	ttn[it]=sqrt(t[it]*t[it]/dt2+pxhxh/dt2);
      }	
      if (adj) ints8r(nt,1.0,0,&d[ihxnt],0.0,0.0,nt,ttn,dint);
      else{
	//	for  (it=1;it<nt;it++) if (ttn[it]<=ttn[it-1]) fprintf(stderr,"ERRR:ttn[%d]=%f,ttn[%d]=%f,t[%d]=%f,t[%d]=%f\n",it-1,ttn[it-1],it,ttn[it],it-1,t[it-1],it,t[it]);
	yxtoxy(nt,1.0,0.0,ttn,nt,1.0,0.0,-nt,nt,tnt);
	ints8r(nt,1.0,0,&m[iqxnt],0.0,0.0,nt,tnt,dint);
      }
      if (adj) for (it=0;it<nt;it++) m[iqxnt+it]+=dint[it];
      else  for (it=0;it<nt;it++) if ((tnt[it]<nt)&&(tnt[it]>=0)) d[ihxnt+it]+=dint[it];
    }            
  }
  free1float(dint);
  free1float(ttn);
  free1float(tnt);  
  return;
}



void radonhyp_crude(float *m,float *t, float *h, float *q, float *d, float **vel,int adj,int nt, int nh, int nq)
{
  register unsigned int it;
  unsigned int ih,iq;
  float ftime;
  float dint;
  float time,hxh,pxhxh;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];
  fprintf(stderr,"Crude approach\n");
  if (adj)  for (it=0;it<nx;it++) m[it]=0;
  else for (it=0;it<ny;it++) d[it]=0;   

  for (it=0;it<nt;it++){
    for (ih=0;ih<nh;ih++){
      hxh=h[ih]*h[ih];
      ihxnt=ih*nt;
      for (iq=0;iq<nq;iq++){    
	pxhxh=hxh/vel[iq][it];
	iqxnt=iq*nt;
        time=sqrt(t[it]*t[it]+pxhxh);
	ftime=time/dt;
	itime=(int) floor(ftime+0.5);
	if (itime<nt){
	  if (adj) m[iqxnt+it]+=d[ihxnt+itime];
          else d[ihxnt+itime]+=m[iqxnt+it];
	}
      }            
    }
  }
  return;
}

void radonhyp_inv(float *m,float *t, float *h, float *q, float *d, float **vel,int adj,int nt, int nh, int nq)
{
  register unsigned int it;
  unsigned int ih,iq;
  float *ttn,*dint,*tnt,hxh,pxhxh;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];
  float dt2=dt*dt;

  if ((ttn=alloc1float(nt))==NULL)
    err("cannot allocate memory for ttn \n");
  if ((tnt=alloc1float(nt))==NULL)
    err("cannot allocate memory for tnt \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  fprintf(stderr,"NEW APPROACH+++++++++++++++++++++++++++++++++++\n");

  if (adj)  for (it=0;it<nx;it++) m[it]=0;
  else for (it=0;it<ny;it++) d[it]=0;   

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=0;it<nt;it++){
	pxhxh=hxh*vel[iq][it];
	ttn[it]=sqrt(t[it]*t[it]/dt2+pxhxh/dt2);
      }	
      if (adj) ints8r(nt,1.0,0,&d[ihxnt],0.0,0.0,nt,ttn,dint);
      else{
	yxtoxy(nt,1.0,0.0,ttn,nt,1.0,0.0,-nt,nt,tnt);
	ints8r(nt,1.0,0,&m[iqxnt],0.0,0.0,nt,tnt,dint);
      }
      if (adj) for (it=0;it<nt;it++) m[iqxnt+it]+=dint[it];
      else  for (it=0;it<nt;it++) d[ihxnt+it]+=dint[it];
    }            
  }
  free1float(dint);
  free1float(ttn);
  free1float(tnt);  
  return;
}

void radonhyp_crude_inv(float *m,float *t, float *h, float *q, float *d, float **vel,int adj,int nt, int nh, int nq)
{
  register unsigned int it;
  unsigned int ih,iq;
  float ftime;
  float dint;
  float time,hxh,pxhxh;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];
  fprintf(stderr,"Crude approach nh=%d, nt=%d, nq=%d\n",nh,nt,nq);

  if (adj) memset((void *) m,(int)'\0',nx*FSIZE);  
  //for (it=0;it<nx;it++) m[it]=0;
  else  memset((void *) d,(int)'\0',ny*FSIZE);  
  //for (it=0;it<ny;it++) d[it]=0;

  for (it=0;it<nt;it++){
    for (ih=0;ih<nh;ih++){
      hxh=h[ih]*h[ih];
      ihxnt=ih*nt;
      for (iq=0;iq<nq;iq++){    
	pxhxh=hxh*vel[iq][it];
	iqxnt=iq*nt;
        time=sqrt(t[it]*t[it]+pxhxh);
	ftime=time/dt;
	itime=(int) floor(ftime+0.5);
	if (itime<nt){
	  if (adj) m[iqxnt+it]+=d[ihxnt+itime];
          else d[ihxnt+itime]+=m[iqxnt+it];
	}
      }            
    }
  }
  return;
}



