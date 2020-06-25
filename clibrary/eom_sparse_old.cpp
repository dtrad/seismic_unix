#include <math.h>
#include "su.h"

void irreg_velaxis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid)
{
  /*   Velocity Grid 
   Irregular velocity grid
   Define irregular spacing along the horizontal
   Here q  is perturbartion around the central velocity law
   dperq is a parameter used to generate the increasing space
   perqmin defines the minimum distance to the perturbation */
     
  float *perv;
  unsigned int it;
  int iq;
  unsigned int nqh=nq/2;
  float vaux;

  if ((perv=alloc1float(nq+1))==NULL) err(" Cannot allocate perv");
  perv[nqh]=0;
  for (iq=nqh;iq<nq;iq++) perv[iq+1]=perv[iq]*(1+dperv)+pervmin;
  for (iq=nqh-1;iq>=0;iq--) perv[iq]=perv[iq+1]*(1+dperv)-pervmin;
  for (iq=0;iq<nq;iq++){
    for (it=0;it<nt;it++){
      vaux=perv[iq]*(1.+0.0*t[it]);
      vgrid[iq][it]=vel[it]*vel[it]+vaux*vaux+2*vel[it]*vaux;
      q[iq]=perv[iq];
    }
  }
  //  for(iq=0;iq<nq;iq++) fprintf(stderr,"q[%d]=%f\n",iq,q[iq]);
  free1float(perv); 
  return;
}

int count_nsparse(float *t, float *h, float *q, float **vgrid, int nt, int nh, 
		  int nq, float t0)
{
  unsigned int it;
  unsigned int ih,iq;
  unsigned int j;
  unsigned int nonzero; 
  float time,ftime,hxh,moveout,slowness2;
  unsigned int iqxnt,ihxnt,itint;
  unsigned int itime;
  unsigned int nx=nt*nq;
  unsigned int ny=nt*nh;
  float dt=t[1]-t[0];
  int it0=(int) (t0/dt+0.5);
  // First we need to find out how big the sparse matrix will be
  j=0;

  for (it=it0;it<nt;it++){
    for (iq=0;iq<nq;iq++){
      slowness2=1./(vgrid[iq][it]);
      iqxnt=iq*nt;
      for (ih=0;ih<nh;ih++){
	hxh=h[ih]*h[ih];
	ihxnt=ih*nt;
	moveout=hxh*slowness2;
        time=2*sqrt(t[it]*t[it]/4.+moveout);
	ftime=time/dt;
	itime=(int) floor(ftime+0.5);
	if (itime<nt) j++;
      }            
    }
  }
  nonzero=j;
  fprintf(stderr,"nonzero=%d\n",nonzero);
  return(nonzero);
}

void build_index(float *t, float *h, float *q, float **vgrid, int nt, int nh, 
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

  for (j=0;j<2;j++) for (it=0;it<nonzero;it++) index[j][it]=0;
  j=0;
  for (it=it0;it<nt;it++){
      for (iq=0;iq<nq;iq++){    
	slowness2=1./(vgrid[iq][it]);
	iqxnt=iq*nt;
	for (ih=0;ih<nh;ih++){
	  hxh=h[ih]*h[ih];
	  ihxnt=ih*nt;
	  moveout=hxh*slowness2;
	  time=2*sqrt(t[it]*t[it]/4.0+moveout);
	  ftime=time/dt;
	  itime=(int) floor(ftime+0.5);
	  if (itime<nt){
	    index[0][j]=ihxnt+itime;
	    index[1][j]=iqxnt+it;
	    j++;
	  }
	}            
      }
  }
  return;
}

void eomig(float *m, float *d, unsigned short **index, int adj, int nt, int nh, int nq, int nsparse)
{
  unsigned long j;
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
  
  return;
}

















