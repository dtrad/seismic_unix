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
      vaux=perv[iq];//*(1.+0.0*t[it]);
      vgrid[iq][it]=vel[it]*vel[it]+vaux*vaux+2*vel[it]*vaux;
    }
    fprintf(stderr,"veltop=%f===>velbot=%f\n",sqrt(vgrid[iq][0]),sqrt(vgrid[iq][nt-1]));
    q[iq]=perv[iq];
  }
  //  for(iq=0;iq<nq;iq++) fprintf(stderr,"q[%d]=%f\n",iq,q[iq]);
  free1float(perv); 
  return;
}

void irreg_velaxis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid, int centralq)
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
  unsigned int nqh=centralq;
  float vaux;

  perv=ealloc1float(nq+1);
  perv[nqh]=0;
  for (iq=nqh;iq<nq;iq++) perv[iq+1]=perv[iq]*(1+dperv)+pervmin;
  for (iq=nqh-1;iq>=0;iq--) perv[iq]=perv[iq+1]*(1+dperv)-pervmin;
  for (iq=0;iq<nq;iq++){
    for (it=0;it<nt;it++){
      vaux=perv[iq];//*(1.+0.0*t[it]);
      vgrid[iq][it]=vel[it]*vel[it]+vaux*vaux+2*vel[it]*vaux;
    }
    q[iq]=perv[iq];
  }
  //  for(iq=0;iq<nq;iq++) fprintf(stderr,"q[%d]=%f\n",iq,q[iq]);
  free1float(perv); 
  return;
}

void irreg_slowness_axis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid, int centralq)
{
  /*   Velocity Grid 
   Irregular velocity grid
   Define irregular spacing along the horizontal
   Here q  is perturbartion around the central velocity law
   dperq is a parameter used to generate the increasing space
   perqmin defines the minimum distance to the perturbation */
     
  float *perv;
  int it;
  int iq;
  int nqh=centralq;
  float vaux;
  float qaux;
  float dqaux;
  int nqleft=0; // This is a experimental variable.
  // Used to generate the regular space for the artifacts
  // It must be set to the same value as in eomigls_sparse
  fprintf(stderr,"pervmin=%f,dperv=%f\n",pervmin,dperv);
  if ((perv=alloc1float(nq+1))==NULL) err(" Cannot allocate perv");
  perv[nqh]=0;
  for (iq=nqh;iq<nq;iq++) perv[iq+1]=perv[iq]*(1+dperv)+pervmin;
  for (iq=nqh-1;iq>=0;iq--) perv[iq]=perv[iq+1]*(1+dperv)-pervmin;
  for (iq=0;iq<nqleft;iq++){
    vaux=50000;
    qaux=1./(vaux*vaux);
    dqaux=1/(0.5*vaux*vaux);
    for (it=0;it<nt;it++) vgrid[iq][it]=qaux+iq*dqaux; 
    q[iq]=(vgrid[iq][0]); 
  }
  for (iq=nqleft;iq<nq;iq++){
    //fprintf(stderr,"perv[%d]=%f\n",iq,perv[iq]);
    vaux=perv[iq];

    for (it=0;it<nt;it++){
      //vgrid[iq][it]=1./(vel[it]*vel[it])+vaux; // --> Correct
      vgrid[iq][it]=1./(vel[nt/2]*vel[nt/2])+vaux;   // --> test
      if (vgrid[iq][it]<0) 
	fprintf(stderr,"vgrid[%d][%d]=%e\n",iq,it,vgrid[iq][it]);
	       //vgrid[iq][it]=0;//vgrid[iq][MAX(it-1,0)];
    }
    //    q[iq]=perv[iq];
    q[iq]=(vgrid[iq][0]);
  }
  for (iq=0;iq<nq;iq++)
    fprintf(stderr,"vtop[%d]=%6.0f<======>,vbot[%d]=%6.0f<======>%e\n",iq,sqrt(fabs(1./vgrid[iq][0])),iq,sqrt(fabs(1./vgrid[iq][nt-1])),q[iq]);
  
  //  for(iq=0;iq<nq;iq++) fprintf(stderr,"q[%d]=%f\n",iq,q[iq]);
  free1float(perv); 
  return;
}

void reg_slowness_axis(int nq, int nt, float dperv, float *t, float *vel,float *q, float **vgrid)
{
  /*   Velocity Grid 
   Irregular velocity grid
   Define irregular spacing along the horizontal
   Here q is perturbartion around the central velocity law
   and dq is the regular spaced perurbation */
     
  unsigned int it;
  int iq;
  float dq=dperv;

  for (iq=0;iq<nq;iq++){
    for (it=0;it<nt;it++) vgrid[iq][it]=1./(vel[it]*vel[it])+(iq-((int) nq/2))*dq;
    q[iq]=(vgrid[iq][0]);
    fprintf(stderr,"dq=%f,vel[10]=%f,q[%d]=%f\n",dq,vel[10],iq,q[iq]);
  }
  return;
}









