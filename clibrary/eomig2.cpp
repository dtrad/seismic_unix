#include <math.h>
#include "su.h"

void eomig1(float **csp, float *he,float *m, float *t, float *x, float *vel) 
{
  int i,k,j,ih,ix;
  register int it;
  float *itime,*dint,*dtemp,*obliq,time,sx,gx,t02;
  float temp, moveout;
  extern int nt,nh,nx,nhmax;
  extern float dt,dh,dx;
  
  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((obliq=alloc1float(nt))==NULL)
    err("cannot allocate memory for obliq \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");
 
  for (ih=0;ih<nhmax;ih++){

    // DSR equation
    for (it=0;it<nt;it++){
      moveout=he[ih]*he[ih]/(vel[it]*vel[it]);
      dtemp[it]=csp[it][ih];
      t02=(t[it]*t[it]/4);
      time=2*sqrt(t02+moveout);
      itime[it]=time/dt;
      //fprintf(stderr,"pow=%e\n",pow((t[it]/time),2));
      if (time>1e-2){
	temp=t[it]/time;//temp=pow(temp,1.5);
        temp=sqrt(temp*temp*temp);
      	obliq[it]=temp;
      }
      else obliq[it]=0;
      if (obliq[it]>1) fprintf(stderr,"obliq[%d]=%e\n",it,obliq[it]);
    }
        
    ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,itime,dint);
    
    for (it=0;it<nt;it++)
      m[it]+=(dint[it]*obliq[it]);                  
  }
  free1float(itime);
  free1float(dtemp);
  free1float(dint);
  free1float(obliq);
  return;
}



















