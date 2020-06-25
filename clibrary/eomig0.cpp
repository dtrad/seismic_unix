#include <math.h>
#include "su.h"
void rjwfilter(float *d,int nt, float dt);
void eomig0(float **csp, float *he,float *m, float *t, float *x, float vel) 
{
  int i,k,j,ih,ix;
  register int it;
  float *itime,*dint,*dtemp,*obliq,time,sx,gx,t02;
  float ivel2, temp, moveout;
  extern int nt,nh,nx,nhmax;
  extern float dt,dh,dx;
  
  ivel2=1.0/(vel*vel);

  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((obliq=alloc1float(nt))==NULL)
    err("cannot allocate memory for obliq \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");
 
  for (ih=0;ih<nhmax;ih++){
    moveout=he[ih]*he[ih]*ivel2;

    // DSR equation
    for (it=0;it<nt;it++){
      dtemp[it]=csp[it][ih];
      t02=(t[it]*t[it]/4);
      time=2*sqrt(t02+moveout);
      itime[it]=time/dt;
      //fprintf(stderr,"pow=%e\n",pow((t[it]/time),2));
      if ((t[it]>1e-3)&&(time>1e-3)){
	temp=t[it]/time;//temp=pow(temp,1.5);
        temp=sqrt(temp*temp*temp);
      	obliq[it]=temp;
      }
      else obliq[it]=0;
      //fprintf(stderr,"obliq[%d]=%e\n",it,obliq[it]);
    }
    
    //rjwfilter(d,nt,dt);    
    ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,itime,dint);
    
    for (it=0;it<nt;it++)
      //m[it]+=dint[it];        
      m[it]+=(dint[it]*obliq[it]);                  
  }
  free1float(itime);
  free1float(dtemp);
  free1float(dint);
  free1float(obliq);
  return;
}



















