#include <math.h>
#include "su.h"
void rjwfilter(float *d,int nt, float dt);
void kmig1(float *d, float cdp, float h,float **m, float *t, float *x, float vel) 
{
  int i,k,j,ih,ix;
  register int it;
  float *itime,*dint,*obliq,time,sx,gx,t02;
  float s2,g2,ivel2,xm, temp;
  extern int nt,nh,nx;
  extern float dt,dh,dx;
  
  ivel2=1.0/(vel*vel);
  xm=cdp;

  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((obliq=alloc1float(nt))==NULL)
    err("cannot allocate memory for obliq \n");
  for (ix=0;ix<nx;ix++){
    xm=fabs(cdp-x[ix]);
    sx=xm-h;
    gx=xm+h;
    s2=sx*sx*ivel2;
    g2=gx*gx*ivel2;
    // DSR equation
    for (it=0;it<nt;it++){
      t02=(t[it]*t[it]/4);
      time=sqrt(t02+s2)+sqrt(t02+g2);
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
    ints8r(nt,1.0,0.,d,0.0,0.0,nt,itime,dint);
    
    for (it=0;it<nt;it++)
      //m[it][ix]+=dint[it];        
      m[it][ix]+=(dint[it]*obliq[it]);                  
  }
  free1float(itime);
  free1float(dint);
  free1float(obliq);
  return;
}















