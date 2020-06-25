#include <math.h>
#include "su.h"

void kmig0(float **d, float **m, float *t, float *h, float *x, float vel, float cdp) 
{
  int i,k,j,it,ih,ix,itau;
  float *itime,*dint,*dtemp,time,sx,gx,t02;
  float s2,g2,ivel2,xm;
  extern int nt,nh,nx;
  extern float dt,dh,dx;

  ivel2=1.0/(vel*vel);
  xm=cdp;

  if ((itime=alloc1float(nt))==NULL)
    err("cannot allocate memory for itime \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  if ((dtemp=alloc1float(nt))==NULL)
    err("cannot allocate memory for dtemp \n");

  for (ix=0;ix<nx;ix++){
    xm=fabs(cdp-x[ix]);
    //fprintf(stderr,"xm=%f,cdp=%f,x[%d]=%f\n",xm,cdp,ix,x[ix]);  
    for (ih=0;ih<nh;ih++){
      sx=xm-h[ih];
      gx=xm+h[ih];
      s2=sx*sx*ivel2;
      g2=gx*gx*ivel2;
      //fprintf(stderr,"xm=%f,s2=%f;g2=%f,vel=%f,ivel2=%f12.8\n",xm,s2,g2,vel,ivel2);   
      for (itau=0;itau<nt;itau++){
        // DSR equation
        t02=(t[itau]*t[itau]);
	time=sqrt(t02+s2)+sqrt(t02+g2);
       	itime[itau]=time/dt;
	//fprintf(stderr,"sx=%f,gx=%f,time=%f,itime[%d]=%f\n",sx,gx,time,itau,itime[itau]);
        dtemp[itau]=d[itau][ih];
        //fprintf(stderr,"dtemp[%d]=%f\n",itau,dtemp[itau]);
        //it=(int) (time/dt); if (it<nt) m[itau][ix]=m[itau][ix]+d[it][ih];
      }
      ints8r(nt,1.0,0.,dtemp,0.0,0.0,nt,itime,dint);
      //for (itau=0;itau<nt;itau++){
	//if (dtemp[itau] > 1e-1) fprintf(stderr,"dtemp[%d]=%f\n",itau,dtemp[itau]);
        //fprintf(stderr,"itime[%d]=%f\n",itau,itime[itau]);
	//if (dint[itau] > 1e-3) fprintf(stderr,"dint[%d]=%f\n",itau,dint[itau]);
      
      //}
      for (itau=0;itau<nt;itau++)        
	m[itau][ix]+=dint[itau];                  
    }
  }
  free1float(dtemp);
  free1float(itime);
  free1float(dint);
  return;
}















