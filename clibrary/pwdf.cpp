#include "su.h"
#include "radonwavelet.h"

#define GO 1
#define BACK -1

/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void pwdf0(float **data, float *t, int nt, float *h, int nh, float **model, float *q, int nq)
{
  int iw;
  float dt=t[1]-t[0];
  float t0=t[0];
  char buf[80];
  int nx=nh;   
  float **wav; // 2D wavelet 
  float sx; // scale in x direction
  float st; // scale in t direction
  float tau; // apix along t (tau)

  sx=2;
  st=2;
  tau=1;

  p=ealloc2float(nx,nt);
  hpwavelet(wav,sx,st,tau,nt,nx,dt,dx);

  plotgather(wav,nh,nt,"wavelet");
   

  free2float(p);

  return;

}

void hpwavelet(float **wav,float sx, float st,float tau,int nt,int nx,float dt,float dx)
{
  float tmp,x,t;
  zero_array(wav,nx,nt);

  for(ix=0;ix<nx;ix++){
    x=(ix-nx/2)*dx;
    for(it=0;it<nt;it++){
      t=it*dt;
      tmp=pow(((t-sqrt((x/sx)*(x/sx)+tau*tau))/st),2);
      wav[it][ix]=((1-tmp)*exp(-tmp/2));
    }
  }
  return;
}






