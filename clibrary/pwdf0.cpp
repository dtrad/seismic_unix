#include "su.h"
#include "pwdf.h"

#define GO 1
#define BACK -1

/*
radontoepf
Input is a hyperbolic gather, for example a CSP or CMP
Output is the same shape after Radon multiple removal.

Daniel Trad - June 9- 2000
*/


void pwdf0(float **dataxt, float *t, int nt, float *h, int nh, float **model, float *q, int nq)
{
  float dt=t[1]-t[0];
  float dx=fabs(h[1]-h[0]);
  //float t0=t[0];
  //char buf[80];
  int nx=nh;   
  
  float sx; // scale in x direction
  float st; // scale in t direction
  float tau; // apix along t (tau)
  int plot=1;

  float **datatotwav;
  complex **datakf;

  float **wavxt; // 2D wavelet 
  complex **wavkf;

  float **datawavxt;
  complex **datawavkf;


  float *w;
  float *k;

  int nw;
  int nk;
  
  float dw;
  float dk;

  float sfft;

  /* determine size factors for prime-factor FFTs */
  fksize(nt,nx,&nw,&nk);
  sfft=1./(2*(nw-1)*nk);



  wavxt=ealloc2float(nt,nx);
  wavkf=alloc2complex(nw,nk);

  datatotwav=alloc2float(nt,nx);  
  datakf=alloc2complex(nw,nk);

  datawavxt=alloc2float(nt,nx);  
  datawavkf=alloc2complex(nw,nk);  
  
  w=ealloc1float(nw);
  k=ealloc1float(nk);


  tau=1;


  zero_array(datatotwav,nh,nt);

  for (sx=600;sx<=3000;sx+=600){

    for (st=0.005;st<=0.015;st+=0.005){

      hpwavelet(wavxt,sx,st,tau,nt,nx,dt,dx);
      if (plot) plotgather(wavxt[0],nt,nh,"wavelet");

      fft2_parameters(nt,dt,nx,dx,&nw,&dw,&nk,&dk,w,k);

      /* data(x,t) -> data(k,f)*/ 
      fft2_xt2kf(dataxt, datakf,nt,nx,nw,nk);
      /* wav(x,t -> wav(k,f) */
      fft2_xt2kf(wavxt,wavkf,nt,nx,nw,nk);

      /* Dot product to find the wavelet coefficients (f,w) */
      Atimes_conjB_elem(datawavkf,datakf,wavkf,nw,nk);
      fft2_kf2xt(datawavxt, datawavkf,nt,nx,nw,nk);
      Atimesx(datawavxt,sfft*sfft,nt,nh);      
      
      // Thresolding 
      
      

      /* Wavelet expansion d_s= (f,w) w  */
      fft2_xt2kf(datawavxt, datawavkf,nt,nx,nw,nk);
      Atimes_B_elem(datawavkf,datawavkf,wavkf,nw,nk);

      Atimesx(datawavxt,sfft*sfft,nt,nh);      
      fft2_kf2xt(datawavxt, datawavkf,nt,nx,nw,nk);
      Atimesx(datawavxt,sfft*sfft,nt,nh);      
      /* Add the scale into the wavelet expansion d=sum(d_s) */ 
      Aplus_equalB(datatotwav,datawavxt,nt,nh);
    }
  }

  AequalB(dataxt,datatotwav,nt,nh);

      

  free1float(k);
  free1float(w);

  free2float(wavxt);
  free2complex(wavkf);

  free2float(datatotwav);
  free2complex(datakf);

  free2float(datawavxt);
  free2complex(datawavkf);


  return;

}

void hpwavelet(float **wav,float sx, float st,float tau,int nt,int nx,float dt,float dx)
{
  float tmp,x,t;
  int it, ix;
  zero_array(wav,nx,nt);

  for(ix=0;ix<nx;ix++){
    x=(ix-nx/2)*dx;
    for(it=0;it<nt;it++){
      t=it*dt;
      tmp=pow(((t-sqrt((x/sx)*(x/sx)+tau*tau))/st),2);
      wav[ix][it]=((1-tmp)*exp(-tmp/2))/(sx*st);
    }
  }
  return;
}






