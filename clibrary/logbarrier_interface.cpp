#include "su.h"
#include "radonlogbar.h"

float logbarrier_interface(complex *d,complex **L,complex *m, int nh, int nq, float ftol, int itercg,float eps2)
{
  float *mr=0;
  float *dr=0;
  float **LR=0;
  int nq2;
  int nh2;
  float J=0;
  int iter_end=1;
  nh2=2*nh;
  nq2=2*nq;

  // Transform the complex system of equation in a real one 
  dr=ealloc1float(nh2);
  mr=ealloc1float(nq2);
  LR=ealloc2float(nq2,nh2);

  complex2real(d,L,dr,LR,nh,nq);

  //Atimesx2(dr,LR,mr,nh2,nq2,1);  

  if (1) radonl1freq_karmarkar(dr,LR,mr,nq2,nh2,iter_end,itercg);  

  //Atimesx2(dr,LR,mr,nh2,nq2,0);
  real2complex(d,dr,m,mr,nh,nq);  

  free2float(LR);
  free1float(dr);
  free1float(mr);

  return(J);

}










































