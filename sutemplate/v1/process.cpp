#include "Complex.h"
#include "su.h"

void fftgo0(int sign,float **d,complex  **m, int nh, int nt, float dt, int *nf0);
int fftback0(int sign,float **d,complex  **m, int nh, int nt, float dt, int nf0);
complex cdot(complex *x,complex *y,int n);
bool process(float** datain, float** dataout, int nt, int  nh, float dt){
  int it, ih;
  //  for (ih=0;ih<nh;ih++)  for (it=0;it<nt;it++) dataout[ih][it]=datain[ih][it];
  int nf0,ifreq,maxfreq, nfreq;
  float scale = 0;
  complex** model = ealloc2complex(nh,nt);
  float fmax=0, df;
  complex czero; 
  czero.r = czero.i=0;

  fftgo0(-1,datain ,model,nh,nt,dt,&nf0);
    
  nfreq=nf0/2;
  df=1/(nf0*dt);
  if (fmax==0) maxfreq=nfreq;
  else         maxfreq=(int) (fmax/df);

  for (ifreq=10;ifreq<maxfreq;ifreq++) 
    for (ih=0;ih<nh;ih++)
      model[ifreq][ih] = czero;
  

  fftback0(1,dataout ,model,nh,nt,dt,nf0);

  
  

  return true;
  
}
