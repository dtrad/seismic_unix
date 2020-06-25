#include "cwp.h"
#include "dan.h"

void plotAmpSpec(float **pk, int nk, int nt, float dt)

{
  int nw;
  float dt;
  float **aspec; // Amplitude spectrum
  aspec=ealloc2float(nw,nk);

  nw = (int) (nt*2/0.6);
  nw = npfao(nw,nw*2);
  dw = 2.0*PI/(nw*dt);
  fw = -PI/dt;

  pp = alloc2complex(nw);
  
  for (ik=0;ik<nk;ik++){
    for (it=0; it<nt; it+=2)
      pp[it] = p[ik][it];
    for (it=1; it<nt; it+=2) {
      pp[it].r = -p[ik][it].r;
      pp[it].i = -p[ik][it].i;
    }
    for (it=nt; it<nw; ++it)
    pp[it].r = pp[it].i = 0.0;
    pfacc(1,nw,pp);
    /* zero -Nyquist frequency for symmetry */
    pp[0] = czero;
    for (it=0; it<nt; it++) pfk[ik][it]=abs(pp[it]);
      
  }
  if (1){
    save_gather(aspec,nk,nw,dt,"temp");
    system("suxwigb < temp perc=100 title=temp "); 
  }

  free2float(aspec);

}

