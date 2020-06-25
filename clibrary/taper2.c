#include "su.h"

int taper(float *data, int nh, int ntaper,int flag)
{
  float *taper;
  int k, ih;
  float s;

  taper=ealloc1float(ntaper);

  for (k=0;k<ntaper;++k){ 
    s=sin(k*PI/(2*ntaper));
    taper[k]=s*s;
  }  
  taper[0]=1e-5; 
  /* Taper at the left end of the data set */
  for (ih = 0; ih < ntaper; ++ih) data[ih] /= taper[ih];

  /* Taper at the right end of the data set */
  for (ih=nh - ntaper ; ih < nh; ++ih)  data[ih] /= taper[nh - ih - 1];

  return EXIT_SUCCESS;

}
