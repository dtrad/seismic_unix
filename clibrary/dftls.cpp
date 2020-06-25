#include "interpfk.h"


void dftls(complex **F, complex **FLST, int nh, int nk, float eps)
{
  int ik;
  complex **FF=0;
  complex **FT=0;
  complex **FLS=0;


  FT=ealloc2complex(nh,nk);
  FF=ealloc2complex(nk,nk);
  FLS=ealloc2complex(nh,nk);

  AHtimesB(FF,F,F,nh,nk,nh,nk);
  // regularization
  for (ik=0;ik<nk;ik++) FF[ik][ik]+=eps;
  transpose(F,FT,nh,nk);

  inverse_matrix_multiply(nk,FF,nh,nk,FT,FLS);  

  transpose(FLS,FLST,nk,nh);

  free2complex(FT);
  free2complex(FF);
  free2complex(FLS);

  return;
}

