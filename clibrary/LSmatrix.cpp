#include "radonline_stack.h"

void LSmatrix(complex **F, complex **FLS, int nh, int nk, float eps)
{
  int ik;
  complex **FF=0;
  complex **FT=0;

  FT=ealloc2complex(nh,nk);
  FF=ealloc2complex(nk,nk);


  AHtimesB(FF,F,F,nh,nk,nh,nk);
  // regularization
  for (ik=0;ik<nk;ik++) FF[ik][ik]+=eps;
  transpose(F,FT,nh,nk);
  
  inverse_matrix_multiply(nk,FF,nh,nk,FT,FLS);  

  free2complex(FT);
  free2complex(FF);


  return;
}










