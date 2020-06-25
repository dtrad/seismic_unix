#include "interpfk.h"


void FTmatrix(complex **F, complex **FLS, float *h, float *k, int nh, int nk, float eps)
{

  zero_array(F,nh,nk);
  dft_matrix(F,h,k,nh,nk);
  dftls(F,FLS,nh,nk,eps); 
  //AequalB(FLS,F,nk,nh);

  return;
}



void dft_matrix(complex **F,float *h,float *k,int nh,int nk)
{

  /******************************************************************** 
     FFT Transformation matrices. F, FH and R= top row of LH*L
     Input parameters:
     Out parameter:
     Notes:
     Daniel Trad- 22-02-99
  *********************************************************************/
   
     int ih, ik;
     complex  co;
     complex  dco;
     complex phase, dphase;
     float dk=k[1]-k[0];


     for (ih=0;ih<nh;ih++){
       phase.r=dphase.r=0;
       phase.i=((h[ih]-h[0])*(k[0]-dk));
       dphase.i=((h[ih]-h[0])*dk);
       co=exp(phase);
       dco=exp(dphase);
       
       for (ik=0;ik<nk;ik++){
	 co*=dco;
	 F[ih][ik]=(co);
       }
     }
  	      
     return;
}

void Atimesx_DFT(complex *y,complex **A,complex *x,int ny,int nx, int adj)
{
/*
  From Claerbout's PVI
  subroutine matmult( adj, bb,nx, x, ny, y)
  integer ix, iy, adj, nx, ny
  real bb(ny,nx), x(nx), y(ny)
  call adjnull(adj, 0,x,nx,  y,ny)
  do ix= 1, nx {
  do iy= 1, ny {
  if( adj == 0 ) y(iy) = y(iy) + bb(iy,ix) * x(ix)
  else x(ix) = x(ix) + bb(iy,ix) * y(iy)
  }}
  return; end

*/

  int ix,iy;
  
  if (!adj) for (iy=0;iy<ny;iy++) y[iy].r=y[iy].i=0.0;
  else for (ix=0;ix<nx;ix++) x[ix].r=x[ix].i=0.0;
     
  for (ix=0;ix<nx;ix++)
    for (iy=0;iy<ny;iy++)
      if (!adj) y[iy]+=A[iy][ix]*x[ix];
      else x[ix]+=conjg(A[iy][ix])*y[iy];
  
  return;
}








