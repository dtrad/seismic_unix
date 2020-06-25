#include "su.h"
#include "segy.h"
#include "Complex.h"

void Atimesx(complex *b,complex **A,complex *x,int nr,int nc)
{
     int i,j;

     for (i=0;i<nr;i++){
         b[i].r=b[i].i=0.0;
         for (j=0;j<nc;j++)
	   b[i]+=(A[i][j]*x[j]);
     }
     return;
}

void xtimesA(complex *b,complex *x, complex **A, int nr, int nc)
{
     int i,j;
     for (j=0;j<nc;j++){
       b[j].r=b[j].i=0.0;
       for (i=0;i<nr;i++)
	 b[j]+=(x[i]*A[i][j]);
     }
     return;

}

void xtimesA(complex *b, float *x, complex **A, int nr, int nc)
{
     int i,j;
     for (j=0;j<nc;j++){
       b[j].r=b[j].i=0.0;
       for (i=0;i<nr;i++)
	 b[j]+=(x[i]*A[i][j]);
     }
     return;
}


void Atimesx(float *b,float **A,float *x,int nr,int nc)
{
     int i,j;

     for (i=0;i<nr;i++){
         b[i]=0;

         for (j=0;j<nc;j++)
	   b[i]+=(A[i][j]*x[j]);
     }
     return;
}

void Atimesx(int nr,int nc,float *b,float **A,float *x)
{
     int i,j;

     for (i=1;i<=nr;i++){
         b[i]=0;
         for (j=1;j<=nc;j++)
	   b[i]+=(A[i][j]*x[j]);
     }
     return;
}

void Atimesx(float *y,float **A,float *x,int adj, int ny,int nx)
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
     
     if (!adj) for (iy=0;iy<ny;iy++) y[iy]=0.0;
     else for (ix=0;ix<nx;ix++) x[ix]=0.0;
     
     for (ix=0;ix<nx;ix++)
       for (iy=0;iy<ny;iy++)
	 if (!adj) y[iy]+=(A[iy][ix]*x[ix]);
         else x[ix]+=(A[iy][ix]*y[iy]);
     
     return;
}

void Atimesx2(float *y,float **A,float *x, int ny,int nx, int adj)
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
     
     if (!adj) for (iy=0;iy<ny;iy++) y[iy]=0.0;
     else for (ix=0;ix<nx;ix++) x[ix]=0.0;
     
     for (ix=0;ix<nx;ix++)
       for (iy=0;iy<ny;iy++)
	 if (!adj) y[iy]+=(A[iy][ix]*x[ix]);
         else x[ix]+=(A[iy][ix]*y[iy]);
     
     return;
}


void Atimesxnr(float *y,float **A,float *x, int ny,int nx, int adj)
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
     
     if (!adj) for (iy=1;iy<=ny;iy++) y[iy]=0.0;
     else for (ix=1;ix<=nx;ix++) x[ix]=0.0;
     
     for (ix=1;ix<nx;ix++)
       for (iy=1;iy<ny;iy++)
	 if (!adj) y[iy]+=(A[iy][ix]*x[ix]);
         else x[ix]+=(A[iy][ix]*y[iy]);
     
     return;
}



void Atimesx(complex *y,complex **A,complex *x,int ny,int nx, int adj)
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
	 if (!adj) y[iy]+=(A[iy][ix]*x[ix]);
         else x[ix]+=(conjg(A[iy][ix])*y[iy]);
     
     return;
}

void Atimesx(complex **A, float x, int n1, int n2)
{
  int i1, i2;
  for (i2=0;i2<n2;i2++)
    for (i1=0;i1<n1;i1++)
      A[i2][i1]*=x;
  return;

}
void Atimesx(float **A, float x, int n1, int n2)
{
  int i1, i2;
  for (i2=0;i2<n2;i2++)
    for (i1=0;i1<n1;i1++)
      A[i2][i1]*=x;
  return;

}

