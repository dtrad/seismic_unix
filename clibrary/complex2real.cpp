#include "su.h"
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)
/* Given a complex system of equations 
            b = A x
with A a matrix of size n x m
this function returns a the matrices for the equivalent real system of equations 
            br= AR x
with twice the size
references : Numerical recipes, pag 50, Chapter 2, 
             Solution of linear system of equations

Daniel Trad - UBC October 2000
*/

void complex2real(complex *b, complex **A, complex *m, float *br, float **AR, float *mr, int nb, int nx)
{
  int ix, ib;

  for (ib=0;ib<nb;ib++){
    br[ib]=b[ib].r;
    br[nb+ib]=b[ib].i;
  } 


  for (ix=0;ix<nx;ix++){
    mr[ix]=m[ix].r;
    mr[nx+ix]=m[ix].i;
  } 


  for (ib=0;ib<nb;ib++) 
    for (ix=0;ix<nx;ix++){ 
      AR[ib][ix]=A[ib][ix].r;
      AR[nb+ib][ix]=A[ib][ix].i;
      AR[ib][nx+ix]=-A[ib][ix].i;
      AR[nb+ib][nx+ix]=A[ib][ix].r;
    }
  return;
}

void complex2real(complex *b, complex **A, float *br, float **AR, int nb, int nx)
{
  int ix, ib;

  for (ib=0;ib<nb;ib++){
    br[ib]=b[ib].r;
    br[nb+ib]=b[ib].i;
  } 

  for (ib=0;ib<nb;ib++) 
    for (ix=0;ix<nx;ix++){ 
      AR[ib][ix]=A[ib][ix].r;
      AR[nb+ib][ix]=A[ib][ix].i;
      AR[ib][nx+ix]=-A[ib][ix].i;
      AR[nb+ib][nx+ix]=A[ib][ix].r;
    }
  return;
}

void real2complex(complex *b, float *br, complex *x, float *xr, int nb, int nx)
{
  int ib, ix;
  for (ib=0;ib<nb;ib++){
    b[ib].r=br[ib];
    b[ib].i=br[nb+ib];
  }

  for (ix=0;ix<nx;ix++){
    x[ix].r=xr[ix];
    x[ix].i=xr[nx+ix];
  }

}

void complex2realnr(complex *b, complex **A, float *br, float **AR, int nb, int nx)
{
  int ix, ib;

  for (ib=1;ib<=nb;ib++){
    br[ib]=b[ib-1].r;
    br[nb+ib]=b[ib-1].i;
  } 
  for (ib=1;ib<=nb;ib++) 
    for (ix=1;ix<=nx;ix++){ 
      AR[ib][ix]=A[ib-1][ix-1].r;
      AR[nb+ib][ix]=A[ib-1][ix-1].i;
      AR[ib][nx+ix]=-A[ib-1][ix-1].i;
      AR[nb+ib][nx+ix]=A[ib-1][ix-1].r;
    }
  return;
}

void real2complexnr(complex *b, float *br, complex *x, float *xr, int nb, int nx)
{
  int ib, ix;
  for (ib=1;ib<=nb;ib++){
    b[ib-1].r=br[ib];
    b[ib-1].i=br[nb+ib];
  }

  for (ix=1;ix<=nx;ix++){
    x[ix-1].r=xr[ix];
    x[ix-1].i=xr[nx+ix];
  }

}


