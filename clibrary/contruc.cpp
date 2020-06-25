#include "su.h"
void advance(int conj,int sum,int jump,int nx,float *xx,int ny, float *yy);
void contran(int conj,int sum, int nx, float *xx, int nb, float *bb,float *yy);
void conjzero(int conj, int add, int nx, float *x, int ny, float *y);

void contruc(int conj,int add, int lag, int nx, float *xx, int nb, float *bb,int ny,float *yy)
{
  // if sum==0 erases the output first
  int ib;
  int ix;
  int ns=nx+nb-1;
  float *ss;

  if ((ss = alloc1float(ns)) == NULL)
    err("cannot allocate memory for ss\n");

  if (conj == 0){ 
    contran(0,0,nx,xx,nb,bb,ss);
    advance(0,add,lag,ns,ss,ny,yy);
  }
  else{
    advance(1,0,lag,ns,ss,ny,yy);
    contran(1,0,nx,xx,nb,bb,ss);
    
  }

  return;
}

void advance(int conj,int sum,int jump,int nx,float *xx,int ny, float *yy)
{
  int ix;
  int iy;
  conjzero(conj,sum,nx,xx,ny,yy);
  for (iy=0;iy<ny;iy++){
    ix=iy+jump;
    if (ix >=0 )
      if (ix < nx)
	if (conj==0)
	  yy[iy]+=xx[ix];
	else
	  xx[ix]+=yy[iy];
  }
  return;
}


