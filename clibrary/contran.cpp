/* Claerbout routines: processing versus inversion, 1992 */
#include "su.h"
float dot(int n, float *a, float *b);
void conjzero(int conj, int add, int nx, float *x, int ny, float *y);
void advance(int conj,int sum,int jump,int nx,float *xx,int ny, float *yy);
void contran(int conj,int sum, int nx, float *xx, int nb, float *bb,float *yy);
void conjzero(int conj, int add, int nx, float *x, int ny, float *y);
void contruc(int conj,int add, int lag, int nx, float *xx, int nb, float *bb,int ny,float *yy);
void convin(int conj,int sum, int nx, float *xx, int nb, float *bb,float *yy);
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)

/* Convolution and correlation: from Claerbout, 1992, Processing vs inversion 
   conj == 0 means convolution
   conj == 1 means correlation 
   add==0 erases the output first
   add==1 do not
*/
 

void contruc_2(int conj,int add, int nx, float *xx, int nb, float *bb,float *yy)
{
  // if add==0 erases the output first
  int lag=nx/2;
  int ny=nb;
  int ns=nx+nb-1;
  float *ss;
  //fprintf(stderr,"contruc_2+++++++++++");
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
  free1float(ss);
  return;
}



void contran(int conj,int sum, int nx, float *xx, int nb, float *bb,float *yy)
{
  // if sum==0 erases the output first
  int ny;
  int ib;
  int ix;
  //int lag=nx/2+1;
  
  ny=nx + nb -1;
  conjzero(conj,sum,nb,bb,ny,yy);
  if (conj == 0){ 
    for(ib=0;ib<nb;ib++) 
      for(ix=0;ix<nx;ix++) 
	yy[ib+ix]+=bb[ib]*xx[ix];
    //for (iy=0;iy<ny-lag;iy++) yy[iy]=y[iy+lag];
    //for (iy=ny-lag;iy<ny;iy++) yy[iy]=0;
  }
  else
    for(ib=0;ib<nb;ib++) 
      for(ix=0;ix<nx;ix++) 
	bb[ib]+=yy[ib+ix]*xx[ix];

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
	  yy[iy]=xx[ix];
	else
	  xx[ix]=yy[iy];
  }
  return;
}


void conjzero(int conj, int add, int nx, float *x, int ny, float *y) 
{
  int ix;
  int iy;
  if (add == 0){
    if(conj == 0){
      for (iy=0;iy<ny;iy++)
	y[iy]=0.;
    }
    else{
      for (ix=0;ix<nx;ix++)
	x[ix]=0;
    }
  }
  return;
}



void contruc(int conj,int add, int lag, int nx, float *xx, int nb, float *bb,int ny,float *yy)
{
  // if sum==0 erases the output first
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



void convin(int conj,int sum,int nb, float *bb, int nx, float *xx,float *yy)
{
  // if sum==0 erases the output first
  int ny;
  int ix;
  int iy;

  ny = nx - nb + 1;
  if (ny<1) err("convin: filter output has negative length");
  conjzero(conj,sum,nx,xx,ny,yy);
  if (conj == 0) 
    for(iy=0;iy<ny;iy++) 
      for(ix=0;ix<nx;ix++) 
	yy[iy]+=xx[ix]*bb[iy-ix+nx];
  else
    for(ix=0;ix<nx;ix++) 
      for(iy=0;iy<ny;iy++) 
	xx[ix]+=yy[iy]*bb[iy-ix+nb];
  return;
}



float testadjop_conv(void (*oper) (int, int, int, float *,int, float *,float *), 
		     int nb, float *bb, int nx, int ny)
{
  float *x1;
  float *y1;
  float *x2;
  float *y2;
  float dp1;
  float dp2;
  int it;
  float test;
  
  x1=ealloc1float(nx);
  y1=ealloc1float(ny);
  x2=ealloc1float(nx);
  y2=ealloc1float(ny);
 
  for (it=0;it<ny;it++) y1[it]=frannor();
  for (it=0;it<nx;it++) x1[it]=frannor();

  oper(0,0,nb,bb,nx,x1,y2);
  oper(1,0,nb,bb,nx,x2,y1); 

  dp1=dot(ny,y1,y2);
  dp2=dot(nx,x1,x2);

  if (dp2!=0) test=dp1/dp2;
  else test=0;

  fprintf(stderr,"++++Test adjoint = %f \n",test);
  return(test);
  
  free1float(y2);
  free1float(x2);
  free1float(y1);
  free1float(x1);

}







