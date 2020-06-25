#include <math.h>
#include "su.h"

/* Kirchhoff Migration for prestack data */

void kmig(float *d, float cdp, float h,float **m, float *t, float *x, float **ivel2, int nt, int nx, float dt, float cdpspace, float aper, int adj) 
{
  
  /* ivel2 is  1/vel^2 */
  
  int ix;
  register int it;
  int itime;
  float time,sx,gx,t02;
  float s2,g2,s2v,g2v,xm;
  float tsource;
  float trec;
  
  for (ix=0;ix<nx;ix++){
    xm=fabs(cdp-x[ix]);
    if (xm < aper){
      xm*=cdpspace;
      sx=xm-h;
      gx=xm+h;
      s2=sx*sx;
      g2=gx*gx;
      // DSR equation
      for (it=0;it<nt;it++){
	t02=(t[it]*t[it])*0.25;
	
	s2v=s2*ivel2[ix][it];
	g2v=g2*ivel2[ix][it];
	
	tsource=sqrt(t02+s2v);
	trec=sqrt(t02+g2v);
	time=tsource+trec;
	
	itime=(int) (time/dt+0.5);
	if (itime < nt){
	  if (adj)
	    m[ix][it]+=d[itime];
	  else
	    d[ix][itime]+=m[it];
	}
      }                  
    }
  }
  return;
}





















