#include "su.h"
#include "radonwavelet.h"

/*
Daniel Trad - June 9- 2000
*/


void radonwindow(float **data, int nh, int nt, float **model, int nq, int rtmethod,float *q, float *h, char *solver)
{
  
  
  //radontd_sparse(t,q,h,m,d,nt,nh,nq,dt,vel,dperv,pervmin,t0,inv_par inv,centralq,dataprec,nw,
  //		 fpeak, typewav, LI, nreg);  

  plotgather(data[0],nt,nh,"data_scale");

  return;

}






















