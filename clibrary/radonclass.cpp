#include <iostream.h>
#include "cwp.h"
#include "radonclass.hpp"
#include "segy.h"

void TGather::plot(char *name, char *key)
{
  char buf[80];
  save_gather0(d,nx,x,nt,dt,name);
  sprintf(buf,"suxwigb < %s key=%s title=%s &",name,key,name);
  system(buf);    
}
  
void TGather::gettaxis(float *to){
  memcpy((void *) to,(const void *) t,nt*sizeof(float));
}

void TGather::getxaxis(float *xo){
  memcpy((void *) xo,(const void *) x,nx*sizeof(float));
}


void TGather::getdata(float **data){
  int ix;
  for (ix=0;ix<nx;ix++)
    memcpy((void *) data[ix],(const void *) d[ix],nt*sizeof(float));
}

void TGather::setdata(float **data)
{
  int ix;
  for (ix=0;ix<nx;ix++)
    memcpy((void *) d[ix],(const void *) data[ix],nt*sizeof(float));
}

void TGather::gettrace(float *data, int ix){ 
  int it;
  for (it=0;it<nt;it++)
    memcpy((void *) data,(const void *) d[ix],nt*sizeof(float));
}

void TGather::settrace(float *data, int ix){
  int it;
  for (it=0;it<nt;it++)
    memcpy((void *) d[ix],(const void *) data,nt*sizeof(float));
}

TGather::TGather(int n1, int n2, float dt, float dx, float t0, float x0)
{
  int i;
  nx=n2;
  nt=n1;
  d=ealloc2float(n1,n2);
  memset( (void *) d[0], (int) '\0', n2 * n1 *FSIZE);
  t=ealloc1float(n1);
  x=ealloc1float(n2);
  t[0]=t0; for (i=1;i<n1;i++) t[i]=t[i-1]+dt;
  x[0]=x0; for (i=1;i<n2;i++) x[i]=x[i-1]+dx;
}


TGather::~TGather()
{
  free2float(d);
  free1float(t);
  free1float(x);
}



void TGather::save_gather0(float **d, int nh, float *h, int nt, float dt, char* name)
{
  
  segy tr;
  //tr=cleansegy(tr);
  int type=1; // Default = t - offset gather
              // option  type=2 : t - q Radon gather
  int itr;
  FILE* fp;
  if (fabs(h[1]-h[0]) < 0.1 ) type=2;

  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  //fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  for (itr=0;itr<nh;itr++){
      memcpy((void *) tr.data,
	     (const void *) d[itr],nt*sizeof(float));
      if (type==1) tr.offset=(int) h[itr];
      else if (type==2) tr.f2=h[itr];

      tr.tracl=itr;
      tr.dt=(int) (dt*1e6);
      tr.ns=(int) nt;
      tr.ntr=(int) nh;
      tr.delrt=0;
      tr.trid=1;
      fputtr(fp,&tr);
  }    

  fflush(fp);
  fclose(fp);

  return;
  
}




