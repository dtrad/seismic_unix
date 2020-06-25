#include "su.h"
#include "stdio.h"
#include "segy.h"

void savesudata(float **d, int nh, int nt, float dt, char *name, int shape)
{
  /* This functions writes a su file given a two dimesnional matrix data
     If data[nt][nh] shape=1
     if data[nh][nt] shape=2
  */
  
    FILE *filep;
    segy tr;
    int ih,it;
    if((filep=fopen(name,"w"))==NULL)
      err("cannot open file=%s\n",name);        
    for (ih=0;ih<nh;ih++){
      tr.tracl=ih;
      tr.ntr=nh;
      tr.dt=(unsigned short) dt*1e6;
      tr.ns=(unsigned short) nt;
      if (shape==1) for (it=0;it<nt;it++) tr.data[it]=d[it][ih];
      else if (shape==2) for (it=0;it<nt;it++) tr.data[it]=d[ih][it];
      fputtr(filep,&tr);
    }
    fclose(filep);
    return;
}

void savesudata(float *d, int nh, int nt, float dt, char *name)
{
  /* This functions writes a su file given a one dimensional matrix of 2D data
     data[ih*nt+it] 
  */
  
    FILE *filep;
    segy tr;
    int ih,it;
    if((filep=fopen(name,"w"))==NULL)
      err("cannot open file=%s\n",name);        
    for (ih=0;ih<nh;ih++){
      tr.tracl=ih;
      tr.ntr=nh;
      tr.dt=(unsigned short) dt*1e6;
      tr.ns=(unsigned short) nt;
      for (it=0;it<nt;it++) tr.data[it]=d[ih*nt+it];
      fputtr(filep,&tr);
    }
    fclose(filep);
    return;
}









