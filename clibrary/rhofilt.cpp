#include "su.h"
#include "clibrarytd.h"

void rhofilt(float *d, float *t, float *h, int nt, int nh, float dt)
{

  int lfilter=51;
  int np2=lfilter/2;	
  float *ptrace, *ptrace2;
  float *rho;
  int it, ih;

  if ((ptrace=alloc1float(nt))==NULL)
    err("ptrace could not be allocated\n");
  if ((ptrace2=alloc1float(nt))==NULL)
    err("ptrace2 could not be allocated\n");   
  if ((rho=ealloc1float(lfilter))==NULL)
    err("rho could not be allocated\n"); 
  rho_filter(lfilter,nt,dt,rho);
  fprintf(stderr,"RHO filter......................\n");
     
  for (int ih=0;ih<nh;ih++){ 
    //for (it=0;it<nt;it++)  ptrace[it]=d[ih*nt+it];
    memcpy(ptrace,&d[ih*nt],nt*sizeof(float));
    conv(nt,-np2,ptrace,lfilter,0,rho,nt,0,ptrace2);
    memcpy(&d[ih*nt],ptrace2,nt*sizeof(float));
    //for (it=0;it<nt;it++) d[ih*nt+it]=ptrace2[it];
  }    
  
  return;
  
}













