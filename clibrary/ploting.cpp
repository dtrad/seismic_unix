#include "radonfk1.h"

void xplotgather(float **d, int nh, int nt, float dt, char *s, int num, char *s2)
{
  char buf[180];
  save_gather(d,nh,nt,dt,s);
  //sprintf(buf,"suximage < %s title=%s%d curve=curve1 npair=5 %s\n",s,s,num,s2);
  
  sprintf(buf,"suximage < %s title=%s%d legend=1 wbox=500 %s \n",s,s,num,s2);
  fprintf(stderr,buf);  
  system(buf);
  return;
}

void xplotdiff(float **d1, float **d2, int nh, int nt, float dt, char *s,  char *s2)
{
  char buf[180];
  int ih,it;
  float **d3=0;
  d3= ealloc2float(nt,nh);
  for (ih=0;ih<nh;ih++) for(it=0;it<nt;it++) d3[ih][it]=d1[ih][it]-d2[ih][it];

  save_gather(d3,nh,nt,dt,s);
  sprintf(buf,"suximage < %s title=%s legend=1 wbox=500 %s \n",s,s,s2);
  system(buf);
  
  free2float(d3);

  return;
}

/* test for static variable used in fputtr */
void plot_after_stretch(float **d, int nh, int nt, int dt, char *s1, char *s2)
{
  xplotgather(d,nh,nt,dt,s1,0,s2);
  return;
}

void ximageplotgather(float **d, int nh, int nt, float dt, char *s, int num, char *s2)
{
  char buf[180];
  save_gather(d,nh,nt,dt,s);
  sprintf(buf,"suximage < %s title=%s%d curve=curve1 npair=5 hbox=900 wbox=700 %s legend=1 \n",s,s,num,s2);
  system(buf);
  return;
}



