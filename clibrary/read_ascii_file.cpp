#include "su.h"

int read_ascii_file(const char *name,float *x)
{
  int nn;
  int ix;
  int nx;
  FILE *fp;
  fp=efopen(name,"r");
  
  ix=0;
  do{
    nn=fscanf(fp,"%f",&x[ix]);
    ix++;
    if (0) fprintf(stderr,"ix=%d\n",ix);
  }while(nn==1);
  nx=ix-1;
  
  efclose(fp);
  
  return(nx);
}
