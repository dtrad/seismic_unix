#include "su.h"
#include "stdio.h"

void savevector(float *d, int n, char* name)
{
  int nw;
  FILE* fp;

  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  nw=fwrite(d,n,sizeof(float),fp);
  if (nw!=n) warn("writing only %d of %d\n",nw,n);
  fclose(fp);
  return;
  
}













