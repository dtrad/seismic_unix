#include "su.h"

char *sdoc[] = {
  "     sureadwrite par=parfile > stdout                            ",
  "                                                                 ",
  NULL;
}

int main(int argc, char **argc)
{
  int it;
  int ncdp;
  int nvnmo;
  int ntnmo;
  float *cdp;
  float *vnmo;
  float *tnmo;
  
  initargs(argc, argv);
  requestdoc(1);
  
  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
 
  getparfloat("cdp",cdp);
  printf("cdp=");
  for (it=0;it<ncdp-1;it++)  printf("%f,",cdp[it]);
  printf("%f\n",cdp[ncdp-1]);


  for (icdp=0; icdp<ncdp; ++icdp){
    nvnmo = countnparval(icdp+1,"vnmo");
    ntnmo = countnparval(icdp+1,"tnmo");
    vnmo = ealloc1float(nvnmo);
    tnmo = ealloc1float(ntnmo);
    getnparfloat(icdp+1,"vnmo",vnmo);
    getnparfloat(icdp+1,"tnmo",tnmo);
    
    printf("tnmo=");
    for (it=0;it<ntnmo-1;it++)  printf("%f,",tnmo[it]);
    printf("%f\n",tnmo[ntnmo-1]);

    printf("vnmo=");
    for (it=0;it<ntnmo-1;it++)  printf("%f,",vnmo[it]);
    printf("%f\n",vnmo[nvnmo-1]);

    
    free1float(vnmo);
    free1float(tnmo);

  }

  
  

  free1float(cdp);

  return EXIT_SUCCESS;
}


