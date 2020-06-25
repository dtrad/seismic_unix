#include "su.h"
#include "cwp.h"

char *sdoc[] = {
  "     sureadwrite par=parfile > stdout                            ",
  "                                                                 ",
  "    Use this program to perform changes in par files             ",
  NULL
};

int main(int argc, char **argv)
{
  int it;
  int icdp;
  int ncdp;
  int nvnmo;
  int ntnmo;
  float *cdp;
  float *vnmo;
  float *tnmo;
  float factor;
  float thorizon;

  initargs(argc, argv);
  requestdoc(1);
  if (!getparfloat("thorizon",&thorizon)) thorizon = 100;
  if (!getparfloat("factor",&factor)) factor = 1;
  
  
  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
 
  getparfloat("cdp",cdp);
  printf("cdp=");
  for (it=0;it<ncdp-1;it++)  printf("%4.0f,",cdp[it]);
  printf("%4.0f\n",cdp[ncdp-1]);


  for (icdp=0; icdp<ncdp; icdp++){
    nvnmo = countnparval(icdp+1,"vnmo");
    ntnmo = countnparval(icdp+1,"tnmo");
    vnmo = ealloc1float(nvnmo);
    tnmo = ealloc1float(ntnmo);
    getnparfloat(icdp+1,"vnmo",vnmo);
    getnparfloat(icdp+1,"tnmo",tnmo);
    /* Apply changes to the velocity file to reflect a known velocity characteristic */
    for (it=0;it<nvnmo;it++){
      fprintf(stderr,"tnmo[%d]=%f,vnmo[%d]=%f\n",it,tnmo[it],it,vnmo[it]);
      if (tnmo[it]>thorizon){
	vnmo[it]=vnmo[it]*(1-factor*sin(it*PI/(2*nvnmo)));
	fprintf(stderr,"----------------->vnmo[%d]=%f\n",it,vnmo[it]);
      }
    }    
    printf("tnmo=");
    for (it=0;it<ntnmo-1;it++)  printf("%5.3f,",tnmo[it]);
    printf("%5.3f\n",tnmo[ntnmo-1]);

    printf("vnmo=");
    for (it=0;it<ntnmo-1;it++)  printf("%5.0f,",vnmo[it]);
    printf("%5.0f\n",vnmo[nvnmo-1]);

    
    free1float(vnmo);
    free1float(tnmo);

  }

  
  

  free1float(cdp);

  return EXIT_SUCCESS;
}


