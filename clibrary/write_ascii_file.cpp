#include "radonfk.h"

int write_curve_mute(mutemask_par par)
{
  
  FILE *fp;
  fp=efopen("curve1","w");
  
  fprintf(fp,"%f %f\n",par.tmin, par.ihmin);
  fprintf(fp,"%f %f\n",par.tmin, par.ihmax);
  fprintf(fp,"%f %f\n",par.tmax, par.ihmax);
  fprintf(fp,"%f %f\n",par.tmax, par.ihmin);
  fprintf(fp,"%f %f\n",par.tmin, par.ihmin);

  efclose(fp);
  
  return(nx);
}
