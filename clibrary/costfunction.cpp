/*
  It computes the cost function 
  Inqut 
          m: model
          nx: number of model traces
          norm: implemented 1 Huber, 0 Cauchy, else L2
          Wm 
  Output
          J: standard deviation of the model
  Daniel Trad- 18 April 2000. UBC- Canada
*/
#include "su.h"
#include "math.h"

int costfunction(float *m, int nx, float *Wm, float rho)

{ 
      int i;
      float J=rho;
      for (i=0;i<nx;i++) J+=((Wm[i]*m[i])*(Wm[i]*m[i]));
      fprintf(stderr,"+++++++++++J=%f\n",J);
      return(J);
}






















































