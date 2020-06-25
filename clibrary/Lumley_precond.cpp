#include "su.h"
#include "math.h"

void Lumley_precond(float *M,float *t, float *h, int nt, int nh, float a, float b, float c)
{
  unsigned short it,ih;
  unsigned int ihxnt; 
  fprintf(stderr,"Inside Lumley nt=%d,nh=%d,a=%f,b=%f,c=%f\n",nt,nh,a,b,c);

  for (ih=0;ih<nh;ih++){
    ihxnt=ih*nt;
    for(it=0;it<nt;it++){
      M[ihxnt+it]=1+c*(1+pow((fabs(h[ih]/1000.)),a))/(pow((1+t[it]),b));
      // M[ihxnt+it]=1./(1+sqrt(fabs(h[ih]/1000.)));
      //M[ihxnt+it]=(1+sqrt(fabs(h[ih]/1000.)))/(1+t[it])/sqrt(nh);
      //fprintf(stderr,"M[%d]=%f\n",ihxnt+it,M[ihxnt+it]);
    }
  }
  return;
}











