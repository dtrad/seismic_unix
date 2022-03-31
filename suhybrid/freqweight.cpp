#include "su.h"

float freqweight(int j, float df, float f1, float f2)
/*******************************************************************
return weight for each frequency
******************************************************************
Function parameters:

int j		freq index
float df	freq increment
float f1	taper off freq
float f2	freq beyond which all components are zero
*******************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*******************************************************************/
{
	float w;
	float f=j*df;
        //fprintf(stderr,"f2=%f,f1=%f\n,f=%f\n",f1,f2,f);
	if(f<=f1) return (1.);
	if(f>=f2) return (0.);

	w = (f2-f)/(f2-f1);
	return (w);
}



void freqweighting(complex **d, int nfreq, int nh, float df, float fmin, float fmax)
{

  float wa;
  int ifmin=(int) (fmin/df);
  int ifreq, i;

  for (ifreq=0 ; ifreq<nfreq ; ifreq++){
    if (ifreq < ifmin) for (i=0;i<nh;i++) d[ifreq][i]*=0;
    else{
      wa=freqweight(ifreq,df,fmax-10,fmax);
      for (i=0 ; i<nh ; i++)  d[ifreq][i]*=wa;
    }
  }
  
} 

