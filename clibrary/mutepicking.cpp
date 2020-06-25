#include "radonfk.h"
//#include <iostream.h>
#define plot 1
#define VERBOSE 1


void mutepicking(float **data, int nh,  int nt, float* mutet, float threshold){
  int ih, it;
  for (ih=0;ih<nh;ih++){
    for (it=0;it<nt;it++){
      if (fabs(data[ih][it]) > threshold) mutetime[ih]=it;
      
    }
  }


}







