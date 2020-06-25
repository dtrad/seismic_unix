#include "su.h"
#include "clibrary.h"
#include <math.h>
void matrix(complex **L, float *h, float *q,int nh,int nq,float w, float dq, int rtmethod)
{

        register int ih;
	register int iq;  
        complex  arg;

        for (iq=0;iq<nq;iq++){
	  for (ih=0;ih<nh;ih++){
              arg.r=0;
              if ((rtmethod==1)||(rtmethod==3)) arg.i=w*h[ih]*q[iq];
              else if (rtmethod==2) arg.i=w*h[ih]*h[ih]*q[iq];
	      
	      L[ih][iq]=exp(-arg);///sqrt(nq*nh);
	  }
	}
	
        return;
}

