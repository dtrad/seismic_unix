#include "su.h"
#include "stddef.h"

void radon_param(float fmax, float *x, int nh, float dx,
float qmin, float *pqmaxt, float *pqmax, float *pdq, int nq, int  rtmethod)
/* Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum, it is under research.
    rtmethod=1 LRT
    rtmethod=2 PRT 
    Daniel Trad- UBC- 16-2-99
*/
{  
   float dq=0, qmax=0, qmaxt=0;
   float xmin=x[0];
   float xmax=x[nh-1];
   if (rtmethod==2) { //PRT
	  dq= 1/(fmax*(xmax-xmin)*(xmax-xmin));
          dq=0.8*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(2*fmax*(fabs(xmax-xmin))*dx);
    }
   else if(rtmethod==1) { //LRT
	  dq= 1/(fmax*(xmax-xmin));
          dq=0.8*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(fmax*dx);  
    }
    *pdq=dq;
    *pqmax=qmax;
    *pqmaxt=qmaxt;  
    return;
}

void radon_param(float fmax, float *x, int nh, float dx,
float qmin, float *pqmaxt, float *pqmax, float *pdq, int nq, int  rtmethod,
float factor)
/*

   Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum, it is under research.
    rtmethod=1 LRT
    rtmethod=2 PRT 
    This version allows to define a different factor to undersample the q space
    For example when we want to go beyond the qmax defined by Nyquist
    Daniel Trad- UBC- 16-2-99

*/
{  
   float dq=0, qmax=0, qmaxt=0;
   float xmin=x[0];
   float xmax=x[nh-1];
   
   if (rtmethod==2 || rtmethod==5) { //PRT
	  dq= 1/(fmax*(xmax-xmin)*(xmax-xmin));
          dq=factor*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(2*fmax*(fabs(xmax-xmin))*dx);
    }
   else if(rtmethod==1 || rtmethod==3) { //LRT
	  dq= 1/(fmax*fabs(xmax-xmin));
          dq=factor*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(fmax*dx);  
    }
    *pdq=dq;
    *pqmax=qmax;
    *pqmaxt=qmaxt;  
    return;
}
	















