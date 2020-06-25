#include <stdlib.h>
#include <stdio.h>
int  stopcriteria(float *gcv, float *rho, float normdata, int nh, int nq, int iter, int stopc)
{
  int num;
  int stop=0;
  int i=iter;

  if ((stopc==0) && (i>0)){
    if (1.1*rho[i-1] < rho[i]){
      stop=1;
      if (0) fprintf(stderr,"stop criteria at iteration %d\n",i+1);
    }
    else stop=0;
    
    return(stop);
  }
  else if (stopc==1){
    num=(nq-i)*(nq-i);
    gcv[i]=(rho[i])/num;
    if (0) fprintf(stderr,"gcv[%d]=%e\n",i,gcv[i]);
    if (i>0){
      if ( gcv[i-1] < gcv[i] ) { 
	if (0) fprintf(stderr,"GCV Criteria, iteration %d\n",i+1);
	stop=1;
	return(stop);
      }
      else{
	stop=0;
	return(stop);
      }
    }
  }
  else{
    stop=0;
    return(stop);
  }
  
  // It should never get here 
  return(stop);

}
 








